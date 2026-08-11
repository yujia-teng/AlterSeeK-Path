"""Atomic text writes.

Every generated file goes through here, so a crash or a full disk mid-write
leaves the previous version intact rather than a half-written POSCAR, MCIF, or
operation list that looks loadable.

Leaf module: stdlib only, so it can be imported from anywhere in the package
(`find_sf_operations` in particular cannot import `io`, which imports it).
"""
import contextlib
import os
import uuid
from io import StringIO

# Scratch files stay beside their targets so os.replace is atomic, but their
# private names do not repeat the possibly long public basename.
_SCRATCH_TOKEN_HEX = 16


def _scratch_token():
    return uuid.uuid4().hex[:_SCRATCH_TOKEN_HEX]


def _atomic_write_text(path, text):
    """Write UTF-8 text beside *path*, then atomically replace the target."""
    target = os.path.abspath(os.fspath(path))
    parent = os.path.dirname(target)
    # O_CREAT|O_EXCL with mode 0666 (instead of tempfile.mkstemp) lets the
    # kernel apply the process umask itself, so a brand-new target gets the
    # same mode a plain open() would create — without mkstemp's 0600 or a
    # process-global os.umask() round-trip that could race other threads.
    temporary = os.path.join(
        parent, f".alterseek-{_scratch_token()}.tmp"
    )
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
    fd = os.open(temporary, flags, 0o666)
    try:
        handle = os.fdopen(fd, "w", encoding="utf-8", newline="\n")
    except Exception:
        # os.fdopen never took ownership, so the raw descriptor is still ours.
        os.close(fd)
        try:
            os.remove(temporary)
        except OSError:
            pass
        raise
    try:
        with handle:
            handle.write(text)
        # Keep an existing target's mode (e.g. group-readable outputs in
        # shared cluster directories); a new target keeps the umask default.
        try:
            os.chmod(temporary, os.stat(target).st_mode & 0o7777)
        except OSError:
            pass
        os.replace(temporary, target)
    except Exception:
        try:
            os.remove(temporary)
        except OSError:
            pass
        raise


def _atomic_write_text_set(outputs):
    """Commit several UTF-8 text files as one rollback-protected set.

    All new contents are staged beside their targets before any public file is
    changed. Existing targets are then moved to private backups and the staged
    files are promoted in insertion order. If an ordinary filesystem exception
    interrupts either phase, every target is restored to its original bytes (or
    removed when it did not previously exist) before the exception is reraised.

    This protects caught failures; it is not crash- or power-loss atomicity.
    """
    entries = []
    seen_targets = set()
    transaction_id = _scratch_token()
    preserve_failed_backups = False

    try:
        for index, (path, text) in enumerate(outputs.items()):
            target = os.path.abspath(os.fspath(path))
            if target in seen_targets:
                raise ValueError(f"duplicate transaction target: {target}")
            seen_targets.add(target)
            parent = os.path.dirname(target)
            staged = os.path.join(
                parent, f".alterseek-{transaction_id}-{index}.stage"
            )
            backup = os.path.join(
                parent, f".alterseek-{transaction_id}-{index}.backup"
            )
            if os.path.exists(staged) or os.path.exists(backup):
                raise FileExistsError(
                    f"scratch path already exists for transaction target {target}"
                )
            entry = {
                "target": target,
                "staged": staged,
                "backup": backup,
                "existed": os.path.exists(target),
                "backed_up": False,
                "promoted": False,
                "staged_owned": False,
            }
            entries.append(entry)
            flags = (
                os.O_WRONLY | os.O_CREAT | os.O_EXCL
                | getattr(os, "O_BINARY", 0)
            )
            fd = os.open(staged, flags, 0o666)
            entry["staged_owned"] = True
            try:
                handle = os.fdopen(fd, "w", encoding="utf-8", newline="\n")
            except Exception:
                os.close(fd)
                raise
            with handle:
                handle.write(text)
            try:
                os.chmod(staged, os.stat(target).st_mode & 0o7777)
            except OSError:
                pass

        for entry in entries:
            if entry["existed"]:
                os.replace(entry["target"], entry["backup"])
                entry["backed_up"] = True
        for entry in entries:
            os.replace(entry["staged"], entry["target"])
            entry["promoted"] = True
    except Exception as exc:
        rollback_errors = []
        for entry in reversed(entries):
            try:
                if entry["backed_up"]:
                    os.replace(entry["backup"], entry["target"])
                    entry["backed_up"] = False
                elif entry["promoted"] and os.path.exists(entry["target"]):
                    os.remove(entry["target"])
            except OSError as rollback_exc:
                rollback_errors.append(str(rollback_exc))
        if rollback_errors:
            preserve_failed_backups = True
            raise RuntimeError(
                "multi-file output failed and rollback was incomplete: "
                + "; ".join(rollback_errors)
            ) from exc
        raise
    finally:
        for entry in entries:
            if entry["staged_owned"]:
                try:
                    os.remove(entry["staged"])
                except OSError:
                    pass
            if entry["backed_up"] and not preserve_failed_backups:
                try:
                    os.remove(entry["backup"])
                except OSError:
                    pass


@contextlib.contextmanager
def _atomic_open_text(path):
    """Collect `.write()` calls and commit them atomically on a clean exit.

    A drop-in replacement for `open(path, "w", encoding="utf-8",
    newline="\\n")` in the writers that build a file line by line. If the body
    raises partway through, nothing is written at all -- previously those
    writers left a truncated file behind.
    """
    buffer = StringIO()
    yield buffer
    _atomic_write_text(path, buffer.getvalue())
