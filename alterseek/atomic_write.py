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

# Scratch files stay beside their targets so os.replace is atomic, but their private names do not repeat the possibly long public basename.
_SCRATCH_TOKEN_HEX = 16


def _scratch_token():
    return uuid.uuid4().hex[:_SCRATCH_TOKEN_HEX]


def _atomic_write_text(path, text):
    """Write UTF-8 text beside *path*, then atomically replace the target."""
    target = os.path.abspath(os.fspath(path))
    parent = os.path.dirname(target)
    # O_CREAT|O_EXCL with mode 0666 (instead of tempfile.mkstemp) lets the kernel apply the process umask itself, so a brand-new target gets the same mode a plain open() would create without mkstemp's 0600 or a process-global os.umask() round-trip that could race other threads.
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
        # Keep an existing target's mode (for example, group-readable outputs in shared cluster directories); a new target keeps the umask default.
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
