"""Atomic UTF-8 text writes that never expose partial output."""
import contextlib
import os
import uuid
from io import StringIO

_SCRATCH_TOKEN_HEX = 16


def _scratch_token():
    return uuid.uuid4().hex[:_SCRATCH_TOKEN_HEX]


def _atomic_write_text(path, text):
    """Write UTF-8 text beside *path*, then atomically replace the target."""
    target = os.path.abspath(os.fspath(path))
    parent = os.path.dirname(target)
    # Mode 0666 lets the process umask set permissions.
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
        # Preserve permissions when replacing an existing file.
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
    """Buffer text and replace *path* atomically after a clean exit."""
    buffer = StringIO()
    yield buffer
    _atomic_write_text(path, buffer.getvalue())
