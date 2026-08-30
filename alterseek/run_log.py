"""Capture the terminal transcript of one successful interactive run."""
import contextlib
import datetime
import os
import re
import sys
from io import StringIO

from .atomic_write import _atomic_write_text


RUN_LOG_FILENAME = "alterseek_run.log"
# Terminal bold/colour codes are unreadable in a text file.
_ANSI = re.compile(r"\x1b\[[0-9;]*m")


class _RunLog:
    """Buffer one run's transcript until the workflow succeeds."""

    def __init__(self, path):
        self._path = path
        self._buffer = StringIO()
        self._successful = False

    def write(self, text):
        self._buffer.write(_ANSI.sub("", text))

    def flush(self):
        pass

    def mark_success(self):
        self._successful = True

    def save(self):
        parent = os.path.dirname(self._path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        footer = f"Run log: {self._path}\n"
        text = (
            f"# AlterSeeK-Path run log, {stamp}\n\n"
            f"{self._buffer.getvalue()}{footer}"
        )
        _atomic_write_text(self._path, text)


class _Tee:
    """Send one output stream to the terminal and to the run log at once."""

    def __init__(self, stream, log):
        self._stream = stream
        self._log = log

    def write(self, text):
        self._stream.write(text)
        self._log.write(text)
        return len(text)

    def flush(self):
        self._stream.flush()
        self._log.flush()

    def __getattr__(self, name):
        return getattr(self._stream, name)


_active_log = None


def read_answer():
    """Read one typed answer, keeping it in the run log."""
    answer = input()
    if _active_log is not None:
        _active_log.write(f"{answer}\n")
        _active_log.flush()
    return answer


@contextlib.contextmanager
def run_log(output_path):
    """Capture a run and save its transcript after it is marked successful."""
    global _active_log
    log = _RunLog(output_path)
    saved_stdout, saved_stderr, saved_log = sys.stdout, sys.stderr, _active_log
    sys.stdout = _Tee(saved_stdout, log)
    sys.stderr = _Tee(saved_stderr, log)
    _active_log = log
    try:
        yield log
    finally:
        sys.stdout, sys.stderr = saved_stdout, saved_stderr
        _active_log = saved_log
        if log._successful:
            try:
                log.save()
            except OSError as exc:
                print(f"[Warning] Could not write the run log: {exc}",
                      file=sys.stderr)
            else:
                print(f"Run log: {output_path}")
