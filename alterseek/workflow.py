"""The complete AlterSeeK-Path workflow as one callable entry point."""
import os

from .constants import OUTPUT_DIR, _normalize_vacuum_axis
from .kpoints import KPathBuilder
from .mode2d.kpoints import KPathBuilder2D
from .run_log import RUN_LOG_FILENAME, run_log


def run_workflow(mode_2d=False, vacuum_axis=None):
    """Run the interactive Step 0-5 workflow and save its run log on success.

    ``vacuum_axis`` is ``"a"``, ``"b"``, or ``"c"``; ``None`` leaves the choice
    to ``alterseek_input.toml``. Returns True on success and False after a
    reported failure.
    """
    axis = _normalize_vacuum_axis(vacuum_axis, allow_none=True)
    builder = (KPathBuilder2D(input_vacuum_axis=axis) if mode_2d
               else KPathBuilder(input_vacuum_axis=axis))
    with run_log(os.path.join(OUTPUT_DIR, RUN_LOG_FILENAME)) as log:
        success = builder.interactive_build()
        if success:
            log.mark_success()
    return success
