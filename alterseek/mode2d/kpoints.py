"""KPOINTS workflow entry for two-dimensional mode."""

from ..kpoints import KPointsModifier


class KPointsModifier2D(KPointsModifier):
    """Run the shared interaction and writers with 2D analysis enabled."""

    def __init__(self, input_vacuum_axis=None):
        super().__init__(
            mode_2d=True,
            input_vacuum_axis=input_vacuum_axis,
        )


def run_2d(input_vacuum_axis=None):
    """Run the interactive 2D workflow."""
    return KPointsModifier2D(
        input_vacuum_axis=input_vacuum_axis
    ).interactive_modify()
