"""Output location and vacuum-axis values shared across the package."""
import numbers

# Keep diagnostics under OUTPUT_DIR; calculation inputs and plot configs stay in the working directory.
OUTPUT_DIR = "alterseek_output"

_VACUUM_AXIS_INDEX = {"a": 0, "b": 1, "c": 2}


def _normalize_vacuum_axis(value, *, allow_none=False, allow_index=False):
    """Turn a vacuum-axis letter into its cell-vector index."""
    if value is None and allow_none:
        return None
    if isinstance(value, str):
        axis = value.strip().lower()
        if axis in _VACUUM_AXIS_INDEX:
            return _VACUUM_AXIS_INDEX[axis]
    elif (allow_index and not isinstance(value, bool)
            and isinstance(value, numbers.Integral) and int(value) in (0, 1, 2)):
        return int(value)
    if allow_index:
        raise ValueError('vacuum_axis must be "a", "b", or "c", or the index 0, 1, or 2')
    raise ValueError('vacuum_axis must be "a", "b", or "c"')
