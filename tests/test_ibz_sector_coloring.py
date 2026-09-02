from alterseek.plotting_common import IBZ_UP_EXTRA_SECTOR_COLORS
from alterseek.symmetry import (
    _doubled_ibz_extra_flags,
    _hp1_four_sector_label_groups,
)


def test_ti_copied_sector_labels_include_zero_subscript_variants():
    ti1_labels = ["Gamma", "M", "X", "P", "Z", "Z_0", "N",
                  "X_A", "P_A"]
    ti2_labels = ["Gamma", "M", "X", "P", "N", "S_0", "S", "R", "G",
                  "N_A", "S_0A", "S_A", "R_A"]

    assert [label for label, is_extra in zip(
        ti1_labels, _doubled_ibz_extra_flags(ti1_labels)
    ) if is_extra] == ["X_A", "P_A"]
    assert [label for label, is_extra in zip(
        ti2_labels, _doubled_ibz_extra_flags(ti2_labels)
    ) if is_extra] == ["N_A", "S_0A", "S_A", "R_A"]


def test_ordinary_hpkot_subscript_two_labels_do_not_trigger_copied_sector():
    labels = ["Gamma", "L", "L_2", "S", "S_2"]

    assert not any(_doubled_ibz_extra_flags(labels))


def test_hp1_minus3_has_four_fixed_red_sector_groups():
    labels = [
        "\u0393", "A", "K", "H", "M", "L",
        "K_A", "H_A", "M_A", "L_A", "M_B", "L_B",
    ]
    groups = _hp1_four_sector_label_groups(labels)

    assert len(groups) == 4
    assert [[labels[index] for index in group] for group in groups] == [
        ["\u0393", "A", "K", "H", "M", "L"],
        ["\u0393", "A", "K_A", "H_A", "M_A", "L_A"],
        ["\u0393", "A", "K", "H", "M_B", "L_B"],
        ["\u0393", "A", "K_A", "H_A", "M", "L"],
    ]
    assert IBZ_UP_EXTRA_SECTOR_COLORS == (
        "#8b0000", "#d6604d", "#f4a582")
    top_view_colors = ("#b22222",) + IBZ_UP_EXTRA_SECTOR_COLORS
    assert top_view_colors[0] == "#b22222"
    assert len(set(top_view_colors)) == 4
