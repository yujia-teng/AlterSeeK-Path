from symmetry import _doubled_ibz_extra_flags


def test_ti_copied_sector_labels_include_zero_subscript_variants():
    ti1_labels = ["Gamma", "M", "X", "P", "Z", "Z_0", "N",
                  "M_A", "N_A", "Z_0A"]
    ti2_labels = ["Gamma", "M", "X", "P", "N", "S_0", "S", "R", "G",
                  "N_A", "S_0A", "S_A", "R_A"]

    assert [label for label, is_extra in zip(
        ti1_labels, _doubled_ibz_extra_flags(ti1_labels)
    ) if is_extra] == ["M_A", "N_A", "Z_0A"]
    assert [label for label, is_extra in zip(
        ti2_labels, _doubled_ibz_extra_flags(ti2_labels)
    ) if is_extra] == ["N_A", "S_0A", "S_A", "R_A"]


def test_ordinary_hpkot_subscript_two_labels_do_not_trigger_copied_sector():
    labels = ["Gamma", "L", "L_2", "S", "S_2"]

    assert not any(_doubled_ibz_extra_flags(labels))
