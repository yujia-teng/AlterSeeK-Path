import matplotlib.pyplot as plt
import numpy as np
import pytest

from plot_alterband import _draw_panel, _validate_plot_config, plot_alterband
from plot_alterband_qe import _validate_plot_config as _validate_qe_plot_config


def test_spin_up_band_is_drawn_above_spin_down_band():
    fig, ax = plt.subplots()
    try:
        _draw_panel(
            ax,
            labels=["GAMMA", "M"],
            positions=[0.0, 1.0],
            tick_labels=["Gamma", "M"],
            kpath=np.array([0.0, 0.5, 1.0]),
            bands_up=np.array([[0.0], [0.1], [0.0]]),
            bands_dw=np.array([[0.0], [0.1], [0.0]]),
            elim=(-1.0, 1.0),
            xlim=(0.0, 1.0),
            gap_half=0.01,
            font_size=12,
            rotate_xtick_labels=False,
            xtick_rotation=45.0,
        )

        spin_down = next(line for line in ax.lines if line.get_zorder() == 2)
        spin_up = next(line for line in ax.lines if line.get_zorder() == 3)
        assert spin_down.get_color() == "red"
        assert spin_down.get_linewidth() == 0.7
        assert spin_up.get_color() == "black"
        assert spin_up.get_linewidth() == 0.7
        assert spin_up.get_zorder() > spin_down.get_zorder()
    finally:
        plt.close(fig)


def test_vaspkit_label_repair_is_in_memory_and_reported(tmp_path, capsys):
    klabels = tmp_path / "KLABELS"
    original = "A 0.0\nC 0.5\nD 1.0\n"
    klabels.write_text(original, encoding="utf-8")
    (tmp_path / "KPOINTS").write_text(
        "Path\n30\nLine-Mode\nReciprocal\n"
        "0.0 0.0 0.0 A\n0.5 0.0 0.0 B\n\n"
        "0.5 0.25 0.0 C\n1.0 0.25 0.0 D\n",
        encoding="utf-8",
    )
    band_text = "k band\n0.0 0.0\n0.5 0.1\n1.0 0.0\n"
    up = tmp_path / "up.dat"
    down = tmp_path / "down.dat"
    up.write_text(band_text, encoding="utf-8")
    down.write_text(band_text, encoding="utf-8")

    plot_alterband(
        klabels=klabels,
        band_up=up,
        band_down=down,
        output=tmp_path / "bands.png",
    )

    assert klabels.read_text(encoding="utf-8") == original
    output = capsys.readouterr().out
    assert "B|C" in output
    assert "KLABELS was not modified" in output


@pytest.mark.parametrize(
    "validator, config",
    [
        (_validate_plot_config, {"emin_typo": -2}),
        (_validate_plot_config, {"rotate_xtick_labels": "yes"}),
        (_validate_qe_plot_config, {"lattice_type": "hP2"}),
        (_validate_qe_plot_config, {"split_panels": 4}),
        (_validate_plot_config, {"split_panels": 0}),
        (_validate_plot_config, {"output": None}),
        (_validate_qe_plot_config, {"fermi_ev": None}),
    ],
)
def test_plot_config_rejects_unknown_or_invalid_settings(validator, config, tmp_path):
    with pytest.raises(ValueError):
        validator(config, tmp_path / "plot.toml")


@pytest.mark.parametrize("module_name", ["plot_alterband", "plot_alterband_qe"])
def test_main_reports_bad_config_without_traceback(module_name, tmp_path):
    import importlib

    main = importlib.import_module(module_name).main
    config = tmp_path / "plot.toml"
    config.write_text("emin_typo = 1\n", encoding="utf-8")
    with pytest.raises(SystemExit) as excinfo:
        main(["--config", str(config)])
    assert "Unknown setting" in str(excinfo.value)


def test_main_reports_domain_invalid_setting_cleanly(tmp_path):
    from plot_alterband import main as vasp_main

    klabels = tmp_path / "KLABELS"
    klabels.write_text("A 0.0\nB 1.0\n", encoding="utf-8")
    (tmp_path / "KPOINTS").write_text(
        "Path\n30\nLine-Mode\nReciprocal\n"
        "0.0 0.0 0.0 A\n1.0 0.0 0.0 B\n",
        encoding="utf-8",
    )
    band_text = "k band\n0.0 0.0\n0.5 0.1\n1.0 0.0\n"
    up = tmp_path / "up.dat"
    down = tmp_path / "down.dat"
    up.write_text(band_text, encoding="utf-8")
    down.write_text(band_text, encoding="utf-8")
    config = tmp_path / "alterband.toml"
    config.write_text(
        f'klabels = "{klabels.as_posix()}"\n'
        f'up = "{up.as_posix()}"\n'
        f'down = "{down.as_posix()}"\n'
        "gap_frac = -1\n",
        encoding="utf-8",
    )
    with pytest.raises(SystemExit) as excinfo:
        vasp_main(["--config", str(config)])
    assert "gap_frac" in str(excinfo.value)


@pytest.mark.parametrize(
    "first, second, merged",
    [("k", "k'", "k|k'"), ("k'", "k", "k'|k")],
)
def test_qe_build_tick_data_merges_gap_pairs_in_either_order(first, second, merged):
    from plot_alterband_qe import _build_tick_data

    kpath = np.array([0.0, 1.0, 3.0, 4.0])
    waypoints = [("GAMMA", 30, 0), (first, 30, 1), (second, 1, 2), ("M", 30, 3)]

    labels, positions = _build_tick_data(waypoints, kpath)

    assert labels == ["GAMMA", merged, "M"]
    assert positions == [0.0, 2.0, 4.0]


def test_qe_build_tick_data_keeps_unpaired_helper_labels():
    from plot_alterband_qe import _build_tick_data

    kpath = np.array([0.0, 1.0, 2.0])
    waypoints = [("k", 1, 0), ("M", 30, 1), ("k'", 30, 2)]

    labels, positions = _build_tick_data(waypoints, kpath)

    assert labels == ["k", "M", "k'"]
    assert positions == [0.0, 1.0, 2.0]


def test_qe_parse_kpoints_truncated_after_card_raises_value_error(tmp_path):
    from plot_alterband_qe import _parse_kpoints_qe

    bad = tmp_path / "KPOINTS_alter_qe"
    bad.write_text("K_POINTS crystal_b\n", encoding="utf-8")
    with pytest.raises(ValueError, match="truncated"):
        _parse_kpoints_qe(bad)


def test_qe_parse_kpoints_rejects_missing_waypoint_rows(tmp_path):
    from plot_alterband_qe import _parse_kpoints_qe

    bad = tmp_path / "KPOINTS_alter_qe"
    bad.write_text(
        "K_POINTS crystal_b\n"
        "2\n"
        "0.0 0.0 0.0 30 ! GAMMA\n",
        encoding="utf-8",
    )
    with pytest.raises(
        ValueError,
        match="declares 2 waypoints but contains only 1",
    ):
        _parse_kpoints_qe(bad)


def test_qe_parse_kpoints_malformed_waypoint_raises_value_error(tmp_path):
    from plot_alterband_qe import _parse_kpoints_qe

    bad = tmp_path / "KPOINTS_alter_qe"
    bad.write_text(
        "K_POINTS crystal_b\n2\n0.0 0.0 0.0 30 !GAMMA\n0.5 0.5 !X\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Malformed waypoint line"):
        _parse_kpoints_qe(bad)
