import numpy as np

from alterseek.compute_centroid_hybrid import _mC2d_centered_rectangular_path


def test_mc2_2d_path_uses_centered_rectangular_vertices():
    b_matrix = np.array([
        [1.21778639, 0.71215544, 0.0],
        [-1.21778639, 0.71215544, 0.0],
        [0.0, 0.0, 0.27956152],
    ])
    direct_lattice = 2 * np.pi * np.linalg.inv(b_matrix).T
    x_2d, points = _mC2d_centered_rectangular_path(direct_lattice, 2)

    assert np.isclose(x_2d, 0.33549631)
    assert points['Y'] == [0.5, 0.5, 0.0]
    assert points['S'] == [0.0, 0.5, 0.0]
    assert np.allclose(
        (np.array(points['C']) + np.array(points['SIGMA'])) / 2,
        points['S'],
    )
