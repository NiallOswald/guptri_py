import numpy as np

from guptri_py import guptri, kcf_blocks


def test_guptri():
    A = np.array(
        [
            [22, 34, 31, 31, 17],
            [45, 45, 42, 19, 29],
            [39, 47, 49, 26, 34],
            [27, 31, 26, 21, 15],
            [38, 44, 44, 24, 30],
        ],
        np.float64,
    )
    B = np.array(
        [
            [13, 26, 25, 17, 24],
            [31, 46, 40, 26, 37],
            [26, 40, 19, 25, 25],
            [16, 25, 27, 14, 23],
            [24, 35, 18, 21, 22],
        ],
        np.float64,
    )

    S, T, P, Q, kstr = guptri(A, B)

    assert np.allclose(A - P @ S @ Q.T.conj(), 0.0)
    assert np.allclose(B - P @ T @ Q.T.conj(), 0.0)

    kb = kcf_blocks(kstr)
    kb_true = np.array(
        [
            [0, 2, 1, 1, 1],
            [1, 2, 1, 1, 0],
        ]
    )
    assert np.allclose(kb, kb_true)


def test_kcf_blocks():
    A = np.array(
        [
            [0, 1, 0],
            [0, 0, 2],
        ],
        dtype=np.float64,
    )
    B = np.array(
        [
            [0, 0, 0],
            [0, 0, 3],
        ],
        dtype=np.float64,
    )

    _, _, _, _, kstr = guptri(A, B)

    kb = kcf_blocks(kstr)
    kb_true = np.array(
        [
            [0, 0, 1, 1, 0],
            [1, 0, 1, 1, 0],
        ]
    )

    assert np.allclose(kb, kb_true)
