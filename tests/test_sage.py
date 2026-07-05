import numpy as np
from numpy.linalg import matrix_rank
import warnings
from guptri_py import guptri, kcf_blocks


def _check_guptri_properties(A, B):
    try:
        from sage.all import matrix, RDF
    except ImportError:
        warnings.warn(UserWarning("Sage is not installed, so sage tests are skipped"))
        return
    S, T, P, Q, kstr = guptri(A, B)
    tol = 1e-12
    assert (P.H * A * Q - S).norm() < tol
    assert (P.H * B * Q - T).norm() < tol
    assert (P.H * P - matrix.identity(RDF, P.ncols())).norm() < tol
    assert (Q.H * Q - matrix.identity(RDF, Q.ncols())).norm() < tol
    kb = kcf_blocks(kstr)
    assert np.all(kb[:, 0] == 0) or kb[0, 0] < kb[1, 0]
    assert np.all(kb[:, -1] == 0) or kb[0, -1] > kb[1, -1]
    assert np.all(kb[0, 1:4] == kb[1, 1:4])

    # test that Y = A X + B X for some of the reducing subspaces
    for k in range(1, 5):
        Y, X = guptri(A, B, part=range(k))[2:4]
        AXBX = (A * X).augment(B * X)
        assert Y.ncols() == matrix_rank(AXBX, tol=1e-12)
        assert Y.ncols() == matrix_rank(AXBX.augment(Y), tol=1e-12)


def test_1():
    try:
        from sage.all import matrix, RDF
    except ImportError:
        warnings.warn(UserWarning("Sage is not installed, so sage tests are skipped"))
        return
    A = matrix(RDF, [[0, 1, 0], [0, 0, 2]])
    B = matrix(RDF, [[0, 0, 0], [0, 0, 3]])
    _check_guptri_properties(A, B)


def test_2():
    try:
        from sage.all import matrix
    except ImportError:
        warnings.warn(UserWarning("Sage is not installed, so sage tests are skipped"))
        return
    A = np.array([[22, 34, 31, 31, 17],
                  [45, 45, 42, 19, 29],
                  [39, 47, 49, 26, 34],
                  [27, 31, 26, 21, 15],
                  [38, 44, 44, 24, 30]], float)
    B = np.array([[13, 26, 25, 17, 24],
                  [31, 46, 40, 26, 37],
                  [26, 40, 19, 25, 25],
                  [16, 25, 27, 14, 23],
                  [24, 35, 18, 21, 22]], float)
    _check_guptri_properties(matrix(A), matrix(B))


def test_3():
    try:
        from sage.all import matrix, CDF
    except ImportError:
        warnings.warn(UserWarning("Sage is not installed, so sage tests are skipped"))
        return
    A = matrix(CDF, [[1+1j, 3e-16j], [2e-16j, 0]])
    B = matrix(CDF, [[1, 1e-16j], [1e-16, 0]])
    _check_guptri_properties(A, B)


def test_4():
    try:
        from sage.doctest.control import DocTestDefaults, DocTestController
    except ImportError:
        warnings.warn(UserWarning("Sage is not installed, so sage tests are skipped"))
        return
    import guptri_py
    # "PendingDeprecationWarning: the matrix subclass is not the recommended way to represent matrices or deal with linear algebra"
    warnings.filterwarnings('ignore', category=PendingDeprecationWarning)  # TODO caused by sage's conversion of matrix to numpy
    dd = DocTestDefaults()
    dc = DocTestController(dd, guptri_py.__path__)
    error_status = dc.run()
    assert error_status == 0, "Doctests failed"
