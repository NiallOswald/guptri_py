#############
``guptri_py``
#############

.. image:: https://readthedocs.org/projects/guptri-py/badge/?version=latest
   :target: https://guptri-py.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

***************************************
A GUPTRI wrapper for NumPy and SageMath
***************************************

This Python package provides Python bindings for the software GUPTRI_ and
can be used with `NumPy <NUMPY_>`_ and, optionally, `SageMath <SAGE_>`_.

GUPTRI is a Fortran library by Jim Demmel and Bo Kågström for robust
computation of generalized eigenvalues of singular matrix pencils.
Standard tools like LAPACK do not reliably handle singular generalized
eigenvalue problems.

GUPTRI solves this by computing a generalized block upper triangular form
(generalized Schur staircase form) of a matrix pencil, revealing the Kronecker
structure of the pencil. For details, see the `documentation <guptri_py_rtd_>`_
and the references therein.

Examples
========

See the examples and documentation at
`https://guptri-py.readthedocs.io <guptri_py_rtd_>`_.

Installation
============

**Requirements**: `NumPy <NUMPY_>`_ and, optionally, `SageMath <SAGE_>`_.

To install with Python and NumPy, run the following command::

    pip install git+https://github.com/mwageringel/guptri_py

Alternatively, for use with Sage, run this command::

    sage -pip install git+https://github.com/mwageringel/guptri_py

Optionally, pass `--user` to install into the Python user install directory (no root access required).

After successful installation, to run the tests::

    git clone https://github.com/mwageringel/guptri_py.git && cd guptri_py
    (cd tests && python -m pytest .)

Installing into a virtual environment (with system packages)::

    python -m venv --system-site-packages ./venv   # assumes that numpy and meson (and optionally sage) are installed system-wide

    ./venv/bin/pip install --upgrade --no-index -v .
    # ./venv/bin/pip install --upgrade --no-index -v --no-build-isolation .   # this alternative may be needed to find meson

    (cd venv/ && ./bin/python -m pytest ../tests/)  # running tests (changing to a different directory is important)

    (cd ./venv && ./bin/python)   # it is important to change to a different directory
    (cd ./venv && ./bin/python -m IPython)   # (or using IPython)

Installing into a virtual environment (without system packages)::

    python -m venv ./venv
    ./venv/bin/pip install git+https://github.com/mwageringel/guptri_py
    ./venv/bin/pip install pytest
    (cd ./venv && ./bin/python -m pytest ../tests/)
    (cd ./venv && ./bin/python -m IPython)

Issues
------

* With NumPy ≤ 1.17, it may be necessary to set::

    export NPY_DISTUTILS_APPEND_FLAGS=1

  to fix a linking problem. See https://github.com/numpy/numpy/issues/12799.

.. _SAGE: https://www.sagemath.org/
.. _GUPTRI: https://web.archive.org/web/20080920172251/https://www8.cs.umu.se/research/nla/singular_pairs/guptri/
.. _NUMPY: https://numpy.org/
.. _guptri_py_gh: https://github.com/mwageringel/guptri_py
.. _guptri_py_rtd: https://guptri-py.readthedocs.io/en/latest/#module-guptri_py
