
======
`pcsr`
======

.. Start of sphinx doc include.
.. start long description.
.. start badges.

.. end badges.

Point Cloud Surface Reconstruction, python interface for reconstructing a (triangle)
mesh from a point-cloud.


.. end long description.


Installation
============

Install from latest github source:

.. code-block:: console

   $ python -m pip install --user git+git://github.com/Shane-J-Latham/pcsr.git#egg=pcsr

or from source directory:

.. code-block:: console

   $ python -m pip install --prefix=/path/to/install/root --no-deps --no-build-isolation .


Requirements
============

Requires:

- python-3 version `>= 3.4`,
- `boost <https://boost.org>`_ version `>= 1.65`
- `numpy <https://www.numpy.org/>`_ version `>= 1.14.5`,
- `CGAL <https://cgal.org/>`_ version `>= 3.5`,
- `PCL <https://pointclouds.org/>`_ version `>= 1.8.1`,

Testing
=======

Run tests (unit-tests) using::

   python -m pcsr.tests


Latest source code
==================

Source at github:

   https://github.com/Shane-J-Latham/pcsr


License information
===================

See the file `LICENSE <https://github.com/Shane-J-Latham/pcsr/blob/dev/LICENSE>`_
for terms & conditions, for usage and a DISCLAIMER OF ALL WARRANTIES.

