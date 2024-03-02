"""
Instances of meshes.

.. rubric:: Models

.. autosummary::
   :toctree: generated/

   cow - A cow mesh model.
   unit_cube - A unit-cube mesh model.


"""

from __future__ import absolute_import
import os as _os
import sys as _sys
from importlib import resources as _resources


class ModelLoader(object):
    """
    Callable which loads a specified mesh model from file.
    """

    def __init__(self, file_name, file_type=None):
        """
        """
        if file_type is None:
            file_type = _os.path.splitext(file_name)[1][1:].tolower()
        self._file_name = file_name
        self._file_type = file_type

    @property
    def model_dir(self):
        """
        Name of directory from which file is loaded.
        """
        return "models/data"

    @property
    def file_name(self):
        """
        Name of file from which mesh is loaded.
        """
        return self._file_name

    @property
    def file_type(self):
        """
        A :obj:`str` indicating the type of file format.
        """
        return self._file_type

    def __call__(self):
        """
        Returns the mesh loaded from the :attr:`file_name` file.
        """
        from trimesh import load_mesh

        if _sys.version_info >= (3, 9):
            with \
                _resources.files("pcsr").joinpath(
                    self.model_dir + "/" + self.file_name
                ).open("rb") as fp:  # noqa: E125

                mesh = load_mesh(fp, self.file_type)
        else:
            with \
                _resources.open_binary(
                    "pcsr.models.data",
                    self.file_name
                ) as fp:  # noqa: E125

                mesh = load_mesh(fp, self.file_type)

        return mesh

    def __repr__(self):
        """
        """
        return \
            (
                "ModelLoader(file_name=%s, file_type=%s)"
                %
                (
                    repr(self.file_name),
                    repr(self.file_type)
                )
            )

    def __str__(self):
        """
        """
        return self.__repr__()


#: Returns model of a cow.
cow = ModelLoader("cow.obj", "obj")

#: Returns model of a unit cube.
unit_cube = ModelLoader("unit_cube.obj", "obj")

__all__ = [s for s in dir() if not s.startswith('_')]
