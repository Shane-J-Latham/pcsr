from __future__ import absolute_import
import os as _os

__doc__ = \
    """
Instances of meshes.

.. rubric:: Models

.. autosummary::
   :toctree: generated/


"""


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
        return _os.path.split(__file__)[0]

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
        import trimesh

        return \
            trimesh.load_mesh(
                _os.path.join(self.model_dir, self.file_name),
                self.file_type
            )

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


def append_model_loader_attr_to_doc_string(doc_string):
    """
    Finds attributes of this module which are of instance :obj:`ModelLoader`
    and appends an entry to the :samp:`{doc_string}` string.

    :type doc_string: :obj:`str`
    :param doc_string: Doc string which gets :obj:`ModelLoader` attributes appended.
    :rtype: :obj:`str`
    :return: The :samp:`{doc_string}` with autosummary attributes appended.
    """
    import sys
    this_module = sys.modules[__name__]
    for attr_name in dir(this_module):
        if isinstance(getattr(this_module, attr_name), ModelLoader):
            doc_string = doc_string + ("   %s - A model.\n" % attr_name)
    return doc_string


#: Returns model of a cow.
cow = ModelLoader("cow.obj", "obj")

#: Returns model of a unit cube.
unit_cube = ModelLoader("unit_cube.obj", "obj")

__doc__ = append_model_loader_attr_to_doc_string(__doc__)

__all__ = [s for s in dir() if not s.startswith('_')]
