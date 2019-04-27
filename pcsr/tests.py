"""
Tests for :mod:`pcsr`.
"""
from __future__ import absolute_import
import numpy as _np
import unittest as _unittest

have_trimesh = False
try:
    import trimesh as _trimesh
    have_trimesh = True
except Exception:
    pass

class SurfaceReconTest(_unittest.TestCase):

    """
    Base class for surface reconstruction.
    """

    def setUp(self):
        """
        """
        import copy
        from . import models

        _np.random.seed(54317953)

        if have_trimesh:
            self.cow_mesh = copy.deepcopy(models.cow())
            # scale the mesh to fit inside the image volume.
            cow_radius = self.cow_mesh.bounding_sphere.primitive.radius
            self.cow_mesh.apply_translation(-_np.asarray(self.cow_mesh.centroid, dtype="float64"))
            self.cow_mesh.apply_scale(2.0 / cow_radius)

            self.unit_cube_mesh = copy.deepcopy(models.unit_cube())

class CgalPoissonSurfaceReconTest(SurfaceReconTest):

    """
    Tests for :func:`pcsr.cgal_poisson_reconstruct`.
    """

    @_unittest.skipUnless(
        have_trimesh,
        "Could not import trimesh module for mesh/point-cloud generation."
    )
    def test_recon(self):
        """
        Unit-test for :func:`pcsr.cgal_poisson_reconstruct`.
        """
        from trimesh.sample import sample_surface_even
        from trimesh.proximity import closest_point
        from trimesh.io.export import export_mesh

        from ._pcsr import _cgal_poisson_reconstruct

        print("Generating point cloud...")
        mesh = self.unit_cube_mesh
        # export_mesh(mesh, "mesh.ply")
        points, fidx = sample_surface_even(mesh, 100000)
        normals = -mesh.face_normals[fidx]

        print("Calling CGAL Poisson surface reconstruction...")
        vertices, faces = _cgal_poisson_reconstruct(points, normals)
        print("vertices=\n%s" % vertices)
        print("faces=\n%s" % faces)
        recon_mesh = _trimesh.Trimesh(vertices=vertices, faces=faces)
        # export_mesh(recon_mesh, "mesh_recon.ply")
        c, d, fidx2 = closest_point(mesh, vertices)
        self.assertTrue(
            _np.allclose(0, d, atol=0.05)
        )

__all__ = [s for s in dir() if not s.startswith('_')]

_unittest.main(__name__)
