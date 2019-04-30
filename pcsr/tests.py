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

def export_point_cloud(
    file_name,
    point_cloud,
    exclude_fields=[],
    convert_zyx_to_xyz=False,
    logger=None
):
    """
    """
    import copy
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk

    if logger is None:
        import logging
        logger = logging.getLogger(str(__name__ + "." + "export_point_cloud"))

    logger.info("Creating vtk.vtkUnstructuredGrid...")
    vtk_ug = vtk.vtkUnstructuredGrid()
    np_points = point_cloud["coordinate"].copy()
    if convert_zyx_to_xyz:
        np_points = np_points[:, ::-1].copy()
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(np_points, deep=1))
    vtk_ug.SetPoints(points)
    # cells = vtk.vtkCellArray()

    vtk_point_data = vtk_ug.GetPointData()
    names = list(point_cloud.dtype.names)
    names.remove("coordinate")
    for name in names:
        if name not in exclude_fields:
            np_ary = point_cloud[name].copy()
            if convert_zyx_to_xyz and (len(np_ary.shape) == 2) and (np_ary.shape[1] == 3):
                # Change vector from z,y,x to x,y,z
                np_ary = np_ary[:, ::-1].copy()
            vtk_ary = numpy_to_vtk(np_ary, deep=True)
            vtk_ary.SetName(name)
            logger.debug("Adding point data %s to vtkPolyData.", str(vtk_ary.GetName()))
            vtk_point_data.AddArray(vtk_ary)
        else:
            logger.debug("Excluded field %s.", str(name))

    logger.info("Creating vtk.vtkUnstructuredGrid...done.")

    logger.info(
        "Writing vtk.vtkUnstructuredGrid to file %s, num points=%8d...",
        file_name,
        np_points.shape[0]
    )
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(file_name)
    writer.SetInputData(vtk_ug)
    writer.Write()
    logger.info(
        "Writing vtk.vtkUnstructuredGrid to file %s, num points=%8d...done.",
        file_name,
        np_points.shape[0]
    )
    del writer, points, np_points, vtk_ug

class SurfaceReconTest(_unittest.TestCase):

    """
    Base class for surface reconstruction.
    """

    def export_point_set(self, file_name, points, normals=None):
        """
        """
        dtyp = [("coordinate", (_np.float64, 3)), ]
        if normals is not None:
            dtyp.append(("normal", (_np.float64, 3)))
        point_cloud_ary = _np.zeros((points.shape[0],), dtype=dtyp)
        point_cloud_ary["coordinate"] = points
        if normals is not None:
            point_cloud_ary["normal"] = normals
        export_point_cloud(file_name, point_cloud_ary)

    def export_mesh(self, file_name, mesh):
        """
        """
        import trimesh
        if hasattr(trimesh, "io"):
            from trimesh.io.export import export_mesh
        else:
            from trimesh.exchange.export import export_mesh
        export_mesh(mesh, file_name)


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

class CgalWlopRegularizationTest(SurfaceReconTest):

    """
    Tests for :func:`pcsr.cgal_poisson_reconstruct`.
    """

    @_unittest.skipUnless(
        have_trimesh,
        "Could not import trimesh module for mesh/point-cloud generation."
    )
    def test_regularize(self):
        """
        Unit-test for :func:`pcsr.cgal_wlop_regularize`.
        """
        from trimesh.sample import sample_surface_even
        from trimesh.proximity import closest_point

        from . import cgal_wlop_regularize

        print("Generating point cloud...")
        mesh = self.unit_cube_mesh
        self.export_mesh("mesh.ply", mesh)
        g_sigma = 0.035
        points, fidx = sample_surface_even(mesh, 10000)
        for i in range(10):
            pts, fi = sample_surface_even(mesh, 10000)
            points = _np.vstack((points, pts))
            fidx = _np.hstack((fidx, fi))
            del pts, fi

        self.export_point_set("mesh_points.vtu", points)
        points += _np.random.normal(loc=0.0, scale=g_sigma, size=points.shape)
        self.export_point_set("mesh_points_noise.vtu", points)

        print("Calling CGAL WLOP regularization...")
        pts, nrms = \
            cgal_wlop_regularize(
                points,
                select_percentage=25.0,
                neighbour_radius=0.5,
                number_of_iterations=35
            )
        self.export_point_set("mesh_points_noise_regularized.vtu", pts, nrms)
        print("pts=\n%s" % pts)
        print("nrms=\n%s" % nrms)
        c, d, fidx2 = closest_point(mesh, pts)
        self.assertTrue(
            _np.allclose(0, d, atol=0.25)
        )

class CgalPoissonSurfaceReconTest(SurfaceReconTest):

    """
    Tests for :func:`pcsr.cgal_poisson_reconstruct`.
    """

    def export_mesh(self, mesh, file_name):
        """
        """
        import trimesh
        if hasattr(trimesh, "io"):
            from trimesh.io.export import export_mesh
        else:
            from trimesh.exchange.export import export_mesh
        export_mesh(mesh, file_name)

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

        from . import cgal_poisson_reconstruct

        print("Generating point cloud...")
        mesh = self.unit_cube_mesh
        # self.export_mesh("mesh.ply", mesh)
        points, fidx = sample_surface_even(mesh, 100000)
        normals = -mesh.face_normals[fidx]

        print("Calling CGAL Poisson surface reconstruction...")
        vertices, faces = cgal_poisson_reconstruct(points, normals)
        print("vertices=\n%s" % vertices)
        print("faces=\n%s" % faces)
        recon_mesh = _trimesh.Trimesh(vertices=vertices, faces=faces)
        # self.export_mesh("mesh_recon.ply", recon_mesh)
        c, d, fidx2 = closest_point(mesh, vertices)
        self.assertTrue(
            _np.allclose(0, d, atol=0.05)
        )

__all__ = [s for s in dir() if not s.startswith('_')]

_unittest.main(__name__)
