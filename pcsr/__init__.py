"""
Point Cloud Surface Reconstruction (PCSR).
"""
from . import _version

__version__ = _version.get_versions()['version']


def cgal_jet_smooth(
    points,
    num_neighbours=18,
    jet_degree_fitting=2,
    degree_monge=2
):
    """
    Wrapper for CGAL::jet_smooth_point_set followed by normal
    estimation and orientation using CGAL::jet_estimate_normals
    and CGAL::mst_orient_normals, respectively.

    :type points: :obj:`numpy.ndarray`
    :param points: :obj:`numpy.ndarray` A :samp:`(N, 3)` shaped array of 3D coordinates.
    :type num_neighbours: :obj:`int`
    :param num_neighbours: Number of nearest neighbours used in the jet fitting and
       the jet normal estimation and the normal orientation estimation.
    :type jet_degree_fitting: :obj:`int`
    :param jet_degree_fitting: Order of polynomial surface fit used in
       the CGAL::jet_smooth_point_set smoothing and CGAL::jet_estimate_normals normal estimation.
    :type degree_monge: :obj:`int`
    :param degree_monge: Order of polynomial surface fit used in
    :rtype: :obj:`tuple`
    :return: A :samp:`(points, normals)` pair, where :samp:`points` is
       a :samp:`(S, 3)` shaped array of point coordinates
       and :samp:`normals` is a :samp:`(S, 3)` shaped array of oriented normals.
    """
    from . import _pcsr_cgal
    return \
        _pcsr_cgal._cgal_jet_smooth(
            points,
            int(num_neighbours),
            int(jet_degree_fitting),
            int(degree_monge)
        )


def cgal_bilateral_smooth(
    points,
    num_neighbours=18,
    sharpness_angle=25.0,
    jet_degree_fitting=2
):
    """
    Wrapper for CGAL::bilateral_smooth_point_set followed by normal
    estimation and orientation using CGAL::jet_estimate_normals
    and CGAL::mst_orient_normals, respectively.

    :type points: :obj:`numpy.ndarray`
    :param points: :obj:`numpy.ndarray` A :samp:`(N, 3)` shaped array of 3D coordinates.
    :type num_neighbours: :obj:`int`
    :param num_neighbours: Number of nearest neighbours used in the bilateral smoothing and
       the jet normal estimation and the normal orientation estimation.
    :type sharpness_angle: :obj:`float`
    :param sharpness_angle: Angle degrees in range :samp:`[0, 90]`
        controlling the sharpness of the result, larger angle gives smoother result.
    :type jet_degree_fitting: :obj:`int`
    :param jet_degree_fitting: Order of polynomial surface fit used in
       the CGAL::jet_estimate_normals normal estimation.
    :rtype: :obj:`tuple`
    :return: A :samp:`(points, normals)` pair, where :samp:`points` is
       a :samp:`(S, 3)` shaped array of point coordinates
       and :samp:`normals` is a :samp:`(S, 3)` shaped array of oriented normals.
    """
    from . import _pcsr_cgal
    return \
        _pcsr_cgal._cgal_bilateral_smooth(
            points,
            int(num_neighbours),
            float(sharpness_angle),
            int(jet_degree_fitting)
        )


def cgal_wlop_regularize(
    points,
    select_percentage=10.0,
    neighbour_radius=-1.0,
    number_of_iterations=100,
    require_uniform_sampling=False,
    num_neighbours=18,
    jet_degree_fitting=2
):
    """
    Wrapper for CGAL::wlop_simplify_and_regularize_point_set followed by normal
    estimation and orientation using CGAL::jet_estimate_normals
    and CGAL::mst_orient_normals, respectively.

    :type points: :obj:`numpy.ndarray`
    :param points: :obj:`numpy.ndarray` A :samp:`(N, 3)` shaped array of 3D coordinates.
    :type select_percentage: :obj:`float`
    :param select_percentage: The percentage of regularized points to return.
    :type neighbour_radius: :obj:`float`
    :param neighbour_radius: Spherical neighbourhood radius for the WLOP algorithm, auto-calculated
       if less than or equal to zero.
    :type require_uniform_sampling: :obj:`bool`
    :param require_uniform_sampling: If :obj:`True` attempts to perform uniform
       point regularization.
    :type num_neighbours: :obj:`int`
    :param num_neighbours: Number of nearest neighbours used in the normal estimation
       and orientation.
    :type jet_degree_fitting: :obj:`int`
    :param jet_degree_fitting: Order of polynomial surface fit used in
       the CGAL::jet_estimate_normals normal estimation.
    :rtype: :obj:`tuple`
    :return: A :samp:`(points, normals)` pair, where :samp:`points` is
       a :samp:`(S, 3)` shaped array of point coordinates
       and :samp:`normals` is a :samp:`(S, 3)` shaped array of oriented normals.
    """
    from . import _pcsr_cgal
    return \
        _pcsr_cgal._cgal_wlop_regularize(
            points,
            float(select_percentage),
            float(neighbour_radius),
            int(number_of_iterations),
            bool(require_uniform_sampling),
            int(num_neighbours),
            int(jet_degree_fitting)
        )


def cgal_poisson_reconstruct(
    points,
    normals
):
    """
    Wrapper for CGAL::poisson_surface_reconstruction_delaunay point set surface reconstruction.

    :type points: :obj:`numpy.ndarray`
    :param points: :obj:`numpy.ndarray` A :samp:`(N, 3)` shaped array of 3D coordinates.
    :type normals: :obj:`numpy.ndarray`
    :param normals: :obj:`numpy.ndarray` A :samp:`(N, 3)` shaped array of oriented 3D normal
       vectors.
    :rtype: :obj:`tuple`
    :return: A :samp:`(vertices, faces)` pair, where :samp:`vertices` is
       a :samp:`(S, 3)` shaped array of mesh vertex coordinates
       and :samp:`faces` is a :samp:`(S, 3)` shaped array of :samp:`vertices` indices,
       such that, :samp:`faces[f]` defines the triangular face with
       vertices :samp:`vertices[faces[f][0]]`, :samp:`vertices[faces[f][1]]`
       and :samp:`vertices[faces[f][2]]`.
    """
    from . import _pcsr_cgal
    return \
        _pcsr_cgal._cgal_poisson_reconstruct(
            points,
            normals
        )
