
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <pcsr/pcsr_cgal.h>

namespace pcsr
{

namespace py = pybind11;

template <typename TStlVec>
std::shared_ptr<TStlVec>
array_to_stlvec(py::object threeDObj)
{
    typedef TStlVec StlVec;
    typedef std::shared_ptr<StlVec> StlVecPtr;
    typedef typename StlVec::value_type ThreeD;

    StlVecPtr stlVecPtr(new StlVec());
    py::object j_objs[] = {py::cast(long(0)), py::cast(long(1)), py::cast(long(2))};
    for (long i = 0; i < py::len(threeDObj); ++i)
    {
        ThreeD p;
        py::object i_obj(py::cast(i));
        for (long j = 0; j < 3; ++j)
        {
            p[j] = py::cast<typename ThreeD::value_type>(threeDObj[i_obj][j_objs[j]]);
        }
        stlVecPtr->push_back(p);
    }
    return stlVecPtr;
}

template <typename TStlVec>
py::object
stlvec_to_array(const TStlVec & stlVec)
{
    typedef TStlVec StlVec;
    typedef typename StlVec::value_type ThreeD;

    py::dtype dtyp(py::dtype::of<typename ThreeD::value_type>());

    size_t shape[2]{stlVec.size(), 3};
    // py::array ary(dtyp, shape);
    auto ary = py::array_t<typename ThreeD::value_type>(shape);

    for (std::size_t i = 0; i < stlVec.size(); i++)
    {
        for (std::size_t j = 0; j < 3; j++)
        {
            ary.mutable_at(i, j) = stlVec[i][j];
        }
    }
    return py::object(ary);
}

PointStlVecPtr array_to_points(py::object pointsObj)
{
    return array_to_stlvec<PointStlVec>(pointsObj);
}

VectorStlVecPtr array_to_vectors(py::object vectorsObj)
{
    return array_to_stlvec<VectorStlVec>(vectorsObj);
}

py::object
trimeshpair_to_tuple(const TriMeshPair & meshPair)
{
    PointStlVecPtr points;
    TriFaceStlVecPtr faces;
    points = meshPair.first;
    faces = meshPair.second;
    py::object pointsObj(stlvec_to_array(*points));
    py::object facesObj(stlvec_to_array(*faces));

    return py::make_tuple(pointsObj, facesObj);
}

py::object
points_normal_pair_to_tuple(const PointVectorStlVecPair & pointNormalPair)
{
    PointStlVecPtr points;
    VectorStlVecPtr normals;
    points = pointNormalPair.first;
    normals = pointNormalPair.second;
    py::object pointsObj(stlvec_to_array(*points));
    py::object normalsObj(stlvec_to_array(*normals));

    return py::make_tuple(pointsObj, normalsObj);
}

py::object cgal_poisson_reconstruct(
    py::object coordinate,
    py::object normal
)
{
    std::cout << "Converting numpy arrays to STL vectors..." << std::endl;
    PointStlVecPtr points(array_to_points(coordinate));
    VectorStlVecPtr normals(array_to_vectors(normal));

    py::object verticesFacesTuple;
    {
        std::cout << "Reconstructing surface..." << std::endl;
        CgalPoissonSurfaceReconstructor surfaceReconstructor;
        TriMeshPair meshPair(surfaceReconstructor.reconstruct(*points, *normals));
        points.reset();
        normals.reset();
        std::cout << "Converting mesh STL vectors to numpy array tuple..." << std::endl;
        verticesFacesTuple = trimeshpair_to_tuple(meshPair);
    }

    return verticesFacesTuple;
}

py::object cgal_jet_smooth(
    py::object coordinate,
    std::size_t num_neighbours=64,
    std::size_t jet_degree_fitting=2,
    std::size_t degree_monge=2
)
{
    std::cout << "Converting numpy array to STL vectors..." << std::endl;
    PointStlVecPtr points(array_to_points(coordinate));

    py::object coordinatesNormalsTuple;
    {
        std::cout << "Jet smoothing point set..." << std::endl;
        CgalJetPointSetSmoother smoother;
        smoother.num_neighbours = num_neighbours;
        smoother.jet_degree_fitting = jet_degree_fitting;
        smoother.degree_monge = degree_monge;
        PointVectorStlVecPair pointsNormalsPair(smoother.smooth(*points));
        points.reset();
        std::cout << "Converting mesh STL vectors to numpy array tuple..." << std::endl;
        coordinatesNormalsTuple = points_normal_pair_to_tuple(pointsNormalsPair);
    }

    return coordinatesNormalsTuple;
}

py::object cgal_bilateral_smooth(
    py::object coordinate,
    std::size_t num_neighbours=64,
    double sharpness_angle=25.0,
    std::size_t jet_degree_fitting=2
)
{
    std::cout << "Converting numpy array to STL vectors..." << std::endl;
    PointStlVecPtr points(array_to_points(coordinate));

    py::object coordinatesNormalsTuple;
    {
        std::cout << "Bilateral smoothing point set..." << std::endl;
        CgalBilateralPointSetSmoother smoother;
        smoother.num_neighbours = num_neighbours;
        smoother.jet_degree_fitting = jet_degree_fitting;
        smoother.sharpness_angle = sharpness_angle;
        PointVectorStlVecPair pointsNormalsPair(smoother.smooth(*points));
        points.reset();
        std::cout << "Converting mesh STL vectors to numpy array tuple..." << std::endl;
        coordinatesNormalsTuple = points_normal_pair_to_tuple(pointsNormalsPair);
    }

    return coordinatesNormalsTuple;
}

py::object cgal_wlop_regularize(
    py::object coordinate,
    const double select_percentage=10.0,
    const double neighbour_radius=-1.0,
    std::size_t number_of_iterations=32,
    bool require_uniform_sampling=false,
    std::size_t num_neighbours=18,
    std::size_t jet_degree_fitting=2
)
{
    std::cout << "Converting numpy array to STL vectors..." << std::endl;
    PointStlVecPtr points(array_to_points(coordinate));

    py::object coordinatesNormalsTuple;
    {
        std::cout << "Regularizing point set..." << std::endl;
        CgalWlopPointSetRegularizer regularizer;
        regularizer.select_percentage = select_percentage;
        regularizer.neighbour_radius = neighbour_radius;
        regularizer.number_of_iterations = number_of_iterations;
        regularizer.require_uniform_sampling = require_uniform_sampling;
        regularizer.num_neighbours = num_neighbours;
        regularizer.jet_degree_fitting = jet_degree_fitting;
        PointVectorStlVecPair pointsNormalsPair(regularizer.regularize(*points));
        points.reset();
        std::cout << "Converting mesh STL vectors to numpy array tuple..." << std::endl;
        coordinatesNormalsTuple = points_normal_pair_to_tuple(pointsNormalsPair);
    }

    return coordinatesNormalsTuple;
}

}

PYBIND11_MODULE(_pcsr_cgal, m)
{
    m.def(
        "_cgal_poisson_reconstruct",
        &pcsr::cgal_poisson_reconstruct
    );

    m.def(
        "_cgal_wlop_regularize",
        &pcsr::cgal_wlop_regularize
    );

    m.def(
        "_cgal_jet_smooth",
        &pcsr::cgal_jet_smooth
    );

    m.def(
        "_cgal_bilateral_smooth",
        &pcsr::cgal_bilateral_smooth
    );
}
