
#include <boost/version.hpp>

#if BOOST_VERSION < 106500
#include <boost/python/numeric.hpp>
#else  // #if BOOST_VERSION < 106500
#include <boost/python/numpy.hpp>
#endif

#include <pcsr/pcsr_cgal.h>
#include <boost/version.hpp>
#include <boost/python.hpp>

namespace pcsr
{

#if BOOST_VERSION < 106500

namespace bpnp = ::boost::python::numeric;
using BpNumpyArray = typename ::boost::python::numeric::array;

inline void set_module_and_type_wrapper()
{
    BpNumpyArray::set_module_and_type("numpy", "ndarray");
}

#else

namespace bpnp = ::boost::python::numpy;
using BpNumpyArray = typename boost::python::numpy::ndarray;

inline void set_module_and_type_wrapper() {}

#endif

template <typename TStlVec>
std::shared_ptr<TStlVec>
array_to_stlvec(boost::python::object threeDObj)
{
    typedef TStlVec StlVec;
    typedef std::shared_ptr<StlVec> StlVecPtr;
    typedef typename StlVec::value_type ThreeD;

    StlVecPtr stlVecPtr(new StlVec());
    for (std::size_t i = 0; i < boost::python::len(threeDObj); ++i)
    {
        ThreeD p;
        for (std::size_t j = 0; j < 3; ++j)
        {
            p[j] = boost::python::extract<typename ThreeD::value_type>(threeDObj[i][j])();
        }
        stlVecPtr->push_back(p);
    }
    return stlVecPtr;
}

template <typename TStlVec>
boost::python::object
stlvec_to_array(const TStlVec & stlVec)
{
    typedef TStlVec StlVec;
    typedef typename StlVec::value_type ThreeD;

    bpnp::dtype dtyp(bpnp::dtype::get_builtin<typename ThreeD::value_type>());
    boost::python::tuple shape(boost::python::make_tuple(stlVec.size(), 3));
    BpNumpyArray ary(bpnp::empty(shape, dtyp));
    for (std::size_t i = 0; i < stlVec.size(); i++)
    {
        for (std::size_t j = 0; j < 3; j++)
        {
            ary[boost::python::make_tuple(i, j)] = boost::python::object(stlVec[i][j]);
        }
    }
    return boost::python::object(ary);
}

PointStlVecPtr array_to_points(boost::python::object pointsObj)
{
    return array_to_stlvec<PointStlVec>(pointsObj);
}

VectorStlVecPtr array_to_vectors(boost::python::object vectorsObj)
{
    return array_to_stlvec<VectorStlVec>(vectorsObj);
}

boost::python::object
trimeshpair_to_tuple(const TriMeshPair & meshPair)
{
    PointStlVecPtr points;
    TriFaceStlVecPtr faces;
    points = meshPair.first;
    faces = meshPair.second;
    boost::python::object pointsObj(stlvec_to_array(*points));
    boost::python::object facesObj(stlvec_to_array(*faces));

    return boost::python::make_tuple(pointsObj, facesObj);
}

boost::python::object
points_normal_pair_to_tuple(const PointVectorStlVecPair & pointNormalPair)
{
    PointStlVecPtr points;
    VectorStlVecPtr normals;
    points = pointNormalPair.first;
    normals = pointNormalPair.second;
    boost::python::object pointsObj(stlvec_to_array(*points));
    boost::python::object normalsObj(stlvec_to_array(*normals));

    return boost::python::make_tuple(pointsObj, normalsObj);
}

boost::python::object cgal_poisson_reconstruct(
    boost::python::object coordinate,
    boost::python::object normal
)
{
    std::cout << "Converting numpy arrays to STL vectors..." << std::endl;
    PointStlVecPtr points(array_to_points(coordinate));
    VectorStlVecPtr normals(array_to_vectors(normal));

    boost::python::object verticesFacesTuple;
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

boost::python::object cgal_jet_smooth(
    boost::python::object coordinate,
    std::size_t num_neighbours=64,
    std::size_t jet_degree_fitting=2,
    std::size_t degree_monge=2
)
{
    std::cout << "Converting numpy array to STL vectors..." << std::endl;
    PointStlVecPtr points(array_to_points(coordinate));

    boost::python::object coordinatesNormalsTuple;
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

boost::python::object cgal_bilateral_smooth(
    boost::python::object coordinate,
    std::size_t num_neighbours=64,
    double sharpness_angle=25.0,
    std::size_t jet_degree_fitting=2
)
{
    std::cout << "Converting numpy array to STL vectors..." << std::endl;
    PointStlVecPtr points(array_to_points(coordinate));

    boost::python::object coordinatesNormalsTuple;
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

boost::python::object cgal_wlop_regularize(
    boost::python::object coordinate,
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

    boost::python::object coordinatesNormalsTuple;
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

BOOST_PYTHON_MODULE(_pcsr_cgal)
{
    pcsr::bpnp::initialize();

    boost::python::def(
        "_cgal_poisson_reconstruct",
        &pcsr::cgal_poisson_reconstruct
    );

    boost::python::def(
        "_cgal_wlop_regularize",
        &pcsr::cgal_wlop_regularize
    );

    boost::python::def(
        "_cgal_jet_smooth",
        &pcsr::cgal_jet_smooth
    );

    boost::python::def(
        "_cgal_bilateral_smooth",
        &pcsr::cgal_bilateral_smooth
    );

}

