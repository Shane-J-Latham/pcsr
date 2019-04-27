
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

}

BOOST_PYTHON_MODULE(_pcsr)
{
    pcsr::bpnp::initialize();

    boost::python::def(
        "_cgal_poisson_reconstruct",
        &pcsr::cgal_poisson_reconstruct
    );
}

