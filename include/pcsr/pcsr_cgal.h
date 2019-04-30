#ifndef PCSR_PCSR_CGAL
#define PCSR_PCSR_CGAL

#include <pcsr/pcsr_defs.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <vector>
#include <utility>

namespace pcsr
{

class CgalPointSetProcessor
{
public:
    // Types
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef Kernel::Point_3 Point3;
    typedef std::vector<Point3> Point3StlVec;
    typedef std::shared_ptr<Point3StlVec> Point3StlVecPtr;
    typedef Kernel::Vector_3 Vector3;
    typedef std::pair<Point3, Vector3> Pwn;
    typedef std::vector<Pwn> PwnStlVec;
    typedef std::shared_ptr<PwnStlVec> PwnStlVecPtr;
    // Concurrency
#ifdef CGAL_LINKED_WITH_TBB
    typedef CGAL::Parallel_tag Concurrency_tag;
#else
    typedef CGAL::Sequential_tag Concurrency_tag;
#endif

    static Point3StlVecPtr to_point3(const PointStlVec & pts);

    static PwnStlVecPtr to_pwn(const PointStlVec & pts, const VectorStlVec & nrms);

    static PwnStlVecPtr to_pwn(const Point3StlVec & pts);

    static PointVectorStlVecPair from_pwn(const PwnStlVec & pwn);

    static
    PwnStlVecPtr
    jet_estimate_normals(
        const Point3StlVec & coordinate,
        std::size_t k = 18,
        std::size_t degree_fitting = 2
    );

    static
    PwnStlVec::iterator
    mst_orient_normals(
        PwnStlVec & pwn,
        std::size_t k = 18
    );
};

class CgalWlopPointSetRegularizer: public CgalPointSetProcessor
{
public:
    // Types
    typedef CgalPointSetProcessor Inherited;
    using Inherited::Kernel;
    using Inherited::Point3;
    using Inherited::Vector3;
    using Inherited::Pwn;
    using Inherited::PwnStlVec;
    using Inherited::PwnStlVecPtr;

    CgalWlopPointSetRegularizer(
        const double select_percentage=10.0,
        const double neighbour_radius=-1.0,
        std::size_t number_of_iterations=32,
        bool require_uniform_sampling=false,
        std::size_t num_neighbours=18,
        std::size_t jet_degree_fitting=2
    );

    Point3StlVecPtr regularize(Point3StlVec & coordinate) const;

    PointVectorStlVecPair regularize(const PointStlVec & coordinate) const;

    double select_percentage;
    double neighbour_radius;
    std::size_t number_of_iterations;
    bool require_uniform_sampling;
    std::size_t num_neighbours;
    std::size_t jet_degree_fitting;

};

class CgalPoissonSurfaceReconstructor: public CgalPointSetProcessor
{
public:
    // Types
    typedef CgalPointSetProcessor Inherited;
    using Inherited::Kernel;
    using Inherited::Point3;
    using Inherited::Vector3;
    using Inherited::Pwn;
    using Inherited::PwnStlVec;
    using Inherited::PwnStlVecPtr;
    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
    typedef std::shared_ptr<Polyhedron> PolyhedronPtr;

    PolyhedronPtr reconstructPolyhedron(const PwnStlVec & points) const;

    TriMeshPair reconstruct(const PointStlVec & coordinate, const VectorStlVec & normal) const;
};

}

#endif
