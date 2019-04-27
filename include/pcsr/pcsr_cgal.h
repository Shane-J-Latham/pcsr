#ifndef PCSR_PCSR_CGAL
#define PCSR_PCSR_CGAL

#include <pcsr/pcsr_defs.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <vector>
#include <utility>

namespace pcsr
{

class CgalPoissonSurfaceReconstructor
{
public:
    // Types
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef Kernel::Point_3 Point3;
    typedef Kernel::Vector_3 Vector3;
    typedef std::pair<Point3, Vector3> Pwn;
    typedef std::vector<Pwn> PwnStlVec;
    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
    typedef std::shared_ptr<Polyhedron> PolyhedronPtr;

    PolyhedronPtr reconstructPolyhedron(const PwnStlVec & points) const;

    TriMeshPair reconstruct(const PointStlVec & coordinate, const VectorStlVec & normal) const;
};

}

#endif
