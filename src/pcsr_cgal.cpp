
#include <pcsr/pcsr_cgal.h>
#include <CGAL/poisson_surface_reconstruction.h>

namespace pcsr
{

CgalPoissonSurfaceReconstructor::PolyhedronPtr
CgalPoissonSurfaceReconstructor::reconstructPolyhedron(const PwnStlVec & points) const
{
    double average_spacing =
        CGAL::compute_average_spacing<CGAL::Sequential_tag>(
            points,
            6,
            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
        );

    PolyhedronPtr meshPtr(new Polyhedron());
    double  sm_angle = 10.0;
    double  sm_radius = 20.0;
    double  sm_distance = 0.25;
    if (
        CGAL::poisson_surface_reconstruction_delaunay(
            points.begin(), points.end(),
            CGAL::First_of_pair_property_map<Pwn>(),
            CGAL::Second_of_pair_property_map<Pwn>(),
            *meshPtr,
            0.5 * average_spacing,
            sm_angle,
            sm_radius,
            sm_distance,
            CGAL::Non_manifold_tag()
        )
    )
    {
    }
    else
    {

    }
    return meshPtr;
}

TriMeshPair
CgalPoissonSurfaceReconstructor::reconstruct(const PointStlVec & coordinate, const VectorStlVec & normal) const
{
    PwnStlVec points;
    auto cit = coordinate.begin();
    auto nit = normal.begin();
    for (; cit != coordinate.end() && nit != normal.end(); ++cit, ++nit)
    {
        const Point3 p((*cit)[0], (*cit)[1], (*cit)[2]);
        const Vector3 n((*nit)[0], (*nit)[1], (*nit)[2]);
        points.push_back(Pwn(p, n));
    }
    PolyhedronPtr polyPtr;
    polyPtr = this->reconstructPolyhedron(points);

    PointStlVecPtr pointsPtr(new PointStlVec());
    TriFaceStlVecPtr triFacesPtr(new TriFaceStlVec());

    for (auto vit = polyPtr->vertices_begin(); vit != polyPtr->vertices_end(); ++vit)
    {
        Point p;
        p[0] = (vit->point())[0];
        p[1] = (vit->point())[1];
        p[2] = (vit->point())[2];
        pointsPtr->push_back(p);
    }
    for (auto fit = polyPtr->facets_begin(); fit != polyPtr->facets_end(); ++fit)
    {
        auto idit = fit->facet_begin();
        CGAL_assertion( CGAL::circulator_size(idit) == 3);
        TriFace face;
        for (std::size_t v = 0; v < 3; ++v, ++idit)
        {
            face[v] = std::distance(polyPtr->vertices_begin(), idit->vertex());
        }
        triFacesPtr->push_back(face);
    }
    TriMeshPair triMeshPair(pointsPtr, triFacesPtr);
    return triMeshPair;
}

}
