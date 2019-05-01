
#include <pcsr/pcsr_cgal.h>
#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/property_map.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>

#include <vector>
#include <fstream>
#include <iostream>
// types
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

namespace pcsr
{

CgalPointSetProcessor::Point3StlVecPtr
CgalPointSetProcessor::to_point3(const PointStlVec & coordinate)
{
    Point3StlVecPtr point3(new Point3StlVec());
    for (auto cit = coordinate.begin(); cit != coordinate.end(); ++cit)
    {
        const Point3 p((*cit)[0], (*cit)[1], (*cit)[2]);
        point3->push_back(p);
    }
    return point3;
}

CgalPointSetProcessor::PwnStlVecPtr
CgalPointSetProcessor::to_pwn(const PointStlVec & coordinate, const VectorStlVec & normal)
{
    PwnStlVecPtr pwn(new PwnStlVec());
    auto cit = coordinate.begin();
    auto nit = normal.begin();
    for (; cit != coordinate.end() && nit != normal.end(); ++cit, ++nit)
    {
        const Point3 p((*cit)[0], (*cit)[1], (*cit)[2]);
        const Vector3 n((*nit)[0], (*nit)[1], (*nit)[2]);
        pwn->push_back(Pwn(p, n));
    }
    return pwn;
}

CgalPointSetProcessor::PwnStlVecPtr
CgalPointSetProcessor::to_pwn(const Point3StlVec & coordinate)
{
    PwnStlVecPtr pwn(new PwnStlVec());
    const Vector3 n(0.0, 0.0, 0.0);
    for (auto cit = coordinate.begin(); cit != coordinate.end(); ++cit)
    {
        const Point3 p((*cit)[0], (*cit)[1], (*cit)[2]);
        pwn->push_back(Pwn(p, n));
    }
    return pwn;
}

PointVectorStlVecPair
CgalPointSetProcessor::from_pwn(const PwnStlVec & pwn)
{
    PointVectorStlVecPair pvp(PointStlVecPtr(new PointStlVec()), VectorStlVecPtr(new VectorStlVec()));
    for (auto pit = pwn.begin(); pit != pwn.end(); ++pit)
    {
        Point p;
        Vector n;
        for (std::size_t j = 0; j < 3; ++j)
        {
            p[j] = pit->first[j];
            n[j] = pit->second[j];
        }
        pvp.first->push_back(p);
        pvp.second->push_back(n);
    }
    return pvp;
}

void
CgalPointSetProcessor::jet_estimate_normals(
    PwnStlVec & pwn,
    std::size_t k,
    std::size_t degree_fitting
)
{
    CGAL::jet_estimate_normals<Concurrency_tag>(
        pwn,
        k,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()).
        normal_map(CGAL::Second_of_pair_property_map<Pwn>()).
        degree_fitting(degree_fitting)
    );
}

CgalPointSetProcessor::PwnStlVecPtr
CgalPointSetProcessor::jet_estimate_normals(
    const Point3StlVec & coordinate,
    std::size_t k,
    std::size_t degree_fitting
)
{
    PwnStlVecPtr pwn(to_pwn(coordinate));
    jet_estimate_normals(*pwn);
    return pwn;
}

CgalPointSetProcessor::PwnStlVecPtr
CgalPointSetProcessor::pca_estimate_normals(
    const Point3StlVec & coordinate,
    std::size_t k
)
{
    PwnStlVecPtr pwn(to_pwn(coordinate));
    CGAL::pca_estimate_normals<Concurrency_tag>(
        *pwn,
        k,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()).
        normal_map(CGAL::Second_of_pair_property_map<Pwn>())
    );
    return pwn;
}

CgalPointSetProcessor::PwnStlVec::iterator
CgalPointSetProcessor::mst_orient_normals(
    PwnStlVec & pwn,
    std::size_t k
)
{
    return
        CGAL::mst_orient_normals(
            pwn,
            k,
            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()).
            normal_map(CGAL::Second_of_pair_property_map<Pwn>())
        );
}

/* ========================================================================= */

CgalJetPointSetSmoother::CgalJetPointSetSmoother(
    std::size_t nn,
    std::size_t jdf,
    std::size_t dm
): Inherited(),
   num_neighbours(nn),
   jet_degree_fitting(jdf),
   degree_monge(dm)
{
}

void
CgalJetPointSetSmoother::smooth(Point3StlVec & coordinate) const
{
    CGAL::jet_smooth_point_set<Concurrency_tag>(
        coordinate,
        this->num_neighbours,
        CGAL::parameters::degree_fitting(this->jet_degree_fitting).
        degree_monge(this->degree_monge)
    );
}

PointVectorStlVecPair
CgalJetPointSetSmoother::smooth(const PointStlVec & coordinate) const
{
    Point3StlVecPtr points(this->to_point3(coordinate));
    this->smooth(*points);
    const std::size_t num_neighbs = this->num_neighbours;
    PwnStlVecPtr
        pwn(
            this->jet_estimate_normals(
                *points,
                num_neighbs,
                this->jet_degree_fitting
            )
        );
    points.reset();
    this->mst_orient_normals(*pwn, num_neighbs);

    PointVectorStlVecPair pointsNormalsPair(this->from_pwn(*pwn));
    pwn.reset();
    return pointsNormalsPair;
}

/* ========================================================================= */

CgalBilateralPointSetSmoother::CgalBilateralPointSetSmoother(
    std::size_t nn,
    double sa,
    std::size_t jdf
): Inherited(),
   num_neighbours(nn),
   sharpness_angle(sa),
   jet_degree_fitting(jdf)
{
}

void
CgalBilateralPointSetSmoother::smooth(PwnStlVec & pwn) const
{
    CGAL::bilateral_smooth_point_set<Concurrency_tag>(
        pwn,
        this->num_neighbours,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()).
         normal_map(CGAL::Second_of_pair_property_map<Pwn>()).
         sharpness_angle(this->sharpness_angle)
    );
}

PointVectorStlVecPair
CgalBilateralPointSetSmoother::smooth(const PointStlVec & coordinate) const
{
    const std::size_t num_neighbs = this->num_neighbours;
    Point3StlVecPtr points(this->to_point3(coordinate));
    std::cout << "Smoothing..." << std::endl;
    PwnStlVecPtr
        pwn(
            this->jet_estimate_normals(
                *points,
                num_neighbs,
                this->jet_degree_fitting
            )
        );
    this->mst_orient_normals(*pwn, num_neighbs);
    points.reset();
    this->smooth(*pwn);
    std::cout << "Estimating normals..." << std::endl;
    this->jet_estimate_normals(
        *pwn,
        num_neighbs,
        this->jet_degree_fitting
    );
    std::cout << "Orienting normals..." << std::endl;
    this->mst_orient_normals(*pwn, num_neighbs);

    PointVectorStlVecPair pointsNormalsPair(this->from_pwn(*pwn));
    pwn.reset();
    return pointsNormalsPair;
}

/* ========================================================================= */

CgalWlopPointSetRegularizer::CgalWlopPointSetRegularizer(
    const double sp,
    const double nr,
    std::size_t noi,
    bool rus,
    std::size_t nn,
    std::size_t jdf
): Inherited(),
   select_percentage(sp),
   neighbour_radius(nr),
   number_of_iterations(noi),
   require_uniform_sampling(rus),
   num_neighbours(nn),
   jet_degree_fitting(jdf)
{
}

CgalWlopPointSetRegularizer::Point3StlVecPtr
CgalWlopPointSetRegularizer::regularize(Point3StlVec & coordinate) const
{
    Point3StlVecPtr output(new Point3StlVec());
    CGAL::wlop_simplify_and_regularize_point_set<Concurrency_tag>(
        coordinate,
        std::back_inserter(*output),
        CGAL::parameters::select_percentage(this->select_percentage).
        neighbor_radius(this->neighbour_radius).
        number_of_iterations(this->number_of_iterations).
        require_uniform_sampling(this->require_uniform_sampling)
    );
    return output;
}

PointVectorStlVecPair
CgalWlopPointSetRegularizer::regularize(const PointStlVec & coordinate) const
{
    Point3StlVecPtr points(this->to_point3(coordinate));
    points = this->regularize(*points);
    const std::size_t num_neighbs = this->num_neighbours;
    PwnStlVecPtr
        pwn(
            this->jet_estimate_normals(
                *points,
                num_neighbs,
                this->jet_degree_fitting
            )
        );
    points.reset();
    this->mst_orient_normals(*pwn, num_neighbs);

    PointVectorStlVecPair pointsNormalsPair(this->from_pwn(*pwn));
    pwn.reset();
    return pointsNormalsPair;
}

/* ========================================================================= */

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
    PwnStlVecPtr pwnPtr(this->to_pwn(coordinate, normal));
    PolyhedronPtr polyPtr;
    polyPtr = this->reconstructPolyhedron(*pwnPtr);
    pwnPtr.reset();

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
