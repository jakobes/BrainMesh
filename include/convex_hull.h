#ifndef CONVEX_HULL__H
#define CONVEX_HULL__H

// STL
#include <vector>

// CGAL
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

// Local

//using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
//using Point_3 = Kernel::Point_3;
//using CGAL_point_vector = std::vector< Point_3 >;
//using Polyhedron = CGAL::Polyhedron_3< Kernel >;

template< typename Surface, typename Point_3 = typename Surface::Point_3>
std::shared_ptr< Surface > convex_hull(std::vector<Point_3 >& point_vector)
{
    typename Surface::Polyhedron poly;
    CGAL::convex_hull_3(point_vector.begin(), point_vector.end(), poly);
    auto result = std::make_shared< Surface >(Surface(poly));
    return result;
}


#endif