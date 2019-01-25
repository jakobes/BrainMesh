#ifndef CGALSurface_H
#define CGALSurface_H

#define BOOST_PARAMETER_MAX_ARITY 12

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "utils.h"

#include <iostream>
#include <vector>

#include <CGAL/Surface_mesh.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

// Isotropic remeshing
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>

// Edge collapse -- Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

// copy_face_graph
#include <CGAL/boost/graph/copy_face_graph.h>

// points inside/outside
#include <CGAL/Side_of_triangle_mesh.h>

// Corefine and compute difference
#include <CGAL/Polygon_mesh_processing/corefinement.h>

// Poisson reconstruction
#include <CGAL/poisson_surface_reconstruction.h>        // Why do I need this one here


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

typedef CGAL::Side_of_triangle_mesh<Mesh, Kernel> inside; //  TODO: Cnange name, unintuitive


class CGALSurface {

    private:
        Mesh mesh;

    public:
        CGALSurface(); // empty constructor

        CGALSurface(const std::string f);

        CGALSurface(Polyhedron &);

        /* CGALSurface( Implicit_function implicit_function, */
        /*              double bounding_sphere_radius, */
        /*              double angular_bound, */
        /*              double radius_bound, */
        /*              double distance_bound); */

        /* template<typename Implicit_function> */
        /* CGALSurface(Implicit_function implicit_function, */
        /*          double bounding_sphere_radius, */
        /*          double angular_bound=30., */
        /*          double radius_bound=0.1 , */
        /*          double distance_bound=0.1   ); */


        template<typename Polyhedron_3>
        void get_polyhedron(Polyhedron_3 &polyhedron_3);

        /* void operator^=(CGALSurface &other); */
        /* void operator+=(CGALSurface &other); */
        /* void operator-=(CGALSurface &other); */
        void surface_intersection(CGALSurface &other);
        void surface_union(CGALSurface &other);
        void surface_difference(CGALSurface &other);


        Mesh& get_mesh();

        int fill_holes();

        bool triangulate_faces();

        void stitch_borders();

    /*     void insert_surface(CGALSurface& surface); */    // TODO

        void isotropic_remeshing(const double, const unsigned int, const bool);

        void adjust_boundary(const double);

        void smooth_laplacian(const double);

        void smooth_taubin(const size_t);

        template<typename InputIterator>
        void adjusting_boundary_region(InputIterator, InputIterator, const double);

        template<typename InputIterator>
        void smooth_laplacian_region(InputIterator, InputIterator, const int);

        std::vector<vertex_descriptor> points_inside(CGALSurface &);

        std::vector<vertex_descriptor> points_outside(CGALSurface& other);

        bool self_intersections();

        int num_self_intersections();

        void save(const std::string);

        int collapse_edges(const double stop_ratio);

        void preprocess(const double, const int);

        void fair(std::vector<vertex_descriptor> &);

        int num_faces() const;

        int num_edges() const;

        int num_vertices() const;

        void split_edges(const double target_edge_length);

        // TODO: What does this do?
        void clear() { mesh.clear(); }

        // TODO: Add this for tomorrow
        /* void fix_close_junctures(double c); */

        /* void make_cylinder( double x0, double y0, double  z0,  double x1, double y1, double z1, double radius,  int number_of_segments=360) ; */

        /* void make_cone( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0 , double r1,  int number_of_segments=360) ; */
        // TODO
        /* void make_cube( double x0, double y0, double  z0,  double x1, double y1, double z1); */
        /* void make_sphere( double x0, double y0, double  z0, double r0); */
        /* void insert_mesh(CGALSurface& surf){mesh+=surf.get_mesh();} */
        /* void insert_points(std::vector<Point_3>& points) ; */

    /*     template<typename Implicit_function> */
    /*     CGALSurface(Implicit_function implicit_function, */
    /*          double bounding_sphere_radius, */
    /*          double angular_bound=30., */
    /*          double radius_bound=0.1 , */
    /*          double distance_bound=0.1   ); */

        ~CGALSurface(){}
};

CGALSurface::CGALSurface(const std::string f) {
    utils::read_off(mesh, f);
}


CGALSurface::CGALSurface(Polyhedron &polyhedron) {
    CGAL::copy_face_graph(polyhedron, mesh);
}


template<typename Polyhedron_3>
void CGALSurface::get_polyhedron(Polyhedron_3 &polyhedron_3) {
    CGAL::copy_face_graph(mesh, polyhedron_3);
}


Mesh& CGALSurface::get_mesh() {
    return mesh;
}

int CGALSurface::fill_holes() {
    unsigned int nb_holes = 0;
    for (auto h: halfedges(mesh)) {
        if(is_border(h, mesh)) {
            std::vector<face_descriptor> patch_facets;
            std::vector<vertex_descriptor> patch_vertices;
            bool success = CGAL::cpp11::get<0>(
            CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(mesh, h, std::back_inserter(patch_facets),
                std::back_inserter(patch_vertices),
                CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).geom_traits(Kernel())) );
            std::cout << "* Number of facets in constructed patch: " << patch_facets.size() << std::endl;
            std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
            std::cout << "  Is fairing successful: " << success << std::endl;
            nb_holes++;
        }
    }
    std::cout << std::endl;
    std::cout << nb_holes << " holes have been filled" << std::endl;
    return nb_holes;
}


bool CGALSurface::triangulate_faces() {
    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    for (auto fit: faces(mesh)) {
        if (next(next(halfedge(fit,  mesh), mesh), mesh) != prev(halfedge(fit, mesh), mesh)) {
            std::cerr << "Error: non-triangular face left in mesh." << std::endl;
        }
    }
    return true;
}


void CGALSurface::stitch_borders() {
    CGAL::Polygon_mesh_processing::stitch_borders(mesh);
}


void CGALSurface::isotropic_remeshing(
    const double target_edge_length,
    const unsigned int nb_iter,
    const bool protect_border) {
    CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length, mesh);
    CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mesh),
            target_edge_length, mesh,
            CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
            .protect_constraints(protect_border));
}


void CGALSurface::adjust_boundary(const double c) {
    Mesh::Vertex_range::iterator vb = mesh.vertices().begin(), ve = mesh.vertices().end();
    CGALSurface::adjusting_boundary_region(vb, ve, c);
}


void CGALSurface::smooth_laplacian(const double c) {
    Mesh::Vertex_range::iterator vb = mesh.vertices().begin(), ve = mesh.vertices().end();
    CGALSurface::smooth_laplacian_region(vb, ve, c);
}


void CGALSurface::smooth_taubin(const size_t nb_iter) {
    for (size_t i = 0; i < nb_iter; ++i) {
        this->smooth_laplacian(0.8);
        this->smooth_laplacian(-0.805);
    }
}


template<typename InputIterator>
void CGALSurface::adjusting_boundary_region(
        InputIterator begin,
        InputIterator end,
        const double c) {
    std::vector<std::pair<vertex_descriptor, Point> > smoothed; //rename
    for ( ; begin != end; ++begin) {
        Vector delta = CGAL::Polygon_mesh_processing::compute_vertex_normal(*begin,mesh);
        Point p = mesh.point(*begin) + c*delta;
        smoothed.push_back(std::make_pair(*begin, p));
    }
    for (auto &s: smoothed) {
        mesh.point(s.first) = s.second;
    }
}


template<typename InputIterator>
void CGALSurface::smooth_laplacian_region(InputIterator begin, InputIterator end, const int c) {
    std::vector<std::pair<vertex_descriptor, Point> > smoothed;
    for ( ; begin != end; ++begin) {
        Point current = mesh.point(*begin);
        Vector delta = CGAL::NULL_VECTOR;
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(*begin), mesh), done(vbegin);
        do {
            delta += Vector(mesh.point(*vbegin) - current);
            *vbegin++;
        }
        while(vbegin != done);
        Point p = current + c*delta/mesh.degree(*begin);
        smoothed.push_back(std::make_pair(*begin, p));
    }
    for (auto &s: smoothed) {
        mesh.point(s.first) = s.second;
    }
}

bool CGALSurface::self_intersections() {
    return CGAL::Polygon_mesh_processing::does_self_intersect(mesh,
                CGAL::Polygon_mesh_processing::parameters::vertex_point_map(
                    get(CGAL::vertex_point, mesh)));
}


int CGALSurface::num_self_intersections() {
    std::vector< std::pair<face_descriptor, face_descriptor> > intersected_tris;
    CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(intersected_tris));
    return intersected_tris.size();  // Could actually return the triangles themselves
}


void CGALSurface::save(const std::string outpath) {
    utils::save_off(mesh, outpath);
};


int CGALSurface::collapse_edges(const double stop_ratio) {
    namespace SMS = CGAL::Surface_mesh_simplification;
    SMS::Count_ratio_stop_predicate<Mesh> stop(stop_ratio);

    const int r = SMS::edge_collapse(
        mesh,
        stop,
        CGAL::parameters::get_cost(SMS::Edge_length_cost <Mesh>())
            .get_placement(SMS::Midpoint_placement<Mesh>()));
            /* .visitor(vis)); */
    return r;
}

void CGALSurface::preprocess(const double target_edge_length, const int nb_iter) {
    CGALSurface::triangulate_faces();
    CGALSurface::isotropic_remeshing(target_edge_length, nb_iter, false);
};


std::vector<vertex_descriptor> CGALSurface::points_inside(CGALSurface &other) {
    std::vector<vertex_descriptor> result;
    inside inside_poly2(other.get_mesh());
    for (const auto &v_it: mesh.vertices()) {
        const CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
        if (res == CGAL::ON_BOUNDED_SIDE or res == CGAL::ON_BOUNDARY) {
            result.push_back(v_it);
        }
    }
    return result;
}


std::vector<vertex_descriptor> CGALSurface::points_outside(CGALSurface& other) {
    std::vector<vertex_descriptor> result;
    inside inside_poly2(other.get_mesh());

    for (const auto &v_it: mesh.vertices()) {
        const CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
        if (res == CGAL::ON_UNBOUNDED_SIDE) {
              result.push_back(v_it);
        }
    }
    return result;
}


void CGALSurface::fair(std::vector<vertex_descriptor> &vector) {
    CGAL::Polygon_mesh_processing::fair(mesh, vector);
}


void CGALSurface::surface_intersection(CGALSurface &other) {
    CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh, other.get_mesh(), mesh);
}


void CGALSurface::surface_union(CGALSurface &other) {       // Probably bad to use union name
    CGAL::Polygon_mesh_processing::corefine_and_compute_union(
            mesh,
            other.get_mesh(),
            mesh);
}


void CGALSurface::surface_difference(CGALSurface &other) {
    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(mesh, other.get_mesh(), mesh);
}


int CGALSurface::num_faces() const {
    return mesh.number_of_faces();
}


int CGALSurface::num_edges() const {
    return mesh.number_of_edges();
}


int CGALSurface::num_vertices() const {
    return mesh.number_of_vertices();
}

void CGALSurface::split_edges(const double target_edge_length) {
     CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length, mesh);
}

#endif
