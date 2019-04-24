#ifndef CGALSlice_H
#define CGALSlice_H


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/IO/Triangulation_off_ostream_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/simplify.h> 
#include <CGAL/Polyline_simplification_2/Stop_below_count_ratio_threshold.h>

#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_2.h>

#include <assert.h>
#include <iterator>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <CGAL/centroid.h>
#include <CGAL/Polygon_2_algorithms.h>

// PYBIND11
#include <pybind11/numpy.h>
namespace py = pybind11;


class CGALSlice
{
    public :
        typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
        typedef Kernel::Point_2 Point_2;

        typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb;
        typedef CGAL::Delaunay_mesh_face_base_2<Kernel> Fb;
        typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
        typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds,CGAL::Exact_predicates_tag> CDT;
        typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
        typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;
        typedef Kernel::FT FT;
        typedef std::vector<Point_2> Polyline_2;
        typedef std::vector<Polyline_2> Polylines_2;

        //typedef CGAL::Polyline_simplification_2::Stop_above_cost_threshold Stop;
        typedef CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold Stop;
        typedef CGAL::Polyline_simplification_2::Squared_distance_cost Cost;

        CGALSlice(){}
        ~CGALSlice(){}

        CGALSlice(CGALSlice &slice) { constraints = slice.get_constraints(); }

        CGALSlice(const Polylines_2 &polylines) ;

        void write_STL(const std::string filename);

        bool check_constraints();

        void save(const std::string outpath);

        void add_constraints(CGALSlice &slice) { add_constraints(slice.get_constraints()); }

        void add_constraints(Polylines_2 &polylines) {
            min_sphere.add_polylines(polylines);
            constraints.insert(constraints.end(), polylines.begin(), polylines.end());
        }

        void clear_costraints() { constraints.clear(); }

        Polylines_2& get_constraints() { return constraints; }

        void keep_component(const size_t next) {
            // 0 is largest polyline
            Polyline_2 temp = constraints[next];
            constraints.clear();
            constraints.push_back(temp);
        };

        void keep_component(const std::vector< size_t > constraint_indices) {
            // 0 is largest polyline
            Polylines_2 constraints_tmp(constraint_indices.size());
            for (auto const idx: constraint_indices)
            {
                constraints_tmp.emplace_back(constraints[idx]);
            }
            constraints.clear();        // Is this necessary?
            constraints = constraints_tmp;
        };

        py::array_t< double > get_constraints_numpy(const size_t);

        void find_holes(const int min_num_edges);

        size_t num_constraints() { return constraints.size(); }

        void create_mesh(const double mesh_resolution);

        void simplify(const double stop_crit); // simplify all polylines

        int subdomain_map(const double x, const double y);

        struct Minimum_sphere
        {
            typedef CGAL::Min_sphere_of_spheres_d_traits_2< Kernel, FT > Traits;
            typedef CGAL::Min_sphere_of_spheres_d< Traits > Min_sphere;
            typedef Traits::Sphere Sphere;

            void add_polylines(const Polylines_2 &polylines)
            {
                 for (const auto &it: polylines)
                 {
                     for (const auto &pit: it)
                          S.push_back(Sphere(pit, 0.0));
                 }
             }

             double get_bounding_sphere_radius()
             {
                 Min_sphere ms(S.begin(), S.end());
                 return CGAL::to_double(ms.radius());
             }

             private:
                 std::vector< Sphere > S;
        };


    private:
        Minimum_sphere min_sphere;
        Polyline_2 seeds;
        Polylines_2 constraints;
        CDT cdt;
        double plane_qe[4]; // Store the plane  equation -
};


CGALSlice::CGALSlice(const Polylines_2 &polylines)
{
    typedef std::vector<Point_2> plist;

    constraints = polylines;
    std::sort(constraints.begin(), constraints.end(),
            [](const plist &a, const plist &b){ return a.size() > b.size(); });

    min_sphere.add_polylines(polylines);
}


bool CGALSlice::check_constraints()
{
    if (constraints[0].size() < 20)
        return false;
    for (const auto &pol: constraints)
    {
        if (pol.empty())
            return false;
    }
    return true;
}


void CGALSlice::find_holes(const int min_num_edges)
{
    if (!check_constraints())
    {
        std::cout  << "Missing constraints, terminating find holes" << std::endl;
        return;
    }
    if (constraints.size() < 2)
        return;

    Polyline_2 temp = constraints[0];
    Point_2 c2 = CGAL::centroid(temp.begin(), temp.end(), CGAL::Dimension_tag<0>());

    if (CGAL::bounded_side_2(temp.begin(), temp.end(), c2, Kernel()) == CGAL::ON_UNBOUNDED_SIDE)
        std::cout << "Bad slice" << std::endl;

    for (auto pol = std::next(constraints.begin()); pol != constraints.end(); )
    {
        const Point_2 c2 = CGAL::centroid(pol->begin(), pol->end(), CGAL::Dimension_tag< 0 >());
        if (CGAL::bounded_side_2(temp.begin(), temp.end(), c2, Kernel()) == CGAL::ON_BOUNDED_SIDE and  // inside the largest polyline 
                CGAL::bounded_side_2(pol->begin(), pol->end(), c2, Kernel()) == CGAL::ON_BOUNDED_SIDE and
                pol->size() > min_num_edges ) // inside its own polyline i.e. have enclosed area
        {
            seeds.push_back(c2);
            ++pol;
        }
        else if (CGAL::bounded_side_2(temp.begin(), temp.end(), c2, Kernel()) == CGAL::ON_UNBOUNDED_SIDE or
                pol->size() < min_num_edges)
        {
            pol = constraints.erase(pol);
        }
        else
        {
            ++pol;
        }
    }
}


void CGALSlice::write_STL(const std::string filename)
{
    std::ofstream file(filename);
    file.precision(6);
    file << "solid "<< filename << std::endl;

    for (CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit )
    {
        file << "facet normal " << 0.0 << " " << 0.0 << " " << 1.0 <<std::endl;
        file << "outer loop"<< std::endl;

        file << "\t" << "vertex " << fit->vertex(0)->point().x() << " " << fit->vertex(0)->point().y()
            << " " << 0.0 << std::endl;
        file << "\t" << "vertex " << fit->vertex(1)->point().x() << " " << fit->vertex(1)->point().y()
            << " " << 0.0 << std::endl;
        file << "\t" << "vertex " << fit->vertex(2)->point().x() << " " << fit->vertex(2)->point().y()
            << " " << 0.0 << std::endl;

        file <<"endloop" << std::endl;
        file <<"endfacet"<< std::endl;
    }
    file << "endsolid" << std::endl;

}


void CGALSlice::simplify(const double stop_crit)
{
    if (!check_constraints())
    {
        std::cout  << "Missing constraints, terminating simplify" << std::endl;
        return;
    }
    Polylines_2 temp;
    for (const auto &pol: constraints)
    {
        Polyline_2 result;
        CGAL::Polyline_simplification_2::simplify(pol.begin(), pol.end(), Cost(), Stop(stop_crit),
                std::back_inserter(result));
        if (result.size() > 2) // a polyline segment i.e. 2 points can cause segmentation dump in attempt to mesh. 
            temp.push_back(result);
    }
    constraints = temp;
}


void CGALSlice::create_mesh(const double mesh_resolution)
{
    if(!check_constraints())
    {
        std::cout << "Missing constraints, terminating create_mesh" << std::endl;
        return;
    }

    for (const auto &pol: constraints)
         cdt.insert_constraint(pol.begin(), pol.end());

    double r = min_sphere.get_bounding_sphere_radius();
    double longest_edge = r/mesh_resolution;

    Mesher mesher(cdt);
    if (!seeds.empty())
        mesher.set_seeds(seeds.begin(), seeds.end());

    mesher.set_criteria(Criteria(0.125, longest_edge), true );
    mesher.refine_mesh(); // error step_by_step_refine_mesh(); 	
    //-------------------------------------------------------------
    // Remove facets outside can lead to errors if not a simple connected closed polyline.
    //-------------------------------------------------------------

    for(CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit)
    {
        if (!fit->is_in_domain())
            cdt.delete_face(fit);
    }
}


void CGALSlice::save(const std::string outpath)
{
    if (cdt.number_of_faces() == 0)
    {
       std::cout <<"Bad slice, will not be saved"<< std::endl;
       return;
    }

    std::string extension = outpath.substr(outpath.find_last_of(".") + 1);
    std::ofstream out(outpath);
    if (extension == "off")
    {
        std::ofstream out(outpath);
        CGAL::export_triangulation_2_to_off(out,cdt);
    }
    else if (extension == "stl")
    {
        std::cout << " only off extension is functional" << std::endl;
        write_STL(outpath);
    }
}


int CGALSlice::subdomain_map(const double x, const double y)
{
    // Return the smalles constraint of which (x, y) is a member.
    if (constraints.size() == 1)
        return 1;

    const auto vertex = Point_2(x, y);

    int subdomain_id = 2;
    for (auto pol = std::next(constraints.begin()); pol != constraints.end(); )
    {
        if (CGAL::bounded_side_2(pol->begin(), pol->end(), vertex, Kernel()) != CGAL::ON_UNBOUNDED_SIDE)
            subdomain_id++;
        ++pol;
    }
    return subdomain_id;
}


py::array_t< double > CGALSlice::get_constraints_numpy(const size_t constraint_index)
{
    if (constraint_index >= num_constraints())
    {
        std::cout << "Constraint index out of bounds." << std::endl;
        assert(false);
    }

    const auto constraint_size = constraints[constraint_index].size();
    /* auto result = py::array_t< double >({constraint_size}, {2}); */

    py::array_t< double > result;
    result.resize(std::vector< ptrdiff_t >{constraint_size, 2});
    /* auto result = py::array_t< double >(constraint_size*2); */

    /* py::array_t< double > result; */
    /* result.resize({constraint_size, 2}); */

    auto array_buf = result.request();

    double *buf_ptr = static_cast< double * >(array_buf.ptr);

    for (size_t idx = 0; idx < 2*constraint_size; idx = idx + 2)
    {
        /* buf_ptr[idx] = constraints[constraint_index][idx].x(); */
        /* buf_ptr[2*idx] = constraints[constraint_index][idx].y(); */
        buf_ptr[idx] = constraints[constraint_index][idx].x();
        buf_ptr[idx + 1] = constraints[constraint_index][idx].y();
    }

    return result;
}


#endif
