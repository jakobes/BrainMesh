#include <pybind11/pybind11.h>
/* #include <pybind11/operators.h> */
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "Slice.h"
#include "Surface.h"
#include "Domain.h"
#include "surface_mesher.h"


namespace py = pybind11;

template <typename... Args>
using overload_cast_ = py::detail::overload_cast_impl<Args...>;

class PyAbstractMap : public AbstractMap{
public:
       using AbstractMap::AbstractMap; /* Inherit constructors */

};

typedef std::function<double(double,double,double)> Surface_implicit_function;


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

// We want to construct Polyline in Python then let polyline translate into cgal Point in c++
// 

Point_3 Wrapper_point_3(double x,double y,double z)
{
       return Point_3(x,y,z);
}

//TODO: set default arguments

PYBIND11_MODULE(SVMTK, m) {

    py::class_<Point_3>(m, "Point_3")
       .def(py::init<double,double,double>())
       .def("x", &Point_3::x)
       .def("y", &Point_3::y)
       .def("z", &Point_3::z);
       
    py::class_<AbstractMap,PyAbstractMap> abstractmap(m,"AbstractMap");

    py::class_<SubdomainMap,AbstractMap>(m, "SubdomainMap")
        .def(py::init<>())
        .def("print",  &SubdomainMap::print)
        .def("add", &SubdomainMap::add);

    py::class_<Slice,std::shared_ptr<Slice>>(m, "Slice")
        .def(py::init<>())
        .def(py::init<Slice&>())
        .def("create_mesh", &Slice::create_mesh) 
        .def("simplify", &Slice::simplify) 
        .def("keep_component", &Slice::keep_component) 
        .def("save", &Slice::save)
        .def("add_constraints",(void (Slice::*)(Slice&,bool)) &Slice::add_constraints);



    py::class_<Surface,std::shared_ptr<Surface>>(m, "Surface")
        .def(py::init<std::string &>())                                             
        .def(py::init<>())

        .def("implicit_surface",  &Surface::implicit_surface<Surface_implicit_function> )  
        .def("triangulate_hole",&Surface::triangulate_hole)
        .def("clip",&Surface::clip )
        .def("intersection", &Surface::surface_intersection)
        .def("union", &Surface::surface_union)
        .def("difference", &Surface::surface_difference)
        .def("slice", &Surface::mesh_slice)
        .def("span", &Surface::span)

        .def("fill_holes", &Surface::fill_holes)
        .def("triangulate_faces", &Surface::triangulate_faces)
        .def("isotropic_remeshing", &Surface::isotropic_remeshing)
        .def("adjust_boundary", &Surface::adjust_boundary)

        .def("smooth_laplacian", &Surface::smooth_laplacian)
        .def("smooth_taubin", &Surface::smooth_taubin)

        // Either use these two for operator overloading, or return the vertices
        //.def("inside", &Surface::inside)
        //.def("outside", &Surface::outside)
        .def("make_cube", &Surface::make_cube)
	.def("make_cone", &Surface::make_cone)
	.def("make_cylinder", &Surface::make_cylinder)
        .def("make_sphere", &Surface::make_sphere)
        .def("num_self_intersections", &Surface::num_self_intersections) 
        .def("collapse_edges", &Surface::collapse_edges)
        .def("save", &Surface::save)
        .def("split_edges", &Surface::split_edges)
        .def("extension", &Surface::cylindric_extension)
        .def("strictly_inside", &Surface::strictly_inside )
        //.def("load", &Surface::load)//
        .def("separate_narrow_gaps", &Surface::seperate_narrow_gaps)
        .def("reconstruct", &Surface::reconstruct)

        .def("num_faces", &Surface::num_faces)
        .def("num_edges", &Surface::num_edges)
        .def("num_vertices", &Surface::num_vertices);

    py::class_<Domain,std::shared_ptr<Domain>>(m, "Domain")
        .def(py::init<Surface &>())
        .def(py::init<std::vector<Surface>>())
        .def(py::init<std::vector<Surface>, AbstractMap&>())
        .def("set_parameters", &Domain::set_parameters) // std::map<std::string, double>
        .def("set_parameter", &Domain::set_parameter)

        .def("create_mesh", (void (Domain::*)()) &Domain::create_mesh) 
        .def("create_mesh", (void (Domain::*)(double)) &Domain::create_mesh)
        .def("refine_mesh", (void (Domain::*)()) &Domain::refine_mesh)
        .def("refine_mesh", (void (Domain::*)(double)) &Domain::refine_mesh)

        .def("get_boundary", &Domain::get_boundary)

        .def("lloyd", &Domain::lloyd,     py::arg("time_limit")=0, py::arg("max_iteration_number")=0, py::arg("convergence")=0.02, py::arg("freeze_bound")=0.01,py::arg("do_freeze")=true)
        .def("odt", &Domain::odt,         py::arg("time_limit")=0, py::arg("max_iteration_number")=0, py::arg("convergence")=0.02, py::arg("freeze_bound")=0.01,py::arg("do_freeze")=true)
        .def("exude", &Domain::exude,     py::arg("time_limit")=0, py::arg("sliver_bound")=0)
        .def("perturb", &Domain::perturb, py::arg("time_limit")=0, py::arg("sliver_bound")=0)

        .def("add_sharp_border_edges", (void (Domain::*)(Surface&,double)) &Domain::add_sharp_border_edges, py::arg("surface") , py::arg("threshold")=60 ) 
        //.def("add_sharp_border_edges", overload_cast_<Surface&,double>()(&Domain::add_sharp_border_edges) ) 

        .def("reset_borders", &Domain::reset_borders)
        .def("remove_subdomain", (void (Domain::*)(std::vector<int>)) &Domain::remove_subdomain)
        .def("remove_subdomain", (void (Domain::*)(int)) &Domain::remove_subdomain)
        .def("add_corners",  &Domain::add_corners) 
        .def("number_of_cells", &Domain::number_of_cells)
        // add, remove, set 
        .def("set_borders", &Domain::set_borders)
        .def("set_features", (void(Domain::*)()) &Domain::set_features) 
        .def("set_features", (void (Domain::*)(Surface&)) &Domain::set_features) 
        .def("add_feature", &Domain::add_feature) 
        .def("save", &Domain::save, py::arg("OutPath"), py::arg("save_1Dfeatures")=false); 
       




       m.def("separate_surfaces",  (void (*)(Surface&,Surface&,Surface&)) &surface_overlapp<Surface> );
       m.def("separate_surfaces",  (void (*)(Surface&,Surface&)) &surface_overlapp<Surface> );
       m.def("morphological_surface_union", &morphological_surface_union<Surface>); 



}
