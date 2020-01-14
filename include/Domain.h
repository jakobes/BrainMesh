// Copyright (C) 2018-2019 Lars Magnus Valnes 
//
// This file is part of Surface Volume Meshing Toolkit (SVM-TK).
//
// SVM-Tk is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SVM-Tk is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SVM-Tk.  If not, see <http://www.gnu.org/licenses/>.
#ifndef __DOMAIN_H


#define __DOMAIN_H



// LOCAL
#include "Surface.h" 
#include "SubdomainMap.h" 
#include "Polyhedral_vector_to_labeled_function_wrapper.h"
//#include "Labeled_mesh_domain_with_features.h"
// STD 
#include <list>
#include <fstream>
#include <memory>
#include <set> 
#include <map>

//CGAL MESH_3
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // CGALSURFACE  CGALSURFACE ? 
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_3/Detect_polylines_in_polyhedra.h>
#include <CGAL/Mesh_3/polylines_to_protect.h>
#include <CGAL/Mesh_3/Mesher_3.h>

#include <CGAL/Mesh_3/C3T3_helpers.h>

#include <CGAL/refine_mesh_3.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// CGAL BOUNDING SPHERE
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/internal/Mesh_3/Boundary_of_subdomain_of_complex_3_in_triangulation_3_to_off.h>





// CURRENT PLAN :
// TODO: Consider void remove_isolated_vertices() as stand alone function or class function
// TODO: Clean up header files

template <class C3T3,
          class Vertex_index_property_map,
          class Facet_index_property_map,
          class Facet_index_property_map_twice,
          class Cell_index_property_map>
void output_to_medit_(std::ostream& os,
                const C3T3& c3t3,
                const Vertex_index_property_map& vertex_pmap,
                const Facet_index_property_map& facet_pmap,
                const Cell_index_property_map& cell_pmap,
                const Facet_index_property_map_twice& facet_twice_pmap = Facet_index_property_map_twice(),
                const bool print_each_facet_twice = false)
{
  // Based on CGAL output_to_medit, but writes all facets to file.
  // i.e. boundary factes and internal facets 
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Facets_in_complex_iterator Facet_iterator;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Weighted_point Weighted_point;

  const Tr& tr = c3t3.triangulation();

  //-------------------------------------------------------
  // File output
  //-------------------------------------------------------

  //-------------------------------------------------------
  // Header
  //-------------------------------------------------------
  os << std::setprecision(17);

  os << "MeshVersionFormatted 1\n"
     << "Dimension 3\n";
  //-------------------------------------------------------
  // Vertices
  //-------------------------------------------------------
  os << "Vertices\n" << tr.number_of_vertices() << '\n';

  boost::unordered_map<Vertex_handle, int> V;
  int inum = 1;
  for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
       vit != tr.finite_vertices_end();
       ++vit)
  {
    V[vit] = inum++;
    Weighted_point p = tr.point(vit);
    os << CGAL::to_double(p.x()) << ' '
       << CGAL::to_double(p.y()) << ' '
       << CGAL::to_double(p.z()) << ' '
       << get(vertex_pmap, vit)
       << '\n';
  }

  //-------------------------------------------------------
  // Facets
  //-------------------------------------------------------
  int number_of_triangles=0;
  for( auto fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)  
  {
      if (c3t3.is_in_complex(fit->first) or  c3t3.is_in_complex(fit->first->neighbor(fit->second) ) ) 
      number_of_triangles++;
  }

  
  os << "Triangles\n" 
     << number_of_triangles << '\n';

  for( auto fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)
  {
    typename C3T3::Facet f = (*fit);
    
    // Apply priority among subdomains, to get consistent facet orientation per subdomain-pair interface.

      // NOTE: We mirror a facet when needed to make it consistent with No_patch_facet_pmap_first/second.
    if (f.first->subdomain_index() > f.first->neighbor(f.second)->subdomain_index())
    {
        f = tr.mirror_facet(f);
    }
    
    // Get facet vertices in CCW order.
    Vertex_handle vh1 = f.first->vertex((f.second + 1) % 4);
    Vertex_handle vh2 = f.first->vertex((f.second + 2) % 4);
    Vertex_handle vh3 = f.first->vertex((f.second + 3) % 4);
    
    // Facet orientation also depends on parity.
    if (f.second % 2 != 0)
      std::swap(vh2, vh3);
    
    if (c3t3.is_in_complex(fit->first) or  c3t3.is_in_complex(fit->first->neighbor(fit->second) ) )  {

    os << V[vh1] << ' ' << V[vh2] << ' ' << V[vh3] << ' '; 
    os << get(facet_pmap, *fit) << '\n';
    }
    // Print triangle again if needed, with opposite orientation

  }

  //-------------------------------------------------------
  // Tetrahedra
  //-------------------------------------------------------
  os << "Tetrahedra\n"
     << c3t3.number_of_cells_in_complex() << '\n';

  for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
       cit != c3t3.cells_in_complex_end() ;
       ++cit )
  {
    for (int i=0; i<4; i++)
      os << V[cit->vertex(i)] << ' ';

    os << get(cell_pmap, cit) << '\n';
  }

  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
  os << "End\n";

}
template<typename Kernel>
struct Minimum_sphere
{
   
           typedef typename CGAL::Min_sphere_of_spheres_d_traits_3<Kernel, typename Kernel::FT> MinSphereTraits;
           typedef typename CGAL::Min_sphere_of_spheres_d<MinSphereTraits> Min_sphere;
           typedef typename MinSphereTraits::Sphere Sphere;

           template< typename MeshPolyhedron_3>  
           void add_polyhedron(const MeshPolyhedron_3 &polyhedron)
           {
                for (typename MeshPolyhedron_3::Vertex_const_iterator it=polyhedron.vertices_begin();it != polyhedron.vertices_end(); ++it)
                {
                    S.push_back(Sphere(it->point(), 0.0));
                }
            } 

           double get_bounding_sphere_radius()
           {
               Min_sphere ms(S.begin(), S.end());
               return CGAL::to_double(ms.radius());
           }
           private:
                 std::vector<Sphere> S;
};




class Surface;

class Domain {
    public :

        //---------------------------------------------------------------
        typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
        typedef Kernel::Point_3 Point_3;
        typedef CGAL::Mesh_polyhedron_3<Kernel>::type Polyhedron; 
        // Used with_features for experimental testing, currently no imporved usage of this class  
        typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel, Polyhedron> Polyhedral_mesh_domain_3; 
        typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper<Polyhedral_mesh_domain_3, Kernel  > Function_wrapper; 

        typedef Function_wrapper::Function_vector Function_vector; 
        typedef CGAL::Labeled_mesh_domain_3<Kernel> Labeled_Mesh_Domain;

        // Implemented for experimental feature, working with bifurcation exmaple 
        // Need to perform test to find number of 1-D features before failure. 
        typedef CGAL::Mesh_domain_with_polyline_features_3<Labeled_Mesh_Domain> Mesh_domain; 

 
        typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

        typedef int Curve_index;
        typedef int Corner_index;

        // experimental adjustsments to class ; Corner_index and Curve_index
        typedef CGAL::Mesh_complex_3_in_triangulation_3< Tr, Corner_index, Curve_index> C3t3;

        typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;


        // ------ Auxiliary, also  experimental----------------
        typedef Tr::Finite_vertices_iterator Finite_vertices_iterator;
        typedef Tr::Locate_type Locate_type;
        typedef Tr::Facet Facet;
        typedef std::vector<Point_3>  Polyline_3;
        typedef std::vector<Polyline_3> Polylines;
        typedef Tr::Weighted_point Weighted_point;

        typedef C3t3::Subdomain_index Subdomain_index;
        typedef C3t3::Cell_handle Cell_handle;
        typedef C3t3::Vertex_handle Vertex_handle;
        typedef C3t3::Surface_patch_index Surface_patch_index;
        typedef C3t3::Cells_in_complex_iterator Cell_iterator;
        typedef C3t3::Facets_in_complex_iterator Facet_iterator;
        typedef std::map<std::string, double> Parameters;
        typedef CGAL::Triple<Cell_handle, int, int> Edge; 
        typedef std::vector<std::size_t>  Face;

        // --------Constructors---------------------
        Domain(Surface& surface);
        Domain(std::vector<Surface> surfaces);
        Domain(std::vector<Surface> surfaces, AbstractMap& map);
        // -----------------------------------------
        ~Domain() {}


        // --------CGAL mesh criteria---------
        void set_parameters(std::map<std::string, double> &new_parameters); // py::dict? 
        void set_parameter(std::string key, double value);

        // ------- Mesh Creation and Refinement
        void create_mesh(const double mesh_resolution );
        void create_mesh();
        void refine_mesh(const double mesh_resolution );
        void refine_mesh();
        // --------CGAL mesh critera settings
        void default_parameters() {
            parameters["mesh_resolution"]=64.0;
            parameters["facet_angle"]    = 25.0;
            parameters["facet_size"]     = 0.1;
            parameters["facet_distance"] =  0.1;
            parameters["cell_radius_edge_ratio"] = 3.0;
            parameters["cell_size"] = 0.1;
            parameters["edge_size"] = 0.1;
        }



        void save(std::string OutPath, bool save_1Dfeatures); 

        Polylines& get_features() {return features; }  

        void add_feature(Polyline_3 polyline) { features.push_back(polyline);} // Check 

        void set_features();



        Polylines& get_borders() {return borders; }

        void set_borders(){domain_ptr.get()->add_features(get_borders().begin(), get_borders().end());} // -> TODO: integrate or REOMVE 

        void add_borders(Polyline_3 polyline) { borders.push_back(polyline);} 
      
        void add_sharp_border_edges(Polyhedron& polyhedron, double threshold);

        void add_sharp_border_edges(Surface& surface, double threshold=60); // default values 
 
        void add_corners(Surface& surface){  Polyline_3 points  = surface.get_corners(); domain_ptr->add_corners(points.begin(),points.end() );}

        void reset_borders(){borders.clear();} // rename

        void remove_subdomain(std::vector<int> tags);

        void set_features(Surface& surface){ set_features(surface.get_features() );} 

        void set_features(Polylines& polylines){ set_features( polylines.begin() , polylines.end() );} 

        template< typename InputIterator> // Wrapping 
        void set_features(InputIterator begin, InputIterator end){domain_ptr->add_features(begin, end);} 

        std::shared_ptr<Slice> mesh_slice(double x1,double x2, double x3 ,double x4);
        std::shared_ptr<Surface> get_boundary(int tag); // # subdomainindex


        // -----Mesh optimization ->Direct wrapping of CGAL funcions
        void lloyd(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true);// default values
        void odt(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true);  // default values
        void exude( double time_limit = 0, double sliver_bound = 0 ){ CGAL::exude_mesh_3(c3t3, sliver_bound=sliver_bound, time_limit=time_limit);}  // default values
        void perturb( double time_limit=0, double sliver_bound=0){CGAL::perturb_mesh_3 ( c3t3, *domain_ptr.get(), time_limit=time_limit, sliver_bound=sliver_bound) ;} // default values


        void remove_subdomain(int tag);              // tags 
        void rebind_facets();
        int number_of_cells(){return c3t3.number_of_cells();}
        std::set<int>  number_of_subdomains();

    private :
        std::unique_ptr<Mesh_domain> domain_ptr;
        Minimum_sphere<Kernel> min_sphere; 
        Parameters parameters; 
        C3t3 c3t3;
        // ------ Experimental variables for 1-D features -----------
        Polylines borders; 
        Polylines features;

};



template<typename C3T3>
int remove_isolated_vertices(C3T3& c3t3, bool remove_domain=false)
{ 
  //FIXME : CLEANUP
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
  typedef typename C3T3::Cell_handle Cell_handle;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;

  std::map<Vertex_handle, bool> vertex_map;


  for( Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin();vit != c3t3.triangulation().finite_vertices_end();++vit)
  { 
        vertex_map[vit] = false ;  
  }

  for(Cell_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  {
    for (std::size_t i = 0; i < 4; ++i) // 
    {

        vertex_map[cit->vertex(i)] = true;
    }
  }

  int before = c3t3.triangulation().number_of_vertices() ;

  for (typename std::map<Vertex_handle, bool>::const_iterator it = vertex_map.begin();it != vertex_map.end(); ++it) // check post or pre increment
  {
    if (!it->second) 
    {
       c3t3.triangulation().remove(it->first);
    }

  }

  int after = c3t3.triangulation().number_of_vertices() ; 



  std::cout<<"Number of isolated vertices removed: "<< before - after << std::endl;
  if ((before - after) > 10 and !remove_domain)   
  {
       // FIXME : Better 
       std::cout<<"The number of removed vertices are substantial."<< std::endl;  
       std::cout<<"This can cause the mesh to be horrid."<< std::endl;  
       std::cout<<"Try isotropic remeshing or increase the mesh resolution."<< std::endl;  
  }

  return (before - after);


}

// --------------------------------
//  Constructor for Surfaces
// --------------------------------
Domain::Domain( std::vector<Surface> surfaces )
{
    Function_vector v;

    for(std::vector<Surface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       Polyhedron *polyhedron = new Polyhedron();
       sit->get_polyhedron(*polyhedron);
       min_sphere.add_polyhedron(*polyhedron);

       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(*polyhedron);
       v.push_back(polyhedral_domain);
    }
    Function_wrapper wrapper(v);
    
    domain_ptr = std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));
    default_parameters();
}

Domain::Domain( std::vector<Surface> surfaces , AbstractMap& map )
{
    Function_vector v;

    for(std::vector<Surface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       Polyhedron *polyhedron = new Polyhedron();
       sit->get_polyhedron(*polyhedron);
       min_sphere.add_polyhedron(*polyhedron);
       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(*polyhedron);
       v.push_back(polyhedral_domain);
      
    }
    Function_wrapper wrapper(v,map);

    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));
    default_parameters();
}




Domain::Domain(Surface &surface) 
{
    Function_vector v;
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);

    min_sphere.add_polyhedron(polyhedron);

    Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);

    v.push_back(polyhedral_domain);
    Function_wrapper wrapper(v);



    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox())); 
    default_parameters();

}

// --------------------------------
//  
// --------------------------------

std::shared_ptr<Surface> Domain::get_boundary(int tag=0)
{
  // Experimental 
  // TODO: CHANGE

   std::vector<Face> faces;
   Polyline_3 points;

   typedef CGAL::Hash_handles_with_or_without_timestamps Hash_fct;
   typedef boost::unordered_map<Vertex_handle, std::size_t, Hash_fct> VHmap;
   VHmap::iterator map_entry;
   VHmap vh_to_ids;
   int inum=0;
   bool is_new;
    
   Subdomain_index useless = Subdomain_index(tag);
   facets_in_complex_3_to_triangle_soup(c3t3, useless, points, faces, true, false);  

   std::shared_ptr<Surface> surf(new  Surface(points, faces)) ;
   return surf;
}

// ---------- CGAL mesh criteria --------------------
void Domain::set_parameters(std::map<std::string, double> &new_parameters)
{
    for (std::map<std::string, double>::iterator pit= new_parameters.begin(); pit!=new_parameters.end(); ++pit )
    {
        parameters[pit->first] = static_cast<double>(pit->second);
    }
}
void Domain::set_parameter(std::string key , double value )
{
   parameters[key] = value;
}


void Domain::create_mesh()
{
    std::cout << "begin_meshing" << std::endl;
 
    Mesh_criteria criteria(CGAL::parameters::facet_angle=parameters["facet_angle"],
                           CGAL::parameters::facet_size =parameters["facet_size"],
                           CGAL::parameters::facet_distance=parameters["facet_distance"],
                           CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
                           CGAL::parameters::cell_size=parameters["cell_size"] );

    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);

    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();


}

void Domain::create_mesh(const double mesh_resolution )
{

    double r = min_sphere.get_bounding_sphere_radius(); 


    const double cell_size = r/mesh_resolution;
    std::cout << "Cell size: " << cell_size << std::endl;
    Mesh_criteria criteria(CGAL::parameters::edge_size = cell_size,
                                       CGAL::parameters::facet_angle = 30.0,
                                       CGAL::parameters::facet_size = cell_size,
                                       CGAL::parameters::facet_distance = cell_size/10.0, 
                                       CGAL::parameters::cell_radius_edge_ratio = 3.0,
                                       CGAL::parameters::cell_size = cell_size);

    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);


    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();
  
}

// direct cell_size argument ?
/*void Domain::create_mesh(const double cell_size )
{

    Mesh_criteria criteria(CGAL::parameters::edge_size = cell_size,
                                       CGAL::parameters::facet_angle = 30.0,
                                       CGAL::parameters::facet_size = cell_size,
                                       CGAL::parameters::facet_distance = cell_size/10.0, 
                                       CGAL::parameters::cell_radius_edge_ratio = 3.0,
                                       CGAL::parameters::cell_size = cell_size);

    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);


    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();
  
}*/


void Domain::save(std::string OutPath,bool save_1Dfeatures)
{
    
 
     // Hack, due to update in FEniCS
    std::ofstream  medit_file(OutPath);
    // Experimental
    typedef CGAL::Mesh_3::Medit_pmap_generator<C3t3,false,false> Generator;
    typedef typename Generator::Cell_pmap Cell_pmap;
    typedef typename Generator::Facet_pmap Facet_pmap;
    typedef typename Generator::Facet_pmap_twice Facet_pmap_twice;
    typedef typename Generator::Vertex_pmap Vertex_pmap;

    //Default_vertex_pmap vertex_pmap(c3t3,cell_pmap,facet_map);

    Cell_pmap cell_pmap(c3t3);
    Facet_pmap facet_pmap(c3t3, cell_pmap);
    Facet_pmap_twice  facet_twice_pmap(c3t3,cell_pmap);
    Vertex_pmap vertex_pmap(c3t3, cell_pmap,facet_pmap);
    //output_to_medit(medit_file,c3t3, vertex_pmap, facet_pmap, cell_pmap, facet_twice_pmap , false) ;

    output_to_medit_(medit_file,c3t3, vertex_pmap, facet_pmap, cell_pmap, facet_twice_pmap , false) ;
    //c3t3.output_to_medit(medit_file,true,true);


    medit_file.close();
   
    if (save_1Dfeatures)
    {
       boost::unordered_map<Vertex_handle, int> V;
       int inum = 0;
       for( Finite_vertices_iterator vit =  c3t3.triangulation().finite_vertices_begin();vit !=  c3t3.triangulation().finite_vertices_end();++vit)
       {
          V[vit] = inum++;
       }
       std::ofstream  textfile("1D_features.txt");

       for(C3t3::Edges_in_complex_iterator eit = c3t3.edges_in_complex_begin();eit != c3t3.edges_in_complex_end(); ++eit)
       {
         Cell_handle c = eit->first;
         Vertex_handle vh1 = c->vertex(eit->second);
         Vertex_handle vh2 = c->vertex(eit->third);
         textfile << V[vh1] << " " << V[vh2]  <<" " << c3t3.curve_index(*eit) << std::endl;
       }
       textfile.close();
    }
}


void Domain::refine_mesh()
{
    Mesh_criteria criteria(CGAL::parameters::facet_angle=parameters["facet_angle"],
                           CGAL::parameters::facet_size =parameters["facet_size"],
                           CGAL::parameters::facet_distance=parameters["facet_distance"],
                           CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
                           CGAL::parameters::cell_size=parameters["cell_size"] );

   refine_mesh_3(c3t3, *domain_ptr.get(), criteria,CGAL::parameters::no_reset_c3t3());

    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();

}
void Domain::refine_mesh(const double mesh_resolution )
{
    double r = min_sphere.get_bounding_sphere_radius(); 
    const double cell_size = r/mesh_resolution;

    Mesh_criteria criteria(CGAL::parameters::edge_size = cell_size,
                                       CGAL::parameters::facet_angle = 30.0,
                                       CGAL::parameters::facet_size = cell_size,
                                       CGAL::parameters::facet_distance = cell_size/10.0, 
                                       CGAL::parameters::cell_radius_edge_ratio = 3.0,
                                       CGAL::parameters::cell_size = cell_size);

   refine_mesh_3(c3t3, *domain_ptr.get(), criteria,CGAL::parameters::no_reset_c3t3());

   if ( remove_isolated_vertices(c3t3) > 0)
       c3t3.rescan_after_load_of_triangulation();

}



void Domain::remove_subdomain(int tag) //rename
{ 

   std::vector<int> temp;
   temp.push_back(tag);
   remove_subdomain(temp);

}




void Domain::remove_subdomain(std::vector<int> tags) //rename
{
  // TODO: simplify
  int before = c3t3.number_of_cells();

  std::map<Cell_handle,int> cell_map;
  std::vector<std::tuple<Cell_handle,int,int>> rebind;

  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  { 
     cell_map[cit]=static_cast<int>( c3t3.subdomain_index(cit) );
  }

  for( auto j = tags.begin(); j!=tags.end(); ++j)
  {
    Subdomain_index temp(*j) ; 
    for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
    {
      for (std::size_t i = 0; i < 4; i++)
      { 
         if(std::find(tags.begin(), tags.end(), static_cast<int>(c3t3.subdomain_index(cit))  ) == tags.end() ) // not a cell to be removed               
         {
           if(c3t3.subdomain_index(cit->neighbor(i))==temp)
           {
               rebind.push_back(std::make_tuple(cit, i,*j)); 
           }
         }
      }
    }
  }
     
  std::cout << rebind.size() << std::endl;

  for(std::vector<int>::iterator j = tags.begin(); j!=tags.end(); ++j)
  {
      for(auto cit =  c3t3.cells_in_complex_begin(Subdomain_index(*j));cit !=  c3t3.cells_in_complex_end(); ++cit)
      {  
          c3t3.remove_from_complex(cit); 
      }     
  }
  std::cout << cell_map.size() << std::endl;


  c3t3.rescan_after_load_of_triangulation(); 
  remove_isolated_vertices(c3t3,true);

  // rebind all, check if needed
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  {
      for(int i=0; i<4; ++i)
      {
        Cell_handle cn = cit;
        Cell_handle cm = cit->neighbor(i);
        Subdomain_index ci = c3t3.subdomain_index(cn);
        Subdomain_index cj = c3t3.subdomain_index(cm);
       
        //Subdomain_index cix = Subdomain_index(cell_map[cn]); 
        //Subdomain_index cjx = Subdomain_index(cell_map[cm]); 
        if( ci!=cj )
        {   
              if ( ci > cj ) 
              {
              c3t3.remove_from_complex(cit,i);
              c3t3.add_to_complex(cit, i, Surface_patch_index(ci,cj) ); 
              }   
              else 
              {
              c3t3.remove_from_complex(cit,i);
              c3t3.add_to_complex(cit, i, Surface_patch_index(cj,ci) ); 
              }
        }
      }
  }
  // rebind boundary
  for ( auto cit = rebind.begin(); cit!=rebind.end(); ++cit)
  {
      Cell_handle cn = std::get<0>(*cit);
      int s =std::get<1>(*cit);
      Cell_handle cm=  cn->neighbor(s);
      Subdomain_index ck = Subdomain_index(std::get<2>(*cit));
      //Subdomain_index cjx = Subdomain_index(cell_map[cm]); 

      //  c3t3.subdomain_index(cm);
             
      c3t3.remove_from_complex(cn,s);
      c3t3.add_to_complex(cn,s,Surface_patch_index(ck,0) );
  }

  int after = c3t3.number_of_cells();

  std::cout << "Number of removed subdomain cells : " << (before -after) << std::endl;
  c3t3.rescan_after_load_of_triangulation();
}



//----------------------   Have overloaded functions or just one????????????????
void Domain::add_sharp_border_edges(Polyhedron& polyhedron, double threshold) // no need to expose ?  
{ 

   Polylines polylinput; 
   typedef boost::property_map<Polyhedron, CGAL::edge_is_feature_t>::type EIF_map;
   EIF_map eif = get(CGAL::edge_is_feature, polyhedron);

   CGAL::Polygon_mesh_processing::detect_sharp_edges(polyhedron,threshold, eif); 
   for( Polyhedron::Edge_iterator he = polyhedron.edges_begin(); he != polyhedron.edges_end() ; ++he)
   {
      if(he->is_feature_edge() ) 
      {
         Polyline_3 polyline;
         polyline.push_back(he->vertex()->point());
         polyline.push_back(he->opposite()->vertex()->point());     
         polylinput.push_back(polyline);
      }    
  }   
  // need to protect internal points of sharp surface

  polylines_to_protect(this->borders, polylinput.begin(),  polylinput.end() ); // borders 
}

void Domain::add_sharp_border_edges(Surface& surface, double threshold) 
{ 
  Polyhedron polyhedron;
  surface.get_polyhedron(polyhedron);
  add_sharp_border_edges(polyhedron, threshold);

}

std::set<int>  Domain::number_of_subdomains()
{
   // TODO: simplify
   std::set<int> sd_indices;
   for(Cell_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
   {
        sd_indices.insert(static_cast<int>(c3t3.subdomain_index(cit)) );
   }
   return sd_indices;
}

std::shared_ptr<Slice> Domain::mesh_slice(double x1,double x2, double x3 ,double x4) 
{ 
   //Experimental 
   Polyline_3 points;
   std::vector<Face> faces;

   std::shared_ptr<Slice > slice(new  Slice()) ;

   for( auto index =  number_of_subdomains().begin() ; index != number_of_subdomains().end() ;++index) 
   {
        points.clear();
        faces.clear();
        facets_in_complex_3_to_triangle_soup(c3t3,  Subdomain_index(*index) , points,faces);
        Surface* surf  = new Surface(points,faces); // implemented for this purpose
        std::shared_ptr<Slice > temp = surf->mesh_slice(x1,x2,x3 ,x4);     // new object ?? 

        //temp operations 1
        slice->keep_component(0);
        slice->add_constraints(*temp.get());                     // add tag to constraint?         
   }
   return slice;
}

void Domain::rebind_facets()
{




}
void Domain::set_features()
{ 
  //Polylines polylines;
  //polylines_to_protect(polylines, features.begin() , features.end()  );
  set_features(features.begin() , features.end() );
}


void Domain::lloyd(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze )
{CGAL::lloyd_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 

void Domain::odt(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze) 
{CGAL::odt_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 


#endif
