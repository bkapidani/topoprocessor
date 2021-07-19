//file lean_cohomology.hpp
#ifndef LEAN_COHOMOLOGY_HPP
#define LEAN_COHOMOLOGY_HPP

#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <utility> 
#include <fstream>
#include <vector>
#include <array>
#include <deque>
#include <cstring>
#include <string>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <map>
#include <set>
#include <limits>
#include <sstream>
#include <future>
#include <chrono>
#include <mutex>
#include "timecounter.h"
#include "sgnint32_t.hpp"
#include "gaussian.h"      //from Pawel
#include "mapped_file.h"
#include "strtot.hpp"

const long double PI = 3.141592653589793238L;

namespace parser {

template<typename T>
std::tuple<T, T, T>
read_gmsh_point_line(const char *str, char **endptr)
{
    T t0, t1, t2, t3;
    
   t0 = strtot<T>(str, endptr); 
   t1 = strtot<T>(*endptr, endptr);
   t2 = strtot<T>(*endptr, endptr);
   t3 = strtot<T>(*endptr, endptr);
    
    //return std::make_tuple(t1/1000.0, t2/1000.0, t3/1000.0);
    return std::make_tuple(t1, t2, t3);
}

template<typename T>
std::tuple<T, T, T>
read_point_line(const char *str, char **endptr)
{
    T t1, t2, t3;
    
    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    
    //return std::make_tuple(t1/1000.0, t2/1000.0, t3/1000.0);
    return std::make_tuple(t1, t2, t3);
}

template<typename T>
std::tuple<T, T, T, T, T>
read_tetrahedron_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4, t5;
    
    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    t5 = strtot<T>(*endptr, endptr);
    
    return std::make_tuple(t1, t2-1, t3-1, t4-1, t5-1);
}

template<typename T>
std::tuple<T, T, T, T, T, T, T, T, T>
read_hex_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4, t5, t6, t7, t8, t9;
    
    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    t5 = strtot<T>(*endptr, endptr);
    t6 = strtot<T>(*endptr, endptr);
    t7 = strtot<T>(*endptr, endptr);
    t8 = strtot<T>(*endptr, endptr);
    t9 = strtot<T>(*endptr, endptr);
    
    return std::make_tuple(t1, t2-1, t3-1, t4-1, t5-1,
                           t6-1, t7-1, t8-1, t9-1);
}

template<typename T>
std::tuple<T, T, T, T, T, T, T, T, T, T>
read_element_line(const char *str, char **endptr)
{
   T t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
   t2 = t3 = t4 = t5 = t6 = t7 = t8 = t9 = T(1);
   // std::ofstream os;
   // os.open("./debug/parser_tets.dat", std::ios::out | std::ios::app );
    
   t0 = strtot<T>(str, endptr); // useless element index
   t1 = strtot<T>(*endptr, endptr); // type of element
      
   t0 = t1;
      
   T n_of_tags = strtot<T>(*endptr, endptr);
   for (size_t k = 0; k < n_of_tags; ++k)
   {
      T t_token = strtot<T>(*endptr, endptr);
      if (k==0)
      {
         t1 = t_token;
      }
   }
   
   if (t0 == 4 || t0 == 3) //tetrehedron or quadrilateral
   {
      t2 = strtot<T>(*endptr, endptr);
      t3 = strtot<T>(*endptr, endptr);
      t4 = strtot<T>(*endptr, endptr);
      t5 = strtot<T>(*endptr, endptr);
      
      // os << t0 << " " << t1 << " " << t2 << " " << t3 << " " << t4 << " " << t5 << std::endl;  
   }
   else if (t0 == 5) //hexahedron
   {
      t2 = strtot<T>(*endptr, endptr);
      t3 = strtot<T>(*endptr, endptr);
      t4 = strtot<T>(*endptr, endptr);
      t5 = strtot<T>(*endptr, endptr);
      t6 = strtot<T>(*endptr, endptr);
      t7 = strtot<T>(*endptr, endptr);
      t8 = strtot<T>(*endptr, endptr);
      t9 = strtot<T>(*endptr, endptr);
       
   }
   else if (t0 == 2) //triangle
   {
      t2 = strtot<T>(*endptr, endptr);
      t3 = strtot<T>(*endptr, endptr);
      t4 = strtot<T>(*endptr, endptr);
   }
   else if (t0 == 1) //physical line
   {
      t2 = strtot<T>(*endptr, endptr);
      t3 = strtot<T>(*endptr, endptr);
   }   
   else if (t0 == 15) //vertex
   {
      t2 = strtot<T>(*endptr, endptr);
   }
   else
      throw std::runtime_error(std::string("Unknown type of element in GMSH file!"));

   // os.close();
    return std::make_tuple(t0, t1, t2-1, t3-1, t4-1, t5-1,
                                   t6-1, t7-1, t8-1, t9-1);
}

template<typename T>
std::tuple<T, T, T, T>
read_triangle_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4;
    
    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    
    return std::make_tuple(t1, t2-1, t3-1, t4-1);
}

template<typename T>
std::tuple<T, T, T, T, T>
read_quad_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4, t5;
    
    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    t5 = strtot<T>(*endptr, endptr);
    
    return std::make_tuple(t1, t2-1, t3-1, t4-1, t5-1);
}

template<typename T>
void
sort_uniq(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}
    
} //namespace priv

using cluster_list         = std::vector<sgnint32_t<int32_t> >; 
using h1_2d_basis          = std::vector<std::vector<int32_t> >;
using thinned_currents     = std::vector<std::vector<int32_t> >;
using volume_type          = std::tuple<uint32_t,uint32_t,uint32_t,uint32_t,
                                        uint32_t,uint32_t,uint32_t,uint32_t>;
//using surface_type         = std::tuple<uint32_t,uint32_t,uint32_t,uint32_t>;
//using edge_type            = std::tuple<uint32_t,uint32_t>;
using surface_type         = std::array<uint32_t,4>;
using edge_type            = std::array<uint32_t,2>;
using label_surface_type   = std::pair<surface_type,uint32_t>;
using label_edge_type      = std::pair<edge_type,uint32_t>;
using label_node_type      = std::pair<uint32_t,uint32_t>;
using tm_tuple             = std::tuple< volume_type, uint32_t >;
using sm_tuple             = std::tuple< surface_type, uint32_t >;

class generic_two_manifold;

class lean_cohomology
{   
   public:
      lean_cohomology(std::string , std::string, const char*, const char*, bool, bool, const char*);
      bool                                      ESTT(std::vector<std::vector<double> >&,
                                                     std::vector<uint32_t>&,
                                                     std::map<uint32_t,uint32_t>&,
                                                     std::map<uint32_t,uint32_t>&);
      size_t                                    volumes_size() { return volumes.size(); }
      size_t                                    surfaces_size() { return surfaces.size(); }
      size_t                                    edges_size() { return edges.size(); }
      size_t                                    nodes_size() { return pts.size(); }
      const std::vector<sgnint32_t<int32_t> >&   vtf(const int32_t& v_id) const;
      const std::vector<sgnint32_t<int32_t> >&   fte(const int32_t& f_id) const;
      const std::vector<sgnint32_t<int32_t> >&   ftv(const int32_t& f_id) const;
      const std::vector<sgnint32_t<int32_t> >&   etf(const int32_t& e_id) const;
      const std::vector<sgnint32_t<int32_t> >&   etn(const int32_t& e_id) const;
      const std::vector<sgnint32_t<int32_t> >&   nte(const int32_t& n_id) const;
      bool                                      is_conductor(const uint32_t& );
      std::string                               print_face(const uint32_t&, const uint32_t&, bool, double, double, double);
      std::string                               print_face(const uint32_t& f, const int32_t& orient);
      std::string                               print_edge(const uint32_t, double, double, double, double, double, double);
      std::string                               print_edge(const uint32_t&, const uint32_t&, bool, double, double, double);
      std::string                               print_dual_edge(const uint32_t& , const uint32_t&, const uint32_t&, bool, double, double, double);
      std::vector<std::string>                  print_dual_face(const uint32_t&, const uint32_t&, bool, double, double, double);
      std::string                               print_face(const uint32_t, bool, double, double, double, double, double, double,                                                                           double, double, double, double, double, double );
      std::vector<double>                       face_barycenter(const uint32_t& );
      std::vector<double>                       edge_barycenter(const uint32_t& );
      std::vector<double>                       vol_barycenter(const uint32_t& );
      std::vector<bool>                         edge_in_conductor, conductor_bool;
      uint32_t                                  f2d,e2d,n2d, tree_element_id;
      
   private:
      std::pair<h1_2d_basis,thinned_currents>                  HomoCoHomo;
      uint32_t                                                 surface_id, n_lazy;
      std::map<uint32_t,std::vector<uint32_t> >                 physical_surfaces;
      std::vector<std::vector<std::pair<uint16_t, int16_t> > >   vect_stt_coeffs;
      std::vector<uint32_t>                                    domains;   
      std::vector<volume_type>                                 volumes;
      std::vector<surface_type>                                surfaces;
      std::vector<edge_type>                                   edges;
      std::vector<std::array<double,3> >                       pts;
      /* node -> cluster of edge IDs around it */
      std::vector<cluster_list>                                _nte_list;
      std::vector<cluster_list>                                _etn_list;   
      /* edge -> cluster of face IDs around it */
      std::vector<cluster_list>                                _etf_list;
      std::vector<cluster_list>                                _fte_list;   
      /* triangle -> cluster of volume IDs around it */
      std::vector<cluster_list>                                _ftv_list;
      std::vector<cluster_list>                                _vtf_list;
      bool                                                     be_lean_not_lazy, sullivan, debuggy;
      generic_two_manifold*                                    g2m;
      
      void MyThrow(const std::runtime_error& e)
      {
         std::cout << " Input file error: " << e.what() << std::endl;
         throw e;
      }

      bool read_mesh(const std::string&, std::vector<uint32_t>&, std::vector<uint32_t>&, std::vector<uint32_t>&);
      bool read_gmesh(const std::string&, std::vector<uint32_t>&, std::vector<uint32_t>&, std::vector<uint32_t>&);
      void unique(std::vector<label_edge_type>&, std::vector<uint32_t>&);
      void unique(std::vector<label_surface_type>&, std::vector<uint32_t>&);   
      pair<uint32_t,sgnint32_t<int32_t> > check_boundary(uint32_t , const std::vector<bool>& );
      std::array<double,3> stdcross(const std::array<double,3>&, const std::array<double,3>&);
      double stddot(const std::array<double,3>&, const std::array<double,3>&);
      void set_boundary(uint32_t, 
              const std::vector<bool>&,
              std::vector<pair<uint32_t,sgnint32_t<int32_t> > >&,
              std::vector<bool>&, 
              std::vector<std::vector<std::pair<uint16_t, 
              int16_t> > >&,std::vector<std::vector<double> >& );
      
      /*Sullivan algorithm functions*/
      void MinCost(const std::vector<int32_t>&, uint32_t);
      bool ResMaxFlow(std::vector<int32_t>&, const uint32_t&, const uint32_t&, const uint32_t&, const int32_t&);
      void UpdateMinCut(std::vector<int32_t>&);
      std::vector<int32_t> flow,face_coeff_in_the_chain, capacities;
      /*stuff for linking number retrieval*/
      
      void RetrieveLoop(std::map<uint32_t,uint32_t>&,std::map<uint32_t,uint32_t>&,
                        const uint32_t& ee, std::vector<std::array<double,3> >&);
      double LinkingNumber(const std::vector<std::array<double,3> >&, const std::vector<std::array<double,3> >&);
      std::vector<double> LinkingNumber(const std::vector<std::array<double,3> >&);
      double SolidAngleQuadrilateral(const std::array<double,3>& a, 
                      const std::array<double,3>& b, 
                      const std::array<double,3>& c,
                      const std::array<double,3>& d);
      double SolidAngleTriangle(const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c);
};

class generic_two_manifold
{
   friend class lean_cohomology;
   public:
      generic_two_manifold(lean_cohomology*);
      uint16_t number_of_gens()     { return this->Ngen; }
      std::uint32_t Nodes()         { return this->n; }
      std::uint32_t Edges()         { return this->e; }
      std::uint32_t Faces()         { return this->f; }
      std::map<uint32_t, std::vector<int16_t> >    rel_abs_cohomology(const uint32_t& );
      std::pair<h1_2d_basis,thinned_currents>      H_to_CoH(const std::vector<uint32_t>&, 
                                                            const std::vector<uint32_t>&, 
                                                            const std::vector<uint32_t>&,
                                                            std::vector<uint32_t>&,
                                                            std::map<uint32_t,uint32_t>&,
                                                            std::map<uint32_t,uint32_t>&);
      //std::vector<int>                      RetrieveCoH(const std::vector<int>& );
      
   private:
      lean_cohomology*             vol_mesh;
      bool                         debuggy=false;
      std::mutex                   write_thinned_current;
      h1_2d_basis                  h1b;
      thinned_currents             tc;
      std::vector<bool>            p_colour, d_colour;
      std::map<uint32_t,uint32_t>  p_parent, p_paredge, p_distance, d_parent, d_paredge, d_distance;
      std::mutex                   p1_pushback_mutex, p2_pushback_mutex;
      std::vector<uint32_t>        physical_surface, physical_edges, remaining_edges;
      uint32_t                     f,e,n;
      uint16_t                     Ngen;
      void                         TwoSidedAlg(const int16_t&, const uint32_t&, const uint32_t&);
      void                         RetrieveGenAndTC(const int16_t&, const uint32_t&);
};

#endif
