//file lean_cohomology.hpp
#ifndef LEAN_COHOMOLOGY_HPP
#define LEAN_COHOMOLOGY_HPP

#include <iostream>
#include <cstdlib>
#include <utility> 
#include <fstream>
#include <vector>
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
#include "gaussian.h" 										//from Pawel
#include "mapped_file.h"
#include "strtot.hpp"

const long double PI = 3.141592653589793238L;


const std::runtime_error 							cond_bnd_error(std::string("ERROR: Conductor on domain boundary!"));

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
std::tuple<T, T, T, T, T, T>
read_element_line(const char *str, char **endptr)
{
    T t0, t1, t2, t3, t4, t5;
	
	// std::ofstream os;
	// os.open("./debug/parser_tets.dat", std::ios::out | std::ios::app );
    
	t0 = strtot<T>(str, endptr);
    t1 = strtot<T>(*endptr, endptr);
	
	if (t1 == 4) //tetrehedron
	{
		t0 = strtot<T>(*endptr, endptr);
		t0 = strtot<T>(*endptr, endptr);
		
		t0 = t1;
		t1 = strtot<T>(*endptr, endptr);
		t2 = strtot<T>(*endptr, endptr);
		t3 = strtot<T>(*endptr, endptr);
		t4 = strtot<T>(*endptr, endptr);
		t5 = strtot<T>(*endptr, endptr);
		
		
		// os << t0 << " " << t1 << " " << t2 << " " << t3 << " " << t4 << " " << t5 << std::endl;  
	}
	else if (t1 == 2) //triangle
	{
		t0 = strtot<T>(*endptr, endptr);
		t0 = strtot<T>(*endptr, endptr);
		
		t0 = t1;
		t1 = strtot<T>(*endptr, endptr);
		t2 = strtot<T>(*endptr, endptr);
		t3 = strtot<T>(*endptr, endptr);
		t4 = strtot<T>(*endptr, endptr);
		t5 = 1;
	}
	else if (t1 == 1) //physical line
	{
		t0 = t1;
		T n_of_tags = strtot<T>(*endptr, endptr);
		for (size_t k = 0; k < n_of_tags; k++)
		{
			T t_token = strtot<T>(*endptr, endptr);
			if (k==0)
			{
				t1 = t_token;
			}
		}
		
		t2 = strtot<T>(*endptr, endptr);
		t3 = strtot<T>(*endptr, endptr);
		t4 = 1;
		t5 = 1;
	}	
	else
	{
		t0 = t1 = 0;
		t2 = t3 = t4 = t5 = 1;
	}
    
	// os.close();
    return std::make_tuple(t0, t1, t2-1, t3-1, t4-1, t5-1);
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
void
sort_uniq(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}
    
} //namespace priv

using cluster_list     		= std::vector<sgnint32_t<int32_t>>; 
using h1_2d_basis 	   		= std::vector<std::vector<int32_t>>;
using thinned_currents 		= std::vector<std::vector<int32_t>>;
using volume_type 			= std::tuple<uint32_t,uint32_t,uint32_t,uint32_t>;
using surface_type 			= std::tuple<uint32_t,uint32_t,uint32_t>;
using label_surface_type 	= std::pair<surface_type,uint32_t>;
using edge_type 			= std::tuple<uint32_t,uint32_t>;
using label_edge_type 		= std::pair<edge_type,uint32_t>;
using label_node_type 		= std::pair<uint32_t,uint32_t>;
using tm_tuple 				= std::tuple< volume_type, uint32_t >;
using sm_tuple 				= std::tuple< surface_type, uint32_t >;

class lean_cohomology
{	
	public:
		lean_cohomology(std::string , std::string, uint32_t, uint32_t, bool);
		bool 									ESTT(std::vector<std::vector<double>>&);
		size_t 									volumes_size() { return volumes.size(); }
		size_t 									surfaces_size() { return surfaces.size(); }
		size_t 									edges_size() { return edges.size(); }
		size_t 									nodes_size() { return pts.size(); }
		const std::vector<sgnint32_t<int32_t>>& vtf(const int32_t& v_id) const;
	    const std::vector<sgnint32_t<int32_t>>& fte(const int32_t& f_id) const;
	    const std::vector<sgnint32_t<int32_t>>& ftv(const int32_t& f_id) const;
		const std::vector<sgnint32_t<int32_t>>& etf(const int32_t& e_id) const;
	    const std::vector<sgnint32_t<int32_t>>& etn(const int32_t& e_id) const;
		const std::vector<sgnint32_t<int32_t>>& nte(const int32_t& n_id) const;
		bool 									is_conductor(const uint32_t& );
		std::string 							print_face(const uint32_t&, const uint32_t&, bool, double, double, double);
		std::string 							print_face(const uint32_t& f, const int32_t& orient);
		std::string 							print_edge(const uint32_t, double, double, double, double, double, double);
		std::string 							print_edge(const uint32_t&, const uint32_t&, bool, double, double, double);
		std::string 							print_dual_edge(const uint32_t& , const uint32_t&, const uint32_t&, bool, double, double, double);
		std::vector<std::string> 				print_dual_face(const uint32_t&, const uint32_t&, bool, double, double, double);
		std::string 							print_face(const uint32_t, bool, double, double, double, double, double, double, double, double, double, double, double, double );
		std::vector<double> 					face_barycenter(const uint32_t& );
		std::vector<double> 					edge_barycenter(const uint32_t& );
		std::vector<double> 					vol_barycenter(const uint32_t& );
		uint32_t 								f2d,e2d,n2d;
		
	private:
		std::pair<h1_2d_basis,thinned_currents> 	HomoCoHomo;
		uint32_t 									insulator_id, conductor_id, surface_id, n_lazy;
		std::map<uint32_t,std::vector<uint32_t>> 	physical_surfaces;
		std::vector<std::vector<std::pair<uint16_t, int16_t>>> 		vect_stt_coeffs;
		std::vector<uint32_t> 						domains;
		std::vector<volume_type> 					volumes;
		std::vector<surface_type> 					surfaces;
		std::vector<edge_type> 						edges;
		std::vector<std::vector<double> > 			pts;
		/* node -> cluster of edge IDs around it */
		std::vector<cluster_list>    				_nte_list;
		std::vector<cluster_list>    				_etn_list;	
		/* edge -> cluster of face IDs around it */
		std::vector<cluster_list>        			_etf_list;
		std::vector<cluster_list>        			_fte_list;	
		/* triangle -> cluster of volume IDs around it */
		std::vector<cluster_list>  					_ftv_list;
		std::vector<cluster_list> 					_vtf_list;
		bool										lean_or_lazy;
		
		
		void MyThrow(const std::runtime_error& e)
		{
			std::cout << "Input file error: " << e.what() << std::endl;
			throw e;
		}

		bool 										read_mesh(const std::string&, std::vector<uint32_t>&, std::vector<uint32_t>&, std::vector<uint32_t>&);
		bool 										read_gmesh(const std::string&, std::vector<uint32_t>&, std::vector<uint32_t>&, std::vector<uint32_t>&);
		void 										unique(std::vector<label_edge_type>&, std::vector<uint32_t>&);
		void 										unique(std::vector<label_surface_type>&, std::vector<uint32_t>&);	
		pair<uint32_t,sgnint32_t<int32_t> > 		check_boundary(uint32_t , const std::vector<bool>& );
		std::vector<double> 						stdcross(const std::vector<double>&, const std::vector<double>&);
		double 										stddot(const std::vector<double>&, const std::vector<double>&);
		void 										set_boundary(uint32_t , const std::vector<bool>&, 
																std::vector<pair<uint32_t,sgnint32_t<int32_t>>>&, 
																std::vector<bool>&, std::vector<std::vector<std::pair<uint16_t, int16_t>>>&,
																std::vector<std::vector<double>>& );
};

class generic_two_manifold
{
	
	public:
		generic_two_manifold(const std::vector<uint32_t>&,const std::vector<uint32_t>&,const std::vector<uint32_t>&, lean_cohomology*);
		uint16_t number_of_gens() 		{ return this->Ngen; }
		std::uint32_t Nodes() 			{ return this->n; }
		std::uint32_t Edges() 			{ return this->e; }
		std::uint32_t Faces() 			{ return this->f; }
		std::map<uint32_t, std::vector<int16_t> > 	rel_abs_cohomology(const uint32_t& );
		std::pair<h1_2d_basis,thinned_currents>   	H_to_CoH(const std::vector<uint32_t>&, const std::vector<uint32_t>&, const std::vector<uint32_t>& );
		//std::vector<int> 							RetrieveCoH(const std::vector<int>& );
		
	private:
		lean_cohomology* 				vol_mesh;
		bool							debuggy=false;
		h1_2d_basis 					h1b;
		thinned_currents 				tc;
		std::vector<bool>  				p_colour, d_colour;
		std::map<uint32_t,uint32_t> 	p_parent, p_paredge, p_distance, d_parent, d_paredge, d_distance;
		std::mutex 						p1_pushback_mutex, p2_pushback_mutex;
		std::vector<uint32_t> 			physical_surface, physical_edges, remaining_edges;
		static const bool 				output_verbose = false;
		uint32_t 						f,e,n;
		uint16_t 						Ngen;
		void 							TwoSidedAlg(const int16_t&, const uint32_t&, const uint32_t&);
		void 							RetrieveGenAndTC(const int16_t&, const uint32_t&);

};

#endif
