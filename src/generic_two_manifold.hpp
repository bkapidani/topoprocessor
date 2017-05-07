//file generic_two_manifold.hpp

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
#include <limits>
#include "gmsh_mesh_parser.hpp"
#include "netgen_mesh_parser.hpp"
#include "mesh.hpp"
#include "geometry_ops.hpp"
#include "colormanip.h"
#include "timecounter.h"
#include "sgnint32_t.hpp"
#include "lean_cohomology.hpp"


const long double PI = 3.141592653589793238L;

// typedef boost::associative_property_map<std::map<uint32_t,std::size_t> > rank;
// typedef boost::associative_property_map<std::map<uint32_t,uint32_t> > parent;

using h1_2d_basis 	   = std::vector<std::vector<int32_t>>;
using thinned_currents = std::vector<std::vector<int32_t>>;

class generic_two_manifold
{
	
	public:
		typedef uint16_t handle_gen_index;
		typedef uint16_t hole_gen_index;
		generic_two_manifold(std::string, std::vector<uint32_t>);
		generic_two_manifold(std::map<uint32_t,std::vector<uint32_t>>, lean_cohomology*);
		std::vector<double> edge_barycenter(const uint32_t& );
		std::vector<double> face_barycenter(const uint32_t& );
		uint16_t number_of_gens() { return this->Ngen; }
		std::uint32_t Nodes() { return this->n; }
		std::uint32_t Edges() { return this->e; }
		std::uint32_t Faces() { return this->f; }
		std::map<uint32_t, std::vector<int16_t> > rel_abs_cohomology(const uint32_t& );
		std::pair<h1_2d_basis,thinned_currents>   H_to_CoH();
{
		
	private:
		lean_cohomology* vol_mesh;
		std::vector<uint32_t> dirichlet_id;
		std::vector<uint32_t> is_dirichlet;
		std::vector<uint16_t> abs_cohomology_labels;
		std::vector<uint16_t> rel_cohomology_labels;
		static const bool output_verbose = true;
		std::string print_face(const uint32_t&, const uint32_t&, bool, double, double, double );
		void find_complex_boundary();
		std::vector<bool> node_is_boundary;
		std::vector<bool> edge_is_boundary;
		uint32_t n,e,f, old_f, old_e;
		std::string print_edge(const uint32_t&, const uint32_t&, const uint32_t&, bool, double, double, double);
		std::string print_edge(const uint32_t&, const uint32_t&, bool, double, double, double);
		std::map<uint32_t, std::set<sgnint32_t<int32_t> > > fte;
		std::map<uint32_t, std::set<sgnint32_t<int32_t> > > etf;
		std::map<uint32_t, std::set<sgnint32_t<int32_t> > > etn;	
		std::map<uint32_t, std::set<sgnint32_t<int32_t> > > nte;
		std::map<uint32_t, std::vector<double> > pts;
		uint16_t Ngen, Nhandles;
};
