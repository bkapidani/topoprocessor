//file lean_cohomology.cpp
#include "lean_cohomology.hpp"

lean_cohomology :: lean_cohomology(std::string meshfile, uint32_t conductor_id, uint32_t insulator_id)
	: insulator_id(insulator_id), conductor_id(conductor_id)
{
	timecounter t_read, t_dk, t_hdk, t_estt, t_gauss, t_lean;
	
	std::vector<uint32_t> intersurface, physical_nodes, physical_edges;
	f2d=e2d=n2d=0;	
	t_read.tic();
    read_mesh(meshfile,intersurface,physical_edges,physical_nodes);
	t_read.toc();
	std::cout << "Loading complex took: " << t_read << " s" << std::endl;

	t_dk.tic();
	generic_two_manifold cond_ins_interface(intersurface,physical_edges,physical_nodes, this);
	t_dk.toc();
	std::cout << "Extracting conductor boundary took: " << t_dk << " s" << std::endl;	

	t_lean.tic();
	t_hdk.tic();
	HomoCoHomo = cond_ins_interface.H_to_CoH(intersurface,physical_edges, physical_nodes);
	n_lazy = cond_ins_interface.number_of_gens();
	t_hdk.toc();
	std::cout << "Computing h^1(dk) and thinned currents took: " << t_hdk << " s" << std::endl;		
	
	t_estt.tic();
	std::vector<int> vg(n_lazy,0);
	std::vector<std::vector<int>> gaussmat(n_lazy,vg);
	if (!ESTT(gaussmat))
		throw std::invalid_argument("Failed to complete ESTT routine!");
	
	t_estt.toc();
	std::cout << "Computing cuts took: " << t_estt << " s" << std::endl;		

	// arma::umat locations;
	// locations.resize(2, gmat.size());
	
	// arma::Col<double> values;
	// values.resize(gmat.size());
	// values.zeros();
	
	
	// for (uint32_t i = 0; i < gmat.size(); i++)
	// {	
		// locations(0,i) = std::get<0>(gmat.at(i));
		// locations(1,i) = std::get<1>(gmat.at(i));
		// values(i)      = std::get<2>(gmat.at(i));
	// }
	
	// gmat.clear();
	// arma::SpMat<double> gaussmat(true, locations, values, n_lazy, n_lazy);
	
	t_gauss.tic();
	// gmatrix<int64_t> gm(gaussmat);
    gmatrix<int> gm(gaussmat);
	// std::ofstream matrix_out;
	// matrix_out.open("./output/matrix.txt");
	// matrix_out << gm << std::endl;
	// matrix_out.close();

    //gm=gm.transpose();
    // gmatrix<double> change_of_basis = gm.gauss_elimination_over_reals();
	gmatrix<int> change_of_basis = gm.gaussElimination();
    std::cerr << "Rank : " << gm.number_nonzero_rows() << std::endl;
    std::cerr << "aaa  : " << gm.number_of_zero_rows() << std::endl;
    std::cerr << "Det  : " << gm.determinant()         << std::endl;
	// std::string sf1("./xRuben/matrix_after_elimination.dat");
	// std::string sf2("./xRuben/change_of_basis.dat");
    // gm.write_to_file(sf1.c_str());
    // change_of_basis.write_to_file(sf2.c_str());
	
	std::vector<int> gen_comb(n_lazy,0);
	uint32_t n_ind_gens=0;
	
	for (uint32_t k=0; k<n_lazy; ++k)
	{
		bool row_is_zero=true;
		std::vector<int> new_gen_comb = gen_comb;
		for (uint32_t j=0; j<n_lazy;++j)
		{
			if (abs(gm.Mat(k,j)) > 1e-12)
			{
				row_is_zero=false;
				break;
			}
			else if (abs(change_of_basis.Mat(j,k)) > 1e-12)
			{
				new_gen_comb[j]+= change_of_basis.Mat(j,k);
				std::cout << new_gen_comb[j] << std::endl;
			}
		}
		
		if (row_is_zero)
		{
			// std::cout << "found true generator" << std::endl;
			n_ind_gens++;
			gen_comb=std::move(new_gen_comb);
		}
		
	}

	t_gauss.toc();
	t_lean.toc();
	
	std::cout << "Computing basis of null space of l.n. matrix took: " << t_gauss << " s" << std::endl;
	std::cout << "Lean cohomology computation took: " << t_lean << " s" << std::endl;
	
	t_lean.tic();
	std::ofstream h1_final;
	h1_final.open("./h1_final.txt");

	for (uint32_t ee=0; ee<edges_size(); ++ee)
	{
		int final_coeff = 0;
		for (auto gg : HomoCoHomo.first[ee])
		{
			if (gen_comb[abs(gg)] != 0)
			{
				if (gg > 0)
					final_coeff += gen_comb[abs(gg)];
				else if (gg < 0)
					final_coeff -= gen_comb[abs(gg)];
			}
		}
		if (final_coeff != 0)
			h1_final << abs(*_etn_list[ee].begin()) << " " 
				 << abs(*_etn_list[ee].rbegin()) << " " << final_coeff << std::endl;
		//h1_final << print_edge(abs(gg),ee,!signbit(gg),0,255,0);
	}
	
	h1_final.close();
	t_lean.toc();
	
	std::cout << "Output of the final H^1 basis to file took: " << t_lean << " s" << std::endl;
}


bool lean_cohomology :: ESTT(std::vector<std::vector<int>>& gmat)
{
	std::vector<uint32_t> 										colour(pts.size()), p_queue;	
	std::vector<bool> 											wired(edges.size(),false), in_stack(surfaces.size());
	std::vector< pair <uint32_t,sgnint32_t<int32_t> > > 		cotree_stack;
	std::vector<std::vector<std::pair<uint16_t, int16_t>>> 		vect_stt_coeffs(edges.size());
	
	// std::vector<uint32_t> 									parent(pts.size()), par_edge(pts.size());
	std::map<uint32_t,uint32_t> 								parent, par_edge;

	pair<uint32_t,sgnint32_t<int32_t> > dummy;
	int32_t adj_f, adj_v;
	uint32_t root, k, nnz, j, qtop;
	sgnint32_t<int32_t>  curr_n;

	// std::ofstream os_dummy;
	// os_dummy.open("./output/cuts.txt");
	// os_dummy.close();

	
	srand (time(NULL));
	// p_queue.push_back(rand() % pts.size());
	
	// p_queue.reserve(pts.size());
	p_queue.push_back(0);
	colour[0]++;
	k=nnz=0;
	
	while (k<p_queue.size())
	{
		qtop=p_queue[k];

		for (auto curr_e : _nte_list[qtop])
		{	
			curr_n=*(_etn_list[abs(curr_e)].begin());

			if (abs(curr_n)==qtop)
				curr_n=*(std::prev(_etn_list[abs(curr_e)].end()));

			if (colour[abs(curr_n)]==0)
			{
				p_queue.push_back(abs(curr_n));
				colour[abs(curr_n)]++;
				parent[abs(curr_n)]=qtop;
				par_edge[abs(curr_n)]=abs(curr_e);
				wired[abs(curr_e)]=true;
				nnz++;
			}
			else;
		}
		  
		colour[qtop]++;
		k++;
	}

	// cotree_stack.reserve(surfaces.size());
	
	for(j=0; j<surfaces.size(); j++)
	{
		dummy=check_boundary(j, wired);
		if (dummy.first==1)
		{
			cotree_stack.push_back(make_pair(j,dummy.second));
			in_stack[j]=true;
		}
	}

	while (nnz < wired.size())
	{
		j=0;
		
		while (j < cotree_stack.size())
		{
			uint32_t curr_e = abs(cotree_stack[j].second);
			
			if (wired[curr_e]==true)
			{
				// std::cout << "closed edge " << curr_f << std::endl;
				for (auto signed_adj_f : etf(curr_e))
				{	
					uint32_t adj_f = abs(signed_adj_f);
					if (!in_stack[adj_f])
					{
						pair<uint32_t,sgnint32_t<int32_t> > dummy=check_boundary(adj_f, wired);
						if (dummy.first==1)
						{
							cotree_stack.push_back(make_pair(adj_f,dummy.second));
							in_stack[adj_f]=true;
						}
					}
				}
				
			}
			else
			{
				wired[curr_e]=true;
				nnz++;
				set_boundary(j, wired, cotree_stack, in_stack, vect_stt_coeffs, gmat);
			}
				
			j++;
		}
		
		// TEMPORARY FAIL - IMPLEMENT LINKING NUMBER INSTEAD
		
		auto discrepanza= wired.size()-nnz;
		
		if (discrepanza>0)
		{
			std::cout << "One cycle over all edges is insufficient: " << discrepanza << " edges still open" << std::endl;
			return false;
		}
	}
	
	// cout << endl << unchecked << endl;
		
	// if (nnz < pts.size())
		// return false;
	// else

	// std::cout << "Terminates!" << std::endl;
	return true;
}

pair<uint32_t,sgnint32_t<int32_t> > lean_cohomology :: check_boundary(uint32_t j, const std::vector<bool>& wired)
{
	auto i = _fte_list[j].begin();
	uint32_t open=0;
	sgnint32_t<int32_t>  ff;
	sgnint32_t<int32_t>  k;
	
	while(  (i != _fte_list[j].end()) && open<2)
	{
		k=*i;
		if (!wired[abs(k)])
		{
			open++;
			if (open==1)
				ff=k;
		}
		
		i++;
	}
	
	return std::make_pair(open,ff);
}

void lean_cohomology :: set_boundary(uint32_t j, const std::vector<bool>& w, std::vector<pair<uint32_t,sgnint32_t<int32_t>>>& cotree_stack, std::vector<bool>& in_stack, std::vector<std::vector<std::pair<uint16_t, int16_t>>>& vect_stt_coeffs, std::vector<std::vector<int>>& gmat)
{
	sgnint32_t<int32_t>  curr_e=cotree_stack[j].second;
	uint32_t 			 curr_f=cotree_stack[j].first;
	
	// std::ofstream os;
	// os.open("./output/cuts.txt", std::ofstream::out | std::ofstream::app);

	
	std::vector<int16_t> sum(n_lazy,0);
	std::vector<int16_t> coeff(n_lazy,0);
	
	for (auto gen : HomoCoHomo.second[curr_f])		
		coeff[abs(gen)-1] += signbit(gen) ? -1 : 1;
	
	for (auto signed_adj_e : _fte_list[curr_f])
	{
		// cout << "shit happens" << endl;
		uint32_t adj_e = abs(signed_adj_e);
		int16_t edge_coeff= signed_adj_e<0 ? -1 : +1;
		
		// if (vect_stt_coeffs[adj_e].size()>0)
		for (auto f_gen : vect_stt_coeffs[adj_e])
			sum[f_gen.first-1]+=f_gen.second*edge_coeff;
			
		for (auto signed_adj_f : _etf_list[adj_e] )
		{
			uint32_t adj_f=abs(signed_adj_f);
			if (!in_stack[adj_f])
			{
				pair<uint32_t,sgnint32_t<int32_t> > dummy=check_boundary(adj_f, w);
				if (dummy.first==1)
				{
					cotree_stack.push_back(make_pair(adj_f,dummy.second));
					in_stack[adj_f]=true;
				}
			}
		}
	}

	// std::ofstream dbg;
	// dbg.open("./output/dbg.txt", std::ofstream::out | std::ofstream::app);
		
	for (uint16_t gen=0; gen < n_lazy; gen++)
	{

		// dbg << "gen: " << gen+1 << "  face: " << abs(curr_f) << "  coeff: " << coeff[gen] << "   sum: " << sum[gen];
		
		int16_t edge_coeff = curr_e<0 ? -(coeff[gen]-sum[gen]) : (coeff[gen]-sum[gen]);

		if (edge_coeff != 0)
		{
			vect_stt_coeffs[abs(curr_e)].push_back(std::make_pair(gen+1,edge_coeff));
			
			// auto dual_face_vector = print_dual_face(gen+1,abs(curr_e),edge_coeff>0,0,255,0);
			// for (auto df : dual_face_vector)
				// os << df;
			
			for (auto val : HomoCoHomo.first[abs(curr_e)])
			{
				double chain_val = signbit(val) ? -1 : 1;
				gmat[uint32_t(gen)][uint32_t(abs(val)-1)]+=chain_val*edge_coeff;
				
				// std::cout << "Gmat(" << uint32_t(gen) << "," << uint32_t(abs(val)-1) << ") += " << chain_val*edge_coeff << std::endl;
			}
		}
		
		// dbg << "   result: " << edge_coeff*curr_e.Sgn() << std::endl;
	}
	
	// dbg.close();
	// os.close();

	return;
}

void lean_cohomology :: unique(std::vector<label_surface_type>& arr, std::vector<uint32_t>& new_labels)
{
	if (!arr.size())
		throw std::invalid_argument("Array must be nonempty!");
	
	uint32_t left = 0;
	uint32_t right = arr.size()-1;

	struct {
		bool operator()(const label_surface_type& t1, const label_surface_type& t2)
		{
			return (t1.first < t2.first);
		}
	} surfcomp;
	
	std::sort(arr.begin(),arr.end(),surfcomp);	
	// std::vector<uint32_t> new_labels(arr.size(),0);

	// std::cout << std::endl;
	
	uint32_t itor = left;
	new_labels[arr[0].second]=itor;
	surfaces.push_back(arr[0].first);
	
	left++;
	while (left<=right)
	{
		if (arr[left].first == arr[left-1].first)
		{
			new_labels[arr[left].second]=itor;
		}
		else
		{
			itor++;
			new_labels[arr[left].second]=itor;
			surfaces.push_back(arr[left].first);
		}
		
		left++;
	}
	
	// arr = new_arr;
	// labels = new_labels;
	
	return;
	
}

void lean_cohomology :: unique(std::vector<label_edge_type>& arr, std::vector<uint32_t>& new_labels)
{
	if (!arr.size())
		throw std::invalid_argument("Array must be nonempty!");
	
	uint32_t left = 0;
	uint32_t right = arr.size()-1;

	struct {
		bool operator()(const label_edge_type& t1, const label_edge_type& t2)
		{
			return (t1.first < t2.first);
		}
	} edgecomp;
	
	std::sort(arr.begin(),arr.end(),edgecomp);	
	// std::vector<uint32_t> new_labels(arr.size(),0);
	
	uint32_t itor = left;
	new_labels[arr[0].second]=itor;
	edges.push_back(arr[0].first);
	
	left++;
	while (left<=right)
	{
		if (arr[left].first == arr[left-1].first)
		{
			new_labels[arr[left].second]=itor;
		}
		else
		{
			itor++;
			new_labels[arr[left].second]=itor;
			edges.push_back(arr[left].first);
		}
		
		left++;
	}
	
	// arr = new_arr;
	// labels = new_labels;
	
	return;
	
}

bool lean_cohomology :: read_mesh(const std::string& _filename, std::vector<uint32_t>& intersurface, std::vector<uint32_t>& physical_edges, std::vector<uint32_t>& physical_nodes)
{	
	timecounter tc, tctot;
	
	/* Open file */
	if (_filename.size() == 0)
	{
		std::cout << "Invalid mesh file name" << std::endl;
		return false;
	}
	
	uint32_t	lines, linecount;
	
	mapped_file mf(_filename);
	
	std::cout << " * * * Reading NETGEN format mesh * * * ";
	std::cout << std::endl;
	
	tctot.tic();
	
	/************************ Read points ************************/
	linecount = 0;
	
	const char *data = mf.mem();
	char *endptr;
	
	lines = strtot<uint32_t>(data, &endptr);
	
	// pts.reserve(lines);
	
	tc.tic();
	while (linecount < lines)
	{
		if ( (linecount%100000) == 0 )
		{
			std::cout << "Reading points: " << linecount;
			std::cout << "/" << lines << "\r";
			std::cout.flush();
		}

		auto t = parser::read_point_line<double>(endptr, &endptr);
		
		
		std::vector<double> point = { std::get<0>(t),std::get<1>(t),std::get<2>(t) };		
		pts.push_back(point);
		
		/* Do something with that point */
		
		linecount++;
	}
	tc.toc();
	
	std::cout << "Reading points: " << linecount;
	std::cout << "/" << lines << " - " << tc << " seconds" << std::endl;
	
	/************************ Read tetrahedra ************************/
	linecount = 0;
	
	lines = strtot<uint32_t>(endptr, &endptr);
	std::vector< tm_tuple > temp_tet;
	temp_tet.reserve(lines);
	
	tc.tic();
	while (linecount < lines)
	{
		if ( (linecount%100000) == 0 )
		{
			std::cout << "Reading tetrahedra: " << linecount;
			std::cout << "/" << lines << "\r";
			std::cout.flush();
		}
		
		auto t = parser::read_tetrahedron_line<uint32_t>(endptr, &endptr);
		
		//auto t = parser::read_array<uint32_t, 5>(mf.get_line());
		std::vector<uint32_t> vectet(4);
		vectet[0] = std::get<1>(t);
		vectet[1] = std::get<2>(t);
		vectet[2] = std::get<3>(t);
		vectet[3] = std::get<4>(t);
		std::sort(vectet.begin(),vectet.end());
		
		uint32_t       p0(vectet[0]);
		uint32_t       p1(vectet[1]);
		uint32_t       p2(vectet[2]);
		uint32_t       p3(vectet[3]);
		uint32_t       d(std::get<0>(t));
		
		auto tuple = std::make_tuple(volume_type(p0, p1, p2, p3), d);
		temp_tet.push_back( tuple );
		
		linecount++;
	}
	tc.toc();
	
	std::cout << "Reading tetrahedra: " << linecount;
	std::cout << "/" << lines  << " - " << tc << " seconds" << std::endl;
	
	/************************ Sort ************************/
	std::cout << "Sorting data...";
	std::cout.flush();
	
	tc.tic();
	
	/* sort tetrahedra, make unique and move them in geometry */
		
	struct {
		bool operator()(const tm_tuple& t1, const tm_tuple& t2)
		{
			return (std::get<0>(t1) < std::get<0>(t2));
		}
	} mycomp;

	std::sort(temp_tet.begin(), temp_tet.end(), mycomp);
	
	// _vtf_list.resize(lines);
	uint32_t tot=0;	

	
	std::vector<label_surface_type> temp_tri0;
	temp_tri0.resize(4*lines);
	std::vector<int32_t> vol_signs;
	
	// vol_signs.reserve(lines);
	// volumes.reserve(lines);
	// domains.reserve(lines);
	
	for (auto tet : temp_tet)
	{
		auto t = std::get<0>(tet);
		volumes.push_back(t);
		
		uint32_t       p0(std::get<0>(t));
		uint32_t       p1(std::get<1>(t));
		uint32_t       p2(std::get<2>(t));
		uint32_t       p3(std::get<3>(t));

		temp_tri0[tot]           = std::make_pair(surface_type(p0, p1, p2),tot);
		temp_tri0[tot+lines]     = std::make_pair(surface_type(p0, p1, p3),tot+lines);
		temp_tri0[tot+2*lines]   = std::make_pair(surface_type(p0, p2, p3),tot+2*lines);
		temp_tri0[tot+3*lines]   = std::make_pair(surface_type(p1, p2, p3),tot+3*lines);
		tot++;
	
		std::vector<double> v1 { pts[p1][0]-pts[p0][0], pts[p1][1]-pts[p0][1],pts[p1][2]-pts[p0][2] };
		std::vector<double> v2 { pts[p2][0]-pts[p0][0], pts[p2][1]-pts[p0][1],pts[p2][2]-pts[p0][2] };
		std::vector<double> v3 { pts[p3][0]-pts[p0][0], pts[p3][1]-pts[p0][1],pts[p3][2]-pts[p0][2] };		

		int32_t sgn  = this->stddot(v1, stdcross(v2, v3))/double(6)>0? 1 : -1;
		vol_signs.push_back(sgn);
		domains.push_back(std::get<1>(tet));
	}
	std::vector<tm_tuple>().swap(temp_tet);
	std::vector<uint32_t> labels(4*lines);
	surfaces.reserve(4*lines);
	unique(temp_tri0, labels); //this also fills the surfaces vector
	std::vector<label_surface_type>().swap(temp_tri0);
	
	_ftv_list.resize(surfaces.size());
	// _vtf_list.reserve(volumes.size());
	
	for (uint32_t k=0; k<lines; k++)
	{
		// std::cout << labels[k] << " " << _ftv_list.size() << std::endl;
		sgnint32_t<int32_t> v1(k,-vol_signs[k]);
		sgnint32_t<int32_t> v2(k,vol_signs[k]);
		
		_ftv_list[labels[k]].push_back(v1); 
		_ftv_list[labels[k+lines]].push_back(v2); 
		_ftv_list[labels[k+2*lines]].push_back(v1);
		_ftv_list[labels[k+3*lines]].push_back(v2);		
		
		sgnint32_t<int32_t> f1(labels[k],-vol_signs[k]);
		sgnint32_t<int32_t> f2(labels[k+lines],vol_signs[k]);
		sgnint32_t<int32_t> f3(labels[k+2*lines],-vol_signs[k]);
		sgnint32_t<int32_t> f4(labels[k+3*lines],vol_signs[k]);
		
		std::vector<sgnint32_t<int32_t>> dummy(4);		
		_vtf_list.push_back(dummy);
		_vtf_list[k][0] = f1; 
		_vtf_list[k][1] = f2; 
		_vtf_list[k][2] = f3;
		_vtf_list[k][3] = f4;
	}
	
	std::vector<uint32_t>().swap(labels);
	lines = surfaces.size();
	intersurface.resize(lines,0);
	// _fte_list.resize(lines);
	tot = 0;
	std::vector<label_edge_type> temp_edge0(3*lines);
	
	for (auto t : surfaces)
	{
		uint32_t       p0(std::get<0>(t));
		uint32_t       p1(std::get<1>(t));
		uint32_t       p2(std::get<2>(t));

		temp_edge0[tot]          =  std::make_pair(edge_type(p0, p1),tot);
		temp_edge0[tot+lines]    =  std::make_pair(edge_type(p0, p2),tot+lines);
		temp_edge0[tot+2*lines]  =  std::make_pair(edge_type(p1, p2),tot+2*lines);
		tot++;
	}
	
	std::vector<uint32_t> e_labels(3*lines);
	unique(temp_edge0, e_labels);
	std::vector<label_edge_type>().swap(temp_edge0);
	_etf_list.resize(edges.size());
	// _fte_list.reserve(surfaces.size());
	// _etn_list.resize(edges.size());
	physical_edges.resize(edges_size(),0);
	
	for (uint32_t k=0; k<lines; k++)
	{
		
		sgnint32_t<int32_t> f1(k, 1);
		sgnint32_t<int32_t> f2(k,-1);
		sgnint32_t<int32_t> f3(k, 1);
		
		_etf_list[e_labels[k]].push_back(f1);
		_etf_list[e_labels[k+lines]].push_back(f2); 
		_etf_list[e_labels[k+2*lines]].push_back(f3);

		// std::cout << labels[k] << " " << labels[k+lines] << " " << labels[k+2*lines] << std::endl;
		
		sgnint32_t<int32_t> e1(e_labels[k],1);
		sgnint32_t<int32_t> e2(e_labels[k+lines],-1);
		sgnint32_t<int32_t> e3(e_labels[k+2*lines],1);
		
		std::vector<sgnint32_t<int32_t>> dummy(3);		
		_fte_list.push_back(dummy);
		
		_fte_list[k][0] = e1;
		_fte_list[k][1] = e2;
		_fte_list[k][2] = e3;
		
		auto vols = _ftv_list[k];
		
		switch (vols.size()) 
		{
			case 2: 
			{
				auto vol1= abs(*vols.begin());
				auto vol2= abs(*(std::prev(vols.end())));
				
				if ((domains[vol1] == insulator_id && domains[vol2] == conductor_id) || 
					(domains[vol1] == conductor_id && domains[vol2] == insulator_id))
				{
					intersurface[k]++;
					for (auto ee : _fte_list[k])
					{
						if (!physical_edges[abs(ee)])
							e2d++;
						physical_edges[abs(ee)]++;
					}
					f2d++;
				}			
				break;
			}
			case 1:
			{
				auto vol1= abs(*vols.begin());
				if (domains[vol1] == conductor_id)
				{
					std::cout << "Possibile?" << std::endl;
					intersurface[k]++;
					for (auto ee : _fte_list[k])
					{
						if (!physical_edges[abs(ee)])
							e2d++;
						physical_edges[abs(ee)]++;
					}
					f2d++;
				}
				break;
			}
			case 0:
			{
				throw std::invalid_argument("Conductor boundary cannot be on mesh boundary!");
				break;
			}
		}	
	}
	
	std::vector<uint32_t>().swap(e_labels);
	_nte_list.resize(pts.size());
	// _etn_list.reserve(edges.size());
	physical_nodes.resize(pts.size(),0);
	
	for (uint32_t k = 0; k < edges.size(); k++)
	{
		auto t = edges[k];
		
		uint32_t       p0(std::get<0>(t));
		uint32_t       p1(std::get<1>(t));
		
		sgnint32_t<int32_t> n1(p0,-1);
		sgnint32_t<int32_t> n2(p1, 1);
		
		std::vector<sgnint32_t<int32_t>> dummy(2);		
		_etn_list.push_back(dummy);
		_etn_list[k][0] = n1;
		_etn_list[k][1] = n2;
		
		sgnint32_t<int32_t> e1(k,-1);
		sgnint32_t<int32_t> e2(k, 1);
		
		_nte_list[p0].push_back(e1);
		_nte_list[p1].push_back(e2);
		
		if (physical_edges[k])
		{
			for (auto nn : _etn_list[k])
			{
				if (!physical_nodes[abs(nn)])
					n2d++;
				physical_nodes[abs(nn)]++;
			}
		}
	}
	
	tc.toc();
	
	std::cout << "done - " << tc << " seconds" << std::endl;
	
	/************************ Read boundary surfaces ************************/
	linecount = 0;
	auto num_of_tets=lines;
	lines = strtot<uint32_t>(endptr, &endptr);
	
	tc.tic();
	while (linecount < lines)
	{
		/*if (ifs.fail())
		{
			std::cout << "Error while reading boundary surfaces" << std::endl;
			return false;
		}*/
		
		if ( (linecount%50000) == 0 )
		{
			std::cout << "Reading triangle: " << linecount;
			std::cout << "/" << lines << "\r";
			std::cout.flush();
		}
		
		auto t = parser::read_triangle_line<uint32_t>(endptr, &endptr);

		
		uint32_t       p0( std::get<1>(t) );
		uint32_t       p1( std::get<2>(t) );
		uint32_t       p2( std::get<3>(t) );
		uint32_t       bid( std::get<0>(t) );
		
		surface_type   tri( p0, p1, p2 );
		
		if (!physical_surfaces[bid].size())
			physical_surfaces[bid].resize(surfaces.size(),0);

		auto itor = std::lower_bound(surfaces.begin(),surfaces.end(),tri);
		physical_surfaces[bid][std::distance(surfaces.begin(),itor)]++;
		
		linecount++;
	}
	tc.toc();
	
	std::cout << "Reading triangle: " << linecount;
	std::cout << "/" << lines  << " - " << tc << " seconds"  << std::endl;

	tctot.toc();
	// std::cout << cyan << "Total time spent in reading mesh: ";
	// std::cout << tctot << " seconds" << nocolor << std::endl;
	
	// for (const auto& nn : _vtf_list)
		// assert(nn.size() == 4);
	
	// for (const auto& nn : _ftv_list)
		// assert(nn.size() == 2 || nn.size() == 1);
	
	// for (const auto& nn : _fte_list)
		// assert(nn.size() == 3);
	
	// for (const auto& nn : _etf_list)
	// {
		// assert(nn.size() >= 2);
	// }
	
	// for (const auto& nn : _nte_list)
		// assert(nn.size() >= 3);
	// for (const auto& nn : _etn_list)
		// assert(nn.size() == 2);	
	
	return true;
}

/*const std::vector<sgnint32_t<int32_t>>& lean_cohomology :: face_in_conductor(sgnint32_t<int32_t>  f) const
{
	bool in_c = true;
	std::vector<sgnint32_t<int32_t> > ret;
	std::vector<sgnint32_t<int32_t> > false_ret;
	
	for ( auto vv : _ftv_list[abs(f)])
	{
		uint32_t dom=domains[abs(vv)];
		if ( dom == conductor_id )
			ret.push_back(vv);
		else
			in_c=false;
	}
	
	if (in_c)
		return ret;
	else
		return false_ret;
}*/

/*void lean_cohomology :: MinCost(std::vector<int32_t> start_s)
{
	uint32_t i, j, k;
	int e_or;
	
	k=0;
	
	// Initial feasible solutions for primal and dual
	std::vector<int32_t> fp(n_e);
	surfaces.size()low=fp;//fp=0;	
	surfaces.size()in_chain = start_s;
	// int maxflow=0;

	// for (auto pippo : _ftv_list.at(315))
		// std::cout << pippo.Val() << "  " << pippo.Sgn() << "  " << surfaces.size()in_chain.at(abs(pippo)) << "    ";

	// std::cout << std::endl;
	
	// for (auto pippo : _vtf_list.at(210))
		// std::cout << pippo.Val() << "  " << pippo.Sgn() << "  " << capacities.at(abs(pippo)) << "    ";
	
	// std::cout << std::endl;
	
	for ( uint32_t it=surfaces.size(); it<n_e; it++)
		if (surfaces.size()in_chain.at(it) != 0)
			std::cout << "Che cazzo sta succedendo alla faccia " << it <<"? " << surfaces.size()in_chain.at(it) << std::endl;
	
	// useful_variables
	uint32_t source, sink;
	int32_t step;
	bool RMF_SUCCESS;
	pair<sgnint32_t<int32_t> ,sgnint32_t<int32_t> > adj_v;
	std::vector<int32_t> new_min_cut(n_e), dummy_cut(n_e);
	
	timecounter tc;
	tc.tic();
	
	while (true)
	{
		while (k < n_e && ( (surfaces.size()in_chain.at(k)==0) || ((surfaces.size()in_chain.at(k)>0 && surfaces.size()low.at(k) == capacities.at(k)) || (surfaces.size()in_chain.at(k)<0 && surfaces.size()low.at(k)==-capacities.at(k))) ) )
			k++;
		
		if (k==n_e)
		{
			// std::cout << "Completed solution!" << std::endl;
			break;
		}
		else
		{			
			adj_v=make_pair(*(ftv(k).begin()),*std::prev((ftv(k).end())));
			
			if (adj_v.first < 0)
				std::swap(adj_v.first,adj_v.second);
			
			if (surfaces.size()in_chain.at(k)<0)
			{
				source=abs(adj_v.first);
				sink=abs(adj_v.second);
				step=1;
		    }
		    else
		    {
				source=abs(adj_v.second);
				sink=abs(adj_v.first);
				step=-1;
			}
		   
			new_min_cut=dummy_cut;
			// cout << "Running maxflow!" << endl; 
			RMF_SUCCESS=ResMaxFlow(new_min_cut, k, source, sink, step);
		   
			if (!RMF_SUCCESS)
				UpdateMinCut(new_min_cut);
			else;
		}

		k=0;
	}
		
	tc.toc();
	std::cout << tc << std::endl;

	ofstream fin_surf;
	fin_surf.open("./debug/dbg_fin_surf.dat");
	
	for ( uint32_t it=0; it<surfaces.size(); it++)
	{
		if (surfaces.size()in_chain.at(it) < 0)
			fin_surf << print_face(0,it,0,255.0, 255.0, 0);
		else if (surfaces.size()in_chain.at(it) > 0)
			fin_surf << print_face(0,it,1,255.0, 255.0, 0);
	}
	
	fin_surf.close();
	return;
}*/

/*bool lean_cohomology :: ResMaxFlow(std::vector<int32_t>& new_min_cut, const uint32_t& k, const uint32_t& s, const uint32_t& t, const int32_t& nmc_step)
{
	// Minimum cost network flow with restriction on already achieved orthogonalities

	std::vector<int> new_fp(n_e);
	new_min_cut.at(k)+=nmc_step;
	// std::map<uint32_t,uint32_t> parent;
	std::vector<uint32_t> parent(nodes.size(),nodes.size());
	std::vector<uint32_t> colour(nodes.size());
	std::vector<uint32_t> par_edge(nodes.size());
	// std::vector<uint32_t> distance(n_n);

	// parent.at(s)=0;
	colour.at(s)=1;
	par_edge.at(s)=k;
	// distance.at(s)=0;

	std::vector<uint32_t> p_queue;
	uint32_t i=0;
	p_queue.push_back(s);
	
	uint32_t qtop, curr_f, curr_v;
	sgnint32_t<int32_t>  adj_v;
	
	// if (k==1)
		// std::cout << "volumi adiacenti: " << *(surfaces.size()tv(k).begin()) << "  " << *std::prev(surfaces.size()tv(k).end()) << std::endl;
	
	while (i<p_queue.size() && colour.at(t)==0)
	{
		qtop=p_queue.at(i);
		// adj_f=edge_list(qtop,:);
		for (auto adj_f : volumes.size()tf(qtop))
		{	
			
			curr_f=abs(adj_f);
			curr_v=abs(*(surfaces.size()tv(curr_f).begin()));
			
			if (curr_v==qtop)
			{
				curr_v=abs(*std::prev(surfaces.size()tv(curr_f).end()));
				adj_v=*std::prev(surfaces.size()tv(curr_f).end());
			}
			else
				adj_v=*(surfaces.size()tv(curr_f).begin());
			
			
			// std::cout << "faccia corrente: " << curr_f << "dim_flow " << surfaces.size()low.size() << " " << capacities.size() << std::endl;
			
			if ((curr_v != t) || (qtop != s))
			{
				// std::cout << curr_f << std::endl;
				if (adj_v<0)
				{
					if (surfaces.size()in_chain.at(curr_f) >= 0 && surfaces.size()low.at(curr_f) < capacities.at(curr_f))
					{
						if (colour.at(curr_v)==0)
						{
							// dimQ=dimQ+1;
							// std::cout << "arrivato al nodo " << curr_v << " dal lato " << curr_f << std::endl;
							p_queue.push_back(curr_v);
							colour.at(curr_v)+=1;
							parent[curr_v]=qtop;
							par_edge.at(curr_v)=curr_f;
							// distance.at(curr_v)=distance.at(qtop)+1;
							new_fp.at(curr_f)=1;
						}	
					}
					
					new_min_cut.at(curr_f)+=1;
				}
				else 
				{
					if (surfaces.size()in_chain.at(curr_f) <= 0 && surfaces.size()low.at(curr_f) > -capacities.at(curr_f))
					{
						if (colour.at(curr_v)==0)
						{
							// dimQ=dimQ+1;
							// std::cout << "arrivato al nodo " << curr_v << " dal lato " << curr_f << std::endl;
							p_queue.push_back(curr_v);
							colour.at(curr_v)+=1;
							parent[curr_v]=qtop;
							par_edge.at(curr_v)=curr_f;
							// distance.at(curr_v)=distance.at(qtop)+1;
							new_fp.at(curr_f)=-1;
						}	
					}
					new_min_cut.at(curr_f)-=1;
				}
			}	
		}
		   
		colour.at(qtop)+=1;
		i++;
	}
	
	if (colour.at(t)==1)
	{
		uint32_t adj_v=t;
		
		while (parent.at(adj_v) < nodes.size())
		{
			// curr_f=FindIntersect(e_l, index_list.at(parent.at(adj_v)), index_list.at(adj_v), index_list.at(adj_v+1));
			curr_f=par_edge.at(adj_v);
			surfaces.size()low.at(curr_f)+=new_fp.at(curr_f);
			adj_v=parent.at(adj_v);
		}
		
		curr_f=k;

		surfaces.size()low.at(curr_f)=surfaces.size()low.at(curr_f)-nmc_step;
	   
		// std::cout << "Aumentato flusso attraverso la faccia " << curr_f << " partendo dal nodo " << t << " al nodo " << s << std::endl;
		// std::cout << "Il flusso alla faccia " << curr_f << " e' ora pari a " << surfaces.size()low.at(curr_f) << " a fronte di una capacita' di " << capacities.at(curr_f) << std::endl << std::endl;
		return 1;
	}
	else
	{
		// std::cout << "Non sono riuscito ad aumentare il flusso attraverso la faccia " << k << " partendo dal nodo " << t << " al nodo " << s << std::endl;
		
		// for (auto pippo : _vtf_list.at(s))
			// std::cout << abs(pippo) << " ";
		// std::cout << std::endl;
		
		// for (auto pippo : _vtf_list.at(t))
			// std::cout << abs(pippo) << " ";
		// std::cout << std::endl;
		
		// for (auto pippo : _ftv_list.at(k))
			// std::cout << abs(pippo) << " ";
		// std::cout << std::endl;		
		
		return 0;
	}
}*/

/*void lean_cohomology :: UpdateMinCut(std::vector<int>& new_min_cut)
{
	assert( surfaces.size()in_chain.size() == new_min_cut.size() );
	
	for (uint32_t i=0; i<surfaces.size()in_chain.size(); i++)
		surfaces.size()in_chain.at(i)+=new_min_cut.at(i);
	
	// std::cout << "Cambiato insieme di taglio!!" << std::endl << std::endl;

	// ofstream fin_surf;
	// fin_surf.open("./debug/dbg_fin_surf.dat");
	
	// for ( uint32_t it=0; it<surfaces.size(); it++)
	// {
		// if (surfaces.size()in_chain.at(it) < 0)
			// fin_surf << print_face(1,it,0,255.0, 255.0, 0);
		// else if (surfaces.size()in_chain.at(it) > 0)
			// fin_surf << print_face(1,it,1,255.0, 255.0, 0);
	// }
	
	// for ( uint32_t it=surfaces.size(); it<n_e; it++)
		// if (surfaces.size()in_chain.at(it) != 0)
		// {
			// for (auto pippo : _ftv_list.at(it))
			// {
				// std::cout << abs(pippo) << "-> ";
				// for (auto pluto : _vtf_list.at(abs(pippo)))
					// std::cout << abs(pluto) << " ";
				// std::cout << std::endl;
			// }
			// std::cout << std::endl;	
		// }
	// fin_surf.close();
	
	return;
}*/

/*double lean_cohomology :: LinkingNumber(const std::vector<std::vector<double> >& gamma0, const std::vector<std::vector<double> >& gamma1)
{
	// (const std::vector<std::vector<double> >& gamma0, const std::vector<std::vector<double> >& gamma1)
	double n = 0;
	
	// uint32_t i=0;
	// uint32_t j=0;
	// ofstream fid;
	// fid.open("./debug/linking_number_debug.txt");
	// while (i<gamma0.size() || j<gamma1.size())
	// {
		// if (i<gamma0.size())
			// fid << print_edge(1,gamma0.at(i)[0],gamma0.at(i)[1],gamma0.at(i)[2],gamma0.at((i+1)%gamma0.size())[0],gamma0.at((i+1)%gamma0.size())[1],gamma0.at((i+1)%gamma0.size())[2]);
		// if (j<gamma1.size())
			// fid << print_edge(2,gamma1.at(j)[0],gamma1.at(j)[1],gamma1.at(j)[2],gamma1.at((j+1)%gamma1.size())[0],gamma1.at((j+1)%gamma1.size())[1],gamma1.at((j+1)%gamma1.size())[2]);
		// i++;
		// j++;
	// }
	
	// fid.close();
	
	std::vector<double> a(3),b(3),c(3),d(3);

	for (uint32_t i = 0; i <gamma0.size(); i++)
	{
		for (uint32_t j = 0; j <gamma1.size(); j++)
		{
			a[0] = gamma1.at(j)[0]+0-gamma0.at(i)[0];
			a[1] = gamma1.at(j)[1]+0-gamma0.at(i)[1];
			a[2] = gamma1.at(j)[2]+0-gamma0.at(i)[2];
			b[0] = gamma1.at(j)[0]-gamma0.at((i+1) % gamma0.size())[0];
			b[1] = gamma1.at(j)[1]-gamma0.at((i+1) % gamma0.size())[1];
			b[2] = gamma1.at(j)[2]-gamma0.at((i+1) % gamma0.size())[2];
			c[0] = gamma1.at((j+1) % gamma1.size())[0]-gamma0.at((i+1) % gamma0.size())[0];
			c[1] = gamma1.at((j+1) % gamma1.size())[1]-gamma0.at((i+1) % gamma0.size())[1];
			c[2] = gamma1.at((j+1) % gamma1.size())[2]-gamma0.at((i+1) % gamma0.size())[2];
			d[0] = gamma1.at((j+1) % gamma1.size())[0]+0-gamma0.at(i)[0];
			d[1] = gamma1.at((j+1) % gamma1.size())[1]+0-gamma0.at(i)[1];
			d[2] = gamma1.at((j+1) % gamma1.size())[2]+0-gamma0.at(i)[2];
			n = n+ SolidAngleQuadrilateral(a, b, c, d);
		}
	}

    return n/(4*PI);
}*/

/*double lean_cohomology :: SolidAngleQuadrilateral(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d)
{
	return SolidAngleTriangle(a, b, c) + SolidAngleTriangle(c, d, a);
}*/

/*double lean_cohomology :: SolidAngleTriangle(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c)
{
	double determ= a[0]*(b[1]*c[2] - b[2]*c[1])+a[1]*(b[2]*c[0]- b[0]*c[2])+a[2]*(b[0]*c[1] - b[1]*c[0]);

	double al = sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2));
	double bl = sqrt(pow(b[0],2)+pow(b[1],2)+pow(b[2],2));
	double cl = sqrt(pow(c[0],2)+pow(c[1],2)+pow(c[2],2));

	double div = al * bl * cl + cl*(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]) + bl*(c[0]*a[0]+c[1]*a[1]+c[2]*a[2]) + al*(c[0]*b[0]+c[1]*b[1]+c[2]*b[2]);

	// % if in(0, div)
	// %     error('atan2 can not computed.')
	// % end

	double at; 
	if (div > 0)
	{
		at = atan2(determ,div);
	}
	else if (determ > 0)
	{
		at = atan2(-determ,-div) + PI;
	}
	else if (determ < 0)
	{
		at = atan2(-determ,-div) - PI;
	}
	else
		return EXIT_FAILURE;

	return 2*at;
}*/

const std::vector<sgnint32_t<int32_t>>& lean_cohomology :: vtf(const int32_t& v_id) const
{
	if ( _vtf_list.size() == 0 )
		throw std::invalid_argument("you have to populate the list first");
	else
		return _vtf_list[uint32_t(abs(v_id))];
}

bool lean_cohomology :: is_conductor(const uint32_t& v)
{
	return domains[v] == conductor_id? true : false;
}

const std::vector<sgnint32_t<int32_t>>& lean_cohomology :: fte(const int32_t& f_id) const
{
	if ( _fte_list.size() == 0 )
		throw std::invalid_argument("you have to populate the list first");
	else
		return _fte_list[uint32_t(abs(f_id))];
}

const std::vector<sgnint32_t<int32_t>>& lean_cohomology :: ftv(const int32_t& f_id) const
{
	if ( _ftv_list.size() == 0 )
		throw std::invalid_argument("you have to populate the list first");
	else
		return _ftv_list[uint32_t(abs(f_id))];
}

const std::vector<sgnint32_t<int32_t>>& lean_cohomology :: etf(const int32_t& e_id) const
{
	if ( _etf_list.size() == 0 )
		throw std::invalid_argument("you have to populate the list first");
	else
		return _etf_list[uint32_t(abs(e_id))];
}

const std::vector<sgnint32_t<int32_t>>& lean_cohomology :: etn(const int32_t& e_id) const
{
	if ( _etn_list.size() == 0 )
		throw std::invalid_argument("you have to populate the list first");
	else
		return _etn_list[uint32_t(abs(e_id))];
}

const std::vector<sgnint32_t<int32_t>>& lean_cohomology :: nte(const int32_t& n_id) const
{
	if ( _nte_list.size() == 0 )
		throw std::invalid_argument("you have to populate the list first");
	else
		return _nte_list[uint32_t(abs(n_id))];
}	

std::vector<double> lean_cohomology :: face_barycenter(const uint32_t& f)
{
    std::vector<uint32_t> nodes;
	std::vector<double> bc(3,0);
	
	for (const auto& signed_ee : _fte_list[f])
	{
		uint32_t ee = signed_ee.Val();
		
		for (const auto& signed_nn : _etn_list[ee])
		{
			uint32_t nn = signed_nn.Val();
			
			if (!std::binary_search(nodes.begin(),nodes.end(),nn))
			{
				nodes.push_back(nn);
				bc[0]+= pts[nn][0];
				bc[1]+= pts[nn][1];
				bc[2]+= pts[nn][2];
			}	
		}
		
	}
	
	bc[0]/=3;
	bc[1]/=3;
	bc[2]/=3;
	
	return bc;
}

std::vector<double> lean_cohomology :: vol_barycenter(const uint32_t& v)
{
    std::vector<uint32_t> nodes;
	std::vector<double> bc(3,0);
	
	for (const auto& signed_ff : _vtf_list[v])
	{
		uint32_t ff = signed_ff.Val();
		
		if (ff<surfaces.size())
		{
			for (const auto& signed_ee : _fte_list[ff])
			{
				uint32_t ee = signed_ee.Val();
				
				for (const auto& signed_nn : _etn_list[ee])
				{
					uint32_t nn = signed_nn.Val();
					
					if (!std::binary_search(nodes.begin(),nodes.end(),nn))
					{
						nodes.push_back(nn);
						bc[0]+= pts[nn][0];
						bc[1]+= pts[nn][1];
						bc[2]+= pts[nn][2];
					}	
				}
				
			}
		}
	}
	
	bc[0]/=4;
	bc[1]/=4;
	bc[2]/=4;
	
	return bc;
}

std::vector<double> lean_cohomology :: edge_barycenter(const uint32_t& e)
{
	std::vector<double> bc(3,0);
	
	for (const auto& signed_nn : _etn_list[e])
	{
		uint32_t nn = signed_nn.Val();
		bc[0]+= pts[nn][0];
		bc[1]+= pts[nn][1];
		bc[2]+= pts[nn][2];	
	}
	
	bc[0]/=2;
	bc[1]/=2;
	bc[2]/=2;
	
	return bc;
}

std::string lean_cohomology :: print_face(const uint32_t& label, const uint32_t& f, bool orient, double r, double g, double b)
{
	std::ostringstream fr;
	std::set<uint32_t> nodes;
	std::vector<std::vector<double> > n;
	
	for (auto ee : fte(f))
	{
		for (auto nn : etn(abs(ee)))
			nodes.insert(abs(nn));
		
		if (nodes.size() == 3)
			break;
	}
	
	fr << "103.000 " << label << " " << r << " " << g << " " << b << " 0.0 ";
	for (auto nn : nodes)
		n.push_back(pts[nn]);
	
	if (!orient)
		std::swap(n[1],n[2]);
	
	fr << n[0][0] << " " << n[0][1] << " " << n[0][2] << " " ;
	fr << n[1][0] << " " << n[1][1] << " " << n[1][2] << " " ;
	fr << n[2][0] << " " << n[2][1] << " " << n[2][2] << " " ;
		
	fr << std::endl;
	
	return fr.str();
}

std::string lean_cohomology :: print_dual_edge(const uint32_t& label, const uint32_t& v, const uint32_t& f, bool orient, double r, double g, double b)
{
	std::ostringstream fr;
	
	std::vector<double> fb =  face_barycenter(f);
	std::vector<double> vb =  vol_barycenter(v);
	
	fr << "102.100 " << label << " " << r << " " << g << " " << b << " 0.0 ";
	if (orient)
		fr << fb[0] << " " << fb[1] << " " << fb[2] << " " << vb[0] << " " << vb[1] << " " << vb[2] << std::endl;
	else
		fr << vb[0] << " " << vb[1] << " " << vb[2] << " " << fb[0] << " " << fb[1] << " " << fb[2] << std::endl;
	
	return fr.str();
}

std::vector<std::string> lean_cohomology :: print_dual_face(const uint32_t& label, const uint32_t& e, bool orient, double r, double g, double b)
{
	// std::ostringstream fr;
	std::vector<std::string> ret;
	std::vector<double> eb =  edge_barycenter(e);
	double x,y,z;
	x=y=z=0;
	int32_t new_or = orient? 1 : -1;
	
	
	for (auto n : _etn_list[e])
	{
		auto segno = n.Sgn();
		x+= new_or*segno*pts[abs(n)][0];
		y+= new_or*segno*pts[abs(n)][1];
		z+= new_or*segno*pts[abs(n)][2];
	}

	std::vector<double> v1 { x, y, z };
	
	for (auto f : _etf_list[e])
	{
		std::vector<double> fb = face_barycenter(abs(f));
		
		for (auto v : _ftv_list[abs(f)])
		{
			orient = new_or*f.Sgn()>0;
			std::vector<double> vb =  vol_barycenter(abs(v));

			std::vector<double> v2 { eb[0]-vb[0], eb[1]-vb[1],eb[2]-vb[2] };
			std::vector<double> v3 { fb[0]-vb[0], fb[1]-vb[1], fb[2]-vb[2] };
			
			if (this->stddot(v1,stdcross(v2,v3))>0)
				orient= true;
			else
				orient = false;
				
			
			ret.push_back(print_face(label,orient,vb[0],vb[1],vb[2],eb[0],eb[1],eb[2],fb[0],fb[1],fb[2],r,g,b));
		}
	}
	
	return ret;
}

std::string lean_cohomology :: print_face(const uint32_t& f, const int32_t& orient)
{
	std::ostringstream fr;
	std::set<uint32_t> nodes;
	std::vector<std::vector<double> > n;
	
	for (auto ee : fte(f))
	{
		for (auto nn : etn(abs(ee)))
			nodes.insert(abs(nn));
		
		if (nodes.size() == 3)
			break;
	}
	
	for (auto nn : nodes)
		fr << nn << " " ;
		
	fr << orient << std::endl;
	
	return fr.str();
}

std::string lean_cohomology :: print_edge(const uint32_t& label, const uint32_t& e, bool orient, double r, double g, double b)
{
	std::ostringstream fr;
	
	const auto& nn = _etn_list[e];
	uint32_t nn_b = abs(*nn.begin());
	uint32_t nn_e = abs(*std::prev(nn.end()));
 	std::vector<double> n1 =  pts[nn_b];
	std::vector<double> n2 =  pts[nn_e];
	
	fr << "102.100 " << label << " " << r << " " << g << " " << b << " 0.0 ";
	if (orient)
		fr << n1[0] << " " << n1[1] << " " << n1[2] << " " << n2[0] << " " << n2[1] << " " << n2[2] << std::endl;
	else
		fr << n2[0] << " " << n2[1] << " " << n2[2] << " " << n1[0] << " " << n1[1] << " " << n1[2] << std::endl;
	
	return fr.str();
}

std::string lean_cohomology :: print_edge(const uint32_t label, double x1, double y1, double z1, double x2, double y2, double z2)
{
	std::ostringstream fr;
	
	fr << "102.000 " << label << " " << 0 << " " << 0 << " " << 0 << " 0.0 ";
	fr << x1 << " " << y1 << " " << z1 << " " << x2 << " " << y2 << " " << z2 << std::endl;
	
	return fr.str();
}

std::string lean_cohomology :: print_face(const uint32_t label, bool orient, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double r, double g, double b)
{
	std::ostringstream fr;
	
	fr << "103.000 " << label << " " << r << " " << g << " " << b << " 0.0 ";
	if (orient)
		fr << x1 << " " << y1 << " " << z1 << " " << x2 << " " << y2 << " " << z2 << " " << x3 << " " << y3 << " " << z3 << std::endl;
	else
		fr << x1 << " " << y1 << " " << z1 << " " << x3 << " " << y3 << " " << z3 << " " << x2 << " " << y2 << " " << z2 << std::endl;
	
	return fr.str();
}

std::vector<double> lean_cohomology :: stdcross(const std::vector<double>& v1, const std::vector<double>& v2)
{
    std::vector<double> ret(3);
    
    ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
    ret[1] = v1[2]*v2[0] - v1[0]*v2[2];
    ret[2] = v1[0]*v2[1] - v1[1]*v2[0];
    
    return ret;
}

double lean_cohomology :: stddot(const std::vector<double>& v1, const std::vector<double>& v2)
{
    double acc=0;
    
    for (uint32_t i = 0; i < 3; i++)
        acc += v1[i] * v2[i];
        
    return acc;
}
