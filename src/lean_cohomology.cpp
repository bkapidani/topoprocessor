//file lean_cohomology.cpp
#include "lean_cohomology.hpp"

lean_cohomology :: lean_cohomology(
                                   std::string meshfile, 
                                   std::string mesher, 
                           const char* conductor_id, 
                           const char* insulator_id, 
                           bool be_lean_not_lazy,
                           bool sullivan,
                           const char* verbose)
   : be_lean_not_lazy(be_lean_not_lazy), sullivan(sullivan)
{
   timecounter t_read, t_dk, t_hdk, t_estt, t_gauss, t_lean;
   // this->be_lean_not_lazy = be_lean_not_lazy;
   
   std::ifstream is_c(conductor_id), is_i(insulator_id);
   
   std::string num_c_str,num_i_str;
   
   if (strcmp("verbose",verbose)==0)
      debuggy=true;
   else
      debuggy=false;
   
   while(std::getline(is_i,num_i_str))
   {
      if (num_i_str.size()>0)
      {
         auto n = std::stoi(num_i_str);
         while (conductor_bool.size()<n+1)
            conductor_bool.push_back(false);
      }
   }
   
   while(std::getline(is_c,num_c_str))
   {
      if (num_c_str.size()>0)
      {
         auto n = std::stoi(num_c_str);
         while (conductor_bool.size()<n+1)
            conductor_bool.push_back(false);
         conductor_bool[n]=true;
      }
   }
   
   if (debuggy)
   {
      std::cout << "    Domain labels: " << std::endl;
      for (uint32_t it=1; it<conductor_bool.size(); ++it)
         std::cout << "{ " << it << " : " 
                   << (conductor_bool[it] ? "conductor" : "insulator") << " } ";
      std::cout << std::endl;
   }
   
   std::vector<uint32_t> intersurface, physical_nodes, physical_edges, queue;
   std::map<uint32_t,uint32_t> parent, distance;
   f2d=e2d=n2d=0;
   t_lean.tic();

   t_read.tic();
   tree_element_id = 0;
   if (mesher == "netgen")
      read_mesh(meshfile,intersurface,physical_edges,physical_nodes);
   else if (mesher == "gmsh")
      read_gmesh(meshfile,intersurface,physical_edges,physical_nodes);
   t_read.toc();
   
   std::cout << "    Loading complex took: " << t_read << " seconds" << std::endl;

   t_hdk.tic();
   generic_two_manifold cond_ins_interface(this);
   
   HomoCoHomo = cond_ins_interface.H_to_CoH(intersurface,
                                            physical_edges,
                                            physical_nodes,
                                            queue,parent,distance);

   n_lazy = cond_ins_interface.number_of_gens();
   this->g2m = &cond_ins_interface;
   t_hdk.toc();
   std::cout << "    Computing h^1(dk) and thinned currents took: " << t_hdk << " seconds" << std::endl;      
   
   t_estt.tic();
   std::vector<double> vg(n_lazy,0);
   std::vector<std::vector<double> > gaussmat(n_lazy,vg);
   if (!ESTT(gaussmat,queue,parent,distance))
      throw std::invalid_argument("Failed to complete ESTT routine!");
   
   t_estt.toc();
   std::cout << "    Computing cuts took: " << t_estt << " seconds" << std::endl;      
   
   // std::vector<double> gen_comb(n_lazy,0);
   std::vector<std::vector<double> > gen_comb;
   if (be_lean_not_lazy)
   {
      t_gauss.tic();
      gmatrix<double> gm(gaussmat);

      //gm=gm.transpose();
      gmatrix<double> change_of_basis = gm.gauss_elimination_over_reals();
      
      if (debuggy)
      {
         std::cout << change_of_basis << std::endl;
         std::cout << gm << std::endl;
      }

      std::ofstream treeos("tree.txt", std::ofstream::out | std::ofstream::app);
      std::ofstream gmshos("tree_debug.txt", std::ofstream::out | std::ofstream::app);
            
      for (uint32_t k=0; k<gm.number_nonzero_rows(); ++k)
      {
         std::vector<double> new_gen_comb(n_lazy,0);
         for (uint32_t j=0; j<n_lazy;++j)
         {
            if (std::fabs(change_of_basis.Mat(j,k)) > 1e-6)
            {
               new_gen_comb[j] = change_of_basis.Mat(j,k);

               auto curr_e = cond_ins_interface.remaining_edges[j];
               treeos << abs(this->etn(curr_e)[0])+1 << " "
                      << abs(this->etn(curr_e)[1])+1 << std::endl;
               gmshos << ++tree_element_id << " " << 1 << " " << 1 << " " << 99 << " "
                      << abs(this->etn(curr_e)[0])+1 << " "
                      << abs(this->etn(curr_e)[1])+1 << std::endl;
            }
         }
         
         gen_comb.push_back(new_gen_comb);
      }

      treeos.close();
      gmshos.close();

      t_gauss.toc();
      std::cout << "    Computing null space of l.n. matrix took: " << t_gauss << " seconds" << std::endl;
      
      std::cout << "    (Kernel dimension: " << gm.number_nonzero_rows() << "  ";
      std::cout << "Determinant: " << gm.determinant() << ")" << std::endl;
   }
   else
   {
      for (uint32_t j=0; j<n_lazy;++j)
      {
         std::vector<double> ones_gen_comb(n_lazy,0);
         ones_gen_comb[j] = 1;
         gen_comb.push_back(ones_gen_comb);
      }
   }

   // cuts_pre_minimization.close();
   t_lean.toc();
   std::cout << "    Lean cohomology computation took: " << t_lean << " seconds" << std::endl;

   t_gauss.tic();
   std::ofstream h1_pre_minimization, h1_post_minimization, cuts_pre_minimization;
   if (debuggy)
      cuts_pre_minimization.open("./cuts_pre_minimization.txt", std::ofstream::out | std::ofstream::app);
   
   if (!sullivan)
   {
      h1_pre_minimization.open("h1.txt");
      h1_pre_minimization  << gen_comb.size() << std::endl;
   }
   else
   {
      h1_post_minimization.open("h1.txt");
      h1_post_minimization << 1 << std::endl;
   }
   
   std::vector<int32_t> start_s(edges_size(),0);
   for (uint32_t indigen=0; indigen<gen_comb.size(); ++indigen)
   {
      std::vector<std::tuple<uint32_t,uint32_t,double> > coefficients;
      
      for (uint32_t ee=0; ee<edges_size(); ++ee)
      {
         double final_coeff = 0;
         for (auto bb : vect_stt_coeffs[ee])
         {
            auto gg = bb.first;
            //~ std::cout << gg << " ";
            
            if (gen_comb[indigen][gg-1] != 0)
            {
               //~ if (gen_comb[indigen][abs(gg)-1] != 1)
                  //~ std::cout << gen_comb[indigen][abs(gg)-1] << "!" << std::endl;
               
               final_coeff += double(bb.second)*double(gen_comb[indigen][gg-1]);
            }
         }
         //~ std::cout << std::endl;
         if (std::fabs(final_coeff) > 1e-12)
         {
            start_s[ee] += std::round(final_coeff);
            if (!edge_in_conductor[ee])
            {
               if (!sullivan)
                  coefficients.push_back(std::make_tuple(abs(*_etn_list[ee].begin()),abs(*_etn_list[ee].rbegin()),final_coeff));
               if (debuggy)
               {
                  auto dual_face_vector = print_dual_face(indigen%9+1,ee,final_coeff>0,255,255,0);
                  for (auto df : dual_face_vector)
                     cuts_pre_minimization << df;
               }
            }
         }
         //~ else
            //~ start_s[ee] = 0;
      }
      
      /*Run Sullivan*/
        //~ if (sullivan)
            //~ MinCost(start_s,indigen);
        //~ else
        //~ {
            //~ std::remove("cuts_post_minimization.txt");
            //~ std::remove("h1.txt");
        //~ }
      if (!sullivan)
      {
         h1_pre_minimization << coefficients.size() << std::endl;
         for (auto cfs : coefficients)
            h1_pre_minimization << std::get<0>(cfs)+1 << "\t" << std::get<1>(cfs)+1 << "\t" << std::get<2>(cfs) << std::endl;
      }
   }
   t_gauss.toc();
   std::cout << "    Output to file took: " << t_gauss << " seconds" << std::endl;   
    
    //~ timecounter t_sullivan;
    //~ t_sullivan.tic();
    /*Run Sullivan*/
    if (sullivan)
    {
        MinCost(start_s,1);
        //~ std::remove("h1.txt");
    }
    else
    {
       std::remove("cuts_post_minimization.txt");
       //~ std::remove("h1.txt");
    }
    
    if (!debuggy)
    {
       std::remove("cuts_post_minimization.txt");
       //~ std::remove("h1.txt");
       std::remove("cuts_pre_minimization.txt");
    }
    
    h1_pre_minimization.close();
    if (debuggy)
      cuts_pre_minimization.close();
}

bool lean_cohomology :: ESTT(std::vector<std::vector<double> >& gmat,
                             std::vector<uint32_t>& queue,
                             std::map<uint32_t,uint32_t>& parent,
                             std::map<uint32_t,uint32_t>& distance)
{
   std::vector<uint32_t>                                    colour(pts.size(),0);
   std::vector<bool>                                        wired(edges.size(),false), in_stack(surfaces.size());
   std::vector< pair <uint32_t,sgnint32_t<int32_t> > >      cotree_stack;
   std::vector<std::vector<std::pair<uint16_t, int16_t> > >   vect_stt_coeffs(edges.size());
   
   // std::vector<uint32_t>                                    queue;
   // std::map<uint32_t,uint32_t>                              parent, distance;

   pair<uint32_t,sgnint32_t<int32_t> > dummy;
   int32_t adj_f, adj_v;
   uint32_t root, k, nnz, j, qtop;
   sgnint32_t<int32_t>  curr_n;

   srand (time(NULL));
   
   for (uint32_t kk = 0; kk < queue.size(); ++kk)
   {
      colour[queue[kk]]++;
      //std::cout << queue[kk] << std::endl;
   }
   k = nnz = 0;

   for (uint32_t kk = 0; kk < edges.size(); ++kk)
   {
      if (colour[abs(*(_etn_list[kk].begin()  ) ) ]>0 &&
          colour[abs(*(_etn_list[kk].rbegin() ) ) ]>0)
      {
         wired[kk] = true;
         nnz++;
      }
   }
   
   std::ofstream treeos("tree.txt", std::ofstream::out | std::ofstream::app);
   std::ofstream gmshos("tree_debug.txt", std::ofstream::out | std::ofstream::app);

   //treeos << "internal tree: " << std::endl;
   while (k<queue.size())
   {
      qtop=queue[k];

      for (auto curr_e : _nte_list[qtop])
      {   
         curr_n=*(_etn_list[abs(curr_e)].begin());

         if (abs(curr_n)==qtop)
            curr_n=*(std::prev(_etn_list[abs(curr_e)].end()));

         if (colour[abs(curr_n)]==0)
         {
            queue.push_back(abs(curr_n));
            colour[abs(curr_n)]++;
            parent[abs(curr_n)]=qtop;
            distance[abs(curr_n)] = distance[qtop]+1;
            wired[abs(curr_e)]=true;
            
            treeos << abs(this->etn(abs(curr_e))[0])+1 << " "
                   << abs(this->etn(abs(curr_e))[1])+1 << std::endl;
            gmshos << ++tree_element_id << " " << 1 << " " << 1 << " " << 99 << " "
                   << abs(this->etn(abs(curr_e))[0])+1 << " "
                   << abs(this->etn(abs(curr_e))[1])+1 << std::endl;
            
            nnz++;
         }
         else;

      }
        
      colour[qtop]++;
      k++;
   }
   
   treeos.close();
   gmshos.close();

   for(j=0; j<surfaces.size(); ++j)
   {
      dummy=check_boundary(j, wired);
      if (dummy.first==1)
      {
         cotree_stack.push_back(make_pair(j,dummy.second));
         in_stack[j]=true;
      }
   }
   
   
   j=0;
   
   while (nnz < wired.size())
   {
      while (j < cotree_stack.size())
      {
         uint32_t curr_e = abs(cotree_stack[j].second);
         
         if (wired[curr_e]==true)
         {
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
      
      auto discrepancy= wired.size()-nnz;
      
      if (discrepancy>0)
      {
         //std::cout << "    One cycle over all edges is insufficient: " 
         //          << discrepancy << " edges still open" << std::endl;
         
         uint32_t it=0;
         while (it<wired.size() && wired[it]) it++;
         
         if (it<wired.size())
         {
            //std::cout << "    Setting random cotree edge value by l.n. computations" << std::endl;
            
            std::vector<std::array<double,3> > loop;
            RetrieveLoop(parent,distance,it,loop);
            std::vector<double> new_coeffs = LinkingNumber(loop);
            
            for (uint16_t gen=0; gen < n_lazy; ++gen)
            {
               auto edge_coeff = int16_t(std::round(new_coeffs[gen]));
               if (edge_coeff != 0)
               {
                  vect_stt_coeffs[it].push_back(std::make_pair(gen+1,edge_coeff));
                  
                  if (be_lean_not_lazy)
                  {
                     for (auto val : HomoCoHomo.first[it])
                     {
                        double chain_val = signbit(val) ? -1 : 1;
                        gmat[uint32_t(gen)][uint32_t(abs(val)-1)]+=chain_val*edge_coeff;
                     }
                  }
               }
            }
            
            wired[it] = true;
            nnz++;
            cotree_stack.push_back(std::make_pair(uint32_t(0),sgnint32_t<int32_t>({it,1})));
         }
         else
         {
            std::cout << "    BUG: No cotree edges left to be set!" << std::endl;
            return false;
         }
      }
   }
   
   /*for (uint32_t l=0;l<surfaces_size();++l)
   {
      std::vector<int16_t> coeff(n_lazy+1,0);
      std::vector<int16_t> sum(n_lazy+1,0);
      for (auto tc : HomoCoHomo.second[l])
      {
         coeff[abs(tc)]+=signbit(tc) ? -1 : 1;
         //~ if (abs(coeff[abs(tc)])>1)
            //~ std::cout << coeff[abs(tc)] << std::endl;
      }
      
      for (auto ee : _fte_list[l])
         for (auto cfs : vect_stt_coeffs[abs(ee)])
            sum[cfs.first] += cfs.second*ee.Sgn();
        
      for (uint32_t i_l=1; i_l<=n_lazy; ++i_l)
         if (coeff[i_l] != sum[i_l])
            std::cout << coeff[i_l] << " != " << sum[i_l] << std::endl;
         //~ else
            //~ std::cout << "Fine!" << std::endl;
   }*/
   
   this->vect_stt_coeffs = std::move(vect_stt_coeffs);

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

void lean_cohomology :: set_boundary(uint32_t j, const std::vector<bool>& w, 
                                     std::vector<pair<uint32_t,sgnint32_t<int32_t> > >& cotree_stack, 
                                     std::vector<bool>& in_stack, 
                                     std::vector<std::vector<std::pair<uint16_t, int16_t> > >& vect_stt_coeffs, 
                                     std::vector<std::vector<double> >& gmat)
{
   sgnint32_t<int32_t>  curr_e=cotree_stack[j].second;
   uint32_t             curr_f=cotree_stack[j].first;
   uint8_t nw=0;

   
   std::vector<int16_t> sum(n_lazy,0);
   std::vector<int16_t> coeff(n_lazy,0);
   
   for (auto gen : HomoCoHomo.second[curr_f])
      coeff[abs(gen)-1] += (signbit(gen) ? -1 : 1);
   
   for (auto signed_adj_e : _fte_list[curr_f])
   {
      uint32_t adj_e = abs(signed_adj_e);
      int16_t edge_coeff= signed_adj_e<0 ? -1 : +1;
      
      if (adj_e != abs(curr_e))
      {
         nw++;
         for (auto f_gen : vect_stt_coeffs[adj_e])
            sum[f_gen.first-1]+=f_gen.second*edge_coeff;
      }
      
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
   
   //if (nw>2)
   //   std::cout << "big trouble" << std::endl;
      
   for (uint16_t gen=0; gen < n_lazy; ++gen)
   {
      int16_t edge_coeff = curr_e<0 ? -(coeff[gen]-sum[gen]) : (coeff[gen]-sum[gen]);

      if (edge_coeff != 0)
      {
         vect_stt_coeffs[abs(curr_e)].push_back(std::make_pair(gen+1,edge_coeff));
         
         if (be_lean_not_lazy)
         {
            for (auto val : HomoCoHomo.first[abs(curr_e)])
            {
               double chain_val = signbit(val) ? -1 : 1;
               gmat[uint32_t(gen)][uint32_t(abs(val)-1)]+=chain_val*edge_coeff;
            }
         }
      }
   }

   return;
}

void lean_cohomology :: RetrieveLoop(std::map<uint32_t,uint32_t>& p_parent,
                                     std::map<uint32_t,uint32_t>& p_distance,
                                     const uint32_t& ee, std::vector<std::array<double,3> >& gamma0)
{
   std::deque<std::array<double,3> > pts_deque;
   auto ff = _etn_list[ee];

   uint32_t u = abs(*(ff.begin()));
   uint32_t v = abs(*(ff.rbegin()));
      
   bool from_small_to_big=true;
   uint32_t small,big;
   
   pts_deque.push_front(pts[u]);
   pts_deque.push_back(pts[v]);
   
   while (true)
   {
      if (p_distance[u]<p_distance[v])
      {
         small = u;
         big = v;
         from_small_to_big = !from_small_to_big;
      }
      else
      {
         small = v;
         big = u;
      }

      
      if (from_small_to_big)
         pts_deque.push_front(pts[p_parent[big]]);
      else
         pts_deque.push_back(pts[p_parent[big]]);
      
      if (p_parent[big] == small)
         break;
      
      u=small;
      from_small_to_big = !from_small_to_big;
      v=p_parent[big];
   }
   
   while (!pts_deque.empty())
   {
      gamma0.push_back(pts_deque.front());
      pts_deque.pop_front();
   }
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
         auto ct1 = t1.first;
         auto ct2 = t2.first;
         std::sort(ct1.begin(),ct1.end());
         std::sort(ct2.begin(),ct2.end());
         return (ct1 < ct2);
      }
   } surfcomp;
   
   std::sort(arr.begin(),arr.end(),surfcomp);
   // std::vector<uint32_t> new_labels(arr.size(),0);

   // std::cout << std::endl;
   
   uint32_t index = left;
   new_labels[arr[0].second]=index;
   surfaces.push_back(arr[0].first);
   
   left++;
   while (left<=right)
   {
      auto f1 = arr[left].first;
      auto f2 = arr[left-1].first;
      std::sort(f1.begin(),f1.end());
      std::sort(f2.begin(),f2.end());

      if (f1 == f2)
      {
         new_labels[arr[left].second]=index;
      }
      else
      {
         index++;
         new_labels[arr[left].second]=index;
         surfaces.push_back(arr[left].first);
      }
      
      left++;
   }
   
   // arr = new_arr;
   // labels = new_labels;
   
   // std::cout << std::endl;
   // for (auto i=0;i<surfaces.size();++i)
   // {
   //    auto pippo = surfaces[i];
   //    std::sort(pippo.begin(),pippo.end());
   //    for (auto j=0; j<4;++j)
   //       std::cout << pippo[j] << " ";
   //    std::cout << std::endl;
   // }

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
         auto ct1 = t1.first;
         auto ct2 = t2.first;
         std::sort(ct1.begin(),ct1.end());
         std::sort(ct2.begin(),ct2.end());
         return (ct1 < ct2);
      }
   } edgecomp;
   
   std::sort(arr.begin(),arr.end(),edgecomp);   
   // std::vector<uint32_t> new_labels(arr.size(),0);
   
   uint32_t index = left;
   new_labels[arr[0].second]=index;
   edges.push_back(arr[0].first);
   
   left++;
   while (left<=right)
   {
      auto f1 = arr[left].first;
      auto f2 = arr[left-1].first;
      std::sort(f1.begin(),f1.end());
      std::sort(f2.begin(),f2.end());
      if (f1 == f2)
      {
         new_labels[arr[left].second]=index;
      }
      else
      {
         index++;
         new_labels[arr[left].second]=index;
         edges.push_back(arr[left].first);
      }
      
      left++;
   }
   
   // std::cout << std::endl;
   // for (auto i=0;i<edges.size();++i)
   // {
   //    auto pippo = edges[i];
   //    std::sort(pippo.begin(),pippo.end());
   //    for (auto j=0; j<2;++j)
   //       std::cout << pippo[j] << " ";
   //    std::cout << std::endl;
   // }
   
   return;
   
}

bool lean_cohomology :: read_gmesh(const std::string&     _filename, 
                                   std::vector<uint32_t>& intersurface,
                                   std::vector<uint32_t>& physical_edges,
                                   std::vector<uint32_t>& physical_nodes)
{
   
   /* Open file */
   if (_filename.size() == 0)
   {
      std::cout << "    Invalid mesh file name" << std::endl;
      return false;
   }
   
   uint32_t   lines, linecount;
   
   mapped_file mf(_filename);
   
   std::stringstream sng;
   uint32_t string_ind = _filename.size();
   while (_filename[--string_ind] != '.');
   for (uint32_t ii = 0; ii < string_ind; ++ii)
      sng << _filename[ii];
   
   sng << ".mesh";
   std::ofstream os(sng.str().c_str());

   //std::cout << sng.str().c_str() << std::endl;
   /************************ Read useless stuff ************************/

   for ( size_t i=0; i<4; i++)
      auto t = mf.get_line();

   /************************ Read points ************************/
   linecount = 0;
   
   const char *data = mf.mem();
   char *endptr;
   
   lines = strtot<uint32_t>(data, &endptr);
   os << lines << std::endl;
   auto dummy = mf.get_line();
   // pts.reserve(lines);
   
   while (linecount < lines)
   {

      auto t = parser::read_gmsh_point_line<double>(endptr, &endptr);
      
      os << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t) << std::endl;
      
      //std::array<double,3> point = { std::get<0>(t),std::get<1>(t),std::get<2>(t) };      
      //pts.push_back(point);
      
      /* Do something with that point */
      dummy = mf.get_line();
      linecount++;
   }
   
   //std::cout << "    Reading points: " << linecount;
   //std::cout << "/" << lines << " - " << tc << " seconds" << std::endl;
   /************************ Read useless stuff ************************/
   
   for ( size_t i=0; i<2; i++)
   {
      auto t = mf.get_line();
      // std::cout << t << std::endl;
   }
   /************************ Read elements ************************/
   linecount = 0;
   const char *new_data = mf.mem();
   lines = strtot<uint32_t>(new_data, &endptr);
   tree_element_id = lines;
   // os << lines << std::endl;

   std::vector< tm_tuple > temp_hex;
   std::vector< surface_type > temp_surf_tri;
   std::vector< uint32_t > temp_surf_lab;
   // std::vector< em_tuple > temp_bnd_edge; 
   temp_hex.reserve(lines);
   temp_surf_lab.reserve(lines);
   temp_surf_tri.reserve(lines);
   uint32_t num_of_hex=0;

   while (linecount < lines)
   {
      dummy = mf.get_line();
      const char *new_data = mf.mem();
      auto t = parser::read_element_line<uint32_t>(new_data, &endptr);
      auto element_label = std::get<0>(t);
      
      if (element_label == 5)
      {
         num_of_hex++;
         
         uint32_t       p0(std::get<2>(t));
         uint32_t       p1(std::get<3>(t));
         uint32_t       p2(std::get<4>(t));
         uint32_t       p3(std::get<5>(t));
         uint32_t       p4(std::get<6>(t));
         uint32_t       p5(std::get<7>(t));
         uint32_t       p6(std::get<8>(t));
         uint32_t       p7(std::get<9>(t));
         uint32_t       d(std::get<1>(t));
         
         auto vol = volume_type(p0, p1, p2, p3,
                                p4, p5, p6, p7);
         
         auto tuple = std::make_tuple(vol, d);
         temp_hex.push_back( tuple );
      }
      else if (element_label == 3)
      {
         // std::cout << "    I'm here with a quad!" << std::endl;
         uint32_t       p0(  std::get<2>(t) );
         uint32_t       p1(  std::get<3>(t) );
         uint32_t       p2(  std::get<4>(t) );
         uint32_t       p3(  std::get<5>(t) );
         uint32_t       bid( std::get<1>(t) );
         
         surface_type   tri{ p0, p1, p2, p3 };
         temp_surf_tri.push_back(tri);
         temp_surf_lab.push_back(bid);
      }
      
      linecount++;
   }
   
   os << num_of_hex << std::endl;
   
   
   std::vector<label_surface_type> temp_tri0;
   temp_tri0.resize(4*lines);

   for (auto tet : temp_hex)
   {
      auto t = std::get<0>(tet);
      //volumes.push_back(t);
      
      uint32_t       p0(std::get<0>(t));
      uint32_t       p1(std::get<1>(t));
      uint32_t       p2(std::get<2>(t));
      uint32_t       p3(std::get<3>(t));
      uint32_t       p4(std::get<4>(t));
      uint32_t       p5(std::get<5>(t));
      uint32_t       p6(std::get<6>(t));
      uint32_t       p7(std::get<7>(t));
      
      os << std::get<1>(tet) << " " << p0+1 << " " << p1+1 << " " << p2+1 << " " << p3+1 
                             << " " << p4+1 << " " << p5+1 << " " << p6+1 << " " << p7+1
         << std::endl;
   }
   
   os << 0 << std::endl;
   os.close();
   
   bool ret_val = read_mesh(sng.str(),intersurface,physical_edges,physical_nodes);
   std::remove(sng.str().c_str());

   return ret_val;
}

bool lean_cohomology :: read_mesh(const std::string& _filename, std::vector<uint32_t>& intersurface, 
                                  std::vector<uint32_t>& physical_edges, std::vector<uint32_t>& physical_nodes)
{   
   timecounter tc, tctot, tprof;
   bool conductor_on_bnd = false;
   
   /* Open file */
   if (_filename.size() == 0)
   {
      std::cout << "    Invalid mesh file name" << std::endl;
      return false;
   }
   
   uint32_t   lines, linecount;
   
   mapped_file mf(_filename);
   
   std::cout << "     * * * Reading formatted mesh file * * * ";
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
         std::cout << "    Reading points: " << linecount;
         std::cout << "/" << lines << "\r";
         std::cout.flush();
      }

      auto t = parser::read_point_line<double>(endptr, &endptr);
      
      
      std::array<double,3> point = { std::get<0>(t),std::get<1>(t),std::get<2>(t) };      
      pts.push_back(point);
      
      /* Do something with that point */
      
      linecount++;
   }
   tc.toc();
   
   std::cout << "    Reading points: " << linecount;
   std::cout << "/" << lines << " - " << tc << " seconds" << std::endl;
   
   /************************ Read tetrahedra ************************/
   linecount = 0;
   
   lines = strtot<uint32_t>(endptr, &endptr);
   std::vector< tm_tuple > temp_hex;
   temp_hex.reserve(lines);
   std::vector<label_surface_type> temp_quad0;
   
   tc.tic();
   while (linecount < lines)
   {
      if ( (linecount%100000) == 0 )
      {
         std::cout << "    Reading hexahedra: " << linecount;
         std::cout << "/" << lines << "\r";
         std::cout.flush();
      }
      
      auto t = parser::read_hex_line<uint32_t>(endptr, &endptr);
      
      temp_quad0.push_back(std::make_pair(surface_type{0,0,0,0},0));
      temp_quad0.push_back(std::make_pair(surface_type{0,0,0,0},0));
      temp_quad0.push_back(std::make_pair(surface_type{0,0,0,0},0));
      temp_quad0.push_back(std::make_pair(surface_type{0,0,0,0},0));
      temp_quad0.push_back(std::make_pair(surface_type{0,0,0,0},0));
      temp_quad0.push_back(std::make_pair(surface_type{0,0,0,0},0));
      
      std::vector<uint32_t> vechex(8);
      vechex[0] = std::get<1>(t);
      vechex[1] = std::get<2>(t);
      vechex[2] = std::get<3>(t);
      vechex[3] = std::get<4>(t);
      vechex[4] = std::get<5>(t);
      vechex[5] = std::get<6>(t);
      vechex[6] = std::get<7>(t);
      vechex[7] = std::get<8>(t);

      //std::sort(vectet.begin(),vectet.end());

      uint32_t       d(std::get<0>(t));
      
      auto tuple = std::make_tuple(volume_type(vechex[0], vechex[1], vechex[2], vechex[3],
                                               vechex[4], vechex[5], vechex[6], vechex[7]), d);
      temp_hex.push_back( tuple );
      
      linecount++;
   }
   tc.toc();
   
   std::cout << "    Reading hexahedra: " << linecount;
   std::cout << "/" << lines  << " - " << tc << " seconds" << std::endl;
   
   /************************ Sort ************************/
   std::cout << "    Sorting data...";
   std::cout.flush();
   
   tc.tic();
   
   /* sort tetrahedra, make unique and move them in geometry */
   // tprof.tic();
   // struct {
   //    bool operator()(const tm_tuple& t1, const tm_tuple& t2)
   //    {
   //       return (std::get<0>(t1) < std::get<0>(t2));
   //    }
   // } mycomp;

   // std::sort(temp_hex.begin(), temp_hex.end(), mycomp);
   
   // _vtf_list.resize(lines);
   uint32_t tot=0;

   std::vector<int32_t> vol_signs;
   
   vol_signs.reserve(lines);
   volumes.reserve(lines);
   domains.reserve(lines);
   size_t hexind = 0;
   for (auto tet : temp_hex)
   {
      auto t = std::get<0>(tet);
      volumes.push_back(t);
      
      uint32_t       p0(std::get<0>(t));
      uint32_t       p1(std::get<1>(t));
      uint32_t       p2(std::get<2>(t));
      uint32_t       p3(std::get<3>(t));
      uint32_t       p4(std::get<4>(t));
      uint32_t       p5(std::get<5>(t));
      uint32_t       p6(std::get<6>(t));
      uint32_t       p7(std::get<7>(t));

      temp_quad0[tot]           = std::make_pair(surface_type{p0, p1, p2, p3},tot);
      temp_quad0[tot+lines]     = std::make_pair(surface_type{p4, p5, p6, p7},tot+lines);
      temp_quad0[tot+2*lines]   = std::make_pair(surface_type{p0, p3, p7, p4},tot+2*lines);
      temp_quad0[tot+3*lines]   = std::make_pair(surface_type{p1, p2, p6, p5},tot+3*lines);
      temp_quad0[tot+4*lines]   = std::make_pair(surface_type{p0, p4, p5, p1},tot+4*lines);
      temp_quad0[tot+5*lines]   = std::make_pair(surface_type{p2, p3, p7, p6},tot+5*lines);
      tot++;

      std::array<double,3> v1 ({ pts[p1][0]-pts[p0][0], pts[p1][1]-pts[p0][1],pts[p1][2]-pts[p0][2] });
      std::array<double,3> v2 ({ pts[p3][0]-pts[p0][0], pts[p3][1]-pts[p0][1],pts[p3][2]-pts[p0][2] });
      std::array<double,3> v3 ({ pts[p4][0]-pts[p0][0], pts[p4][1]-pts[p0][1],pts[p4][2]-pts[p0][2] });
      
      int32_t sgn  = this->stddot(v1, stdcross(v2, v3))/double(6)>0? 1 : -1;
      vol_signs.push_back(sgn);
      domains.push_back(std::get<1>(tet));
   }

   std::vector<tm_tuple>().swap(temp_hex);
   std::vector<uint32_t> labels(6*lines);
   surfaces.reserve(6*lines);
   
   unique(temp_quad0, labels); //this also fills the surfaces vector
   std::vector<label_surface_type>().swap(temp_quad0);
   
   _ftv_list.resize(surfaces.size());
   _vtf_list.reserve(volumes.size());
   
   for (uint32_t k=0; k<lines; ++k)
   {
      // std::cout << labels[k] << " " << _ftv_list.size() << std::endl;
      sgnint32_t<int32_t> v1(k,-vol_signs[k]);
      sgnint32_t<int32_t> v2(k,vol_signs[k]);
      
      _ftv_list[labels[k]].push_back(v1); 
      _ftv_list[labels[k+  lines]].push_back(v2); 
      _ftv_list[labels[k+2*lines]].push_back(v1);
      _ftv_list[labels[k+3*lines]].push_back(v2); 
      _ftv_list[labels[k+4*lines]].push_back(v1);
      _ftv_list[labels[k+5*lines]].push_back(v2);     
      
      sgnint32_t<int32_t> f1(labels[k],-vol_signs[k]);
      sgnint32_t<int32_t> f2(labels[k+  lines], vol_signs[k]);
      sgnint32_t<int32_t> f3(labels[k+2*lines],-vol_signs[k]);
      sgnint32_t<int32_t> f4(labels[k+3*lines], vol_signs[k]);
      sgnint32_t<int32_t> f5(labels[k+4*lines],-vol_signs[k]);
      sgnint32_t<int32_t> f6(labels[k+5*lines], vol_signs[k]);
      
      std::vector<sgnint32_t<int32_t> > dummy(6);
      _vtf_list.push_back(dummy);
      _vtf_list[k][0] = f1;
      _vtf_list[k][1] = f2; 
      _vtf_list[k][2] = f3;
      _vtf_list[k][3] = f4;
      _vtf_list[k][4] = f5;
      _vtf_list[k][5] = f6;
   }
   
   std::vector<uint32_t>().swap(labels);
   lines = surfaces.size();
   intersurface.resize(lines,0);
   // _fte_list.resize(lines);
   tot = 0;
   std::vector<label_edge_type> temp_edge0(4*lines);
   
   for (auto t : surfaces)
   {
      uint32_t       p0(std::get<0>(t));
      uint32_t       p1(std::get<1>(t));
      uint32_t       p2(std::get<2>(t));
      uint32_t       p3(std::get<3>(t));

      temp_edge0[tot]          =  std::make_pair(edge_type{p0, p1},tot);
      temp_edge0[tot+lines]    =  std::make_pair(edge_type{p1, p2},tot+lines);
      temp_edge0[tot+2*lines]  =  std::make_pair(edge_type{p2, p3},tot+2*lines);
      temp_edge0[tot+3*lines]  =  std::make_pair(edge_type{p0, p3},tot+3*lines);

      tot++;
   }
   
   std::vector<uint32_t> e_labels(4*lines);

   unique(temp_edge0, e_labels); // also fills the edges vector

   std::vector<label_edge_type>().swap(temp_edge0);
   _etf_list.resize(edges.size());
   physical_edges.resize(edges_size(),0);
   edge_in_conductor.resize(edges_size());
   
   std::ofstream cuts_pre_minimization, cuts_post_minimization;
   if (debuggy)
   {
      cuts_pre_minimization.open("./cuts_pre_minimization.txt");
      cuts_post_minimization.open("./cuts_post_minimization.txt");
   }
   
   for (uint32_t k=0; k<lines; ++k)
   {
      sgnint32_t<int32_t> f1(k, 1);
      sgnint32_t<int32_t> f2(k, 1);
      sgnint32_t<int32_t> f3(k, 1);
      sgnint32_t<int32_t> f4(k,-1);
      
      _etf_list[e_labels[k]].push_back(f1);
      _etf_list[e_labels[k+lines]].push_back(f2); 
      _etf_list[e_labels[k+2*lines]].push_back(f3);
      _etf_list[e_labels[k+3*lines]].push_back(f4);
      
      sgnint32_t<int32_t> e1(e_labels[k],         1);
      sgnint32_t<int32_t> e2(e_labels[k+lines],   1);
      sgnint32_t<int32_t> e3(e_labels[k+2*lines], 1);
      sgnint32_t<int32_t> e4(e_labels[k+3*lines],-1);
      
      std::vector<sgnint32_t<int32_t> > dummy(4);      
      _fte_list.push_back(dummy);
      
      _fte_list[k][0] = e1;
      _fte_list[k][1] = e2;
      _fte_list[k][2] = e3;
      _fte_list[k][3] = e4;
      
      auto vols = _ftv_list[k];
      
      switch (vols.size()) 
      {
         case 2: 
         {
            auto vol1= abs(*vols.begin());
            auto vol2= abs(*(std::prev(vols.end())));
            
            if (is_conductor(vol1) !=  is_conductor(vol2))
            {
               intersurface[k]++;
               
               if (debuggy)
               {
                  //~ std::cout << __FILE__ << ":" << __LINE__ << std::endl;
                  if (is_conductor(vol1))
                  {
                     cuts_pre_minimization << print_face(0,k,(*vols.begin()).Sgn()>0,122, 122, 122);
                     cuts_post_minimization << print_face(0,k,(*vols.begin()).Sgn()>0,122, 122, 122);
                  }
                  else// if (domains[vol2] == conductor_id )
                  {
                     cuts_pre_minimization  << print_face(0,k,(*vols.rbegin()).Sgn()>0,122, 122, 122);
                     cuts_post_minimization << print_face(0,k,(*vols.rbegin()).Sgn()>0,122, 122, 122);
                  }
               }
               
               for (auto ee : _fte_list[k])
               {
                  if (!physical_edges[abs(ee)])
                     e2d++;
                  physical_edges[abs(ee)]++;
                  edge_in_conductor[abs(ee)] = false;
               }
               f2d++;
            }
            else if (is_conductor(vol1) && is_conductor(vol2))
            {
               for (auto ee : _fte_list[k])
                  if (!physical_edges[abs(ee)])
                     edge_in_conductor[abs(ee)] = true;
            }
            else if (!is_conductor(vol1) && !is_conductor(vol2))
            {
               for (auto ee : _fte_list[k])
                  edge_in_conductor[abs(ee)] = false;
            }
            break;
         }
         case 1:
         {
            auto vol1= abs(*vols.begin());
            if (is_conductor(vol1))
            {
               conductor_on_bnd = true;
               intersurface[k]++;
               
               if (debuggy)
               {
                  cuts_pre_minimization  << print_face(0,k,(*vols.begin()).Sgn()>0,122, 122, 122);
                  cuts_post_minimization << print_face(0,k,(*vols.begin()).Sgn()>0,122, 122, 122);
               }
               
               for (auto ee : _fte_list[k])
               {
                  if (!physical_edges[abs(ee)])
                     e2d++;
                  physical_edges[abs(ee)]++;
                  edge_in_conductor[abs(ee)] = false;
               }
               f2d++;

               for (auto ee : _fte_list[k])
                  if (!physical_edges[abs(ee)])
                     edge_in_conductor[abs(ee)] = true;
            }
            else
               for (auto ee : _fte_list[k])
                  edge_in_conductor[abs(ee)] = false;
            break;
         }
         case 0:
         {
            throw std::invalid_argument("    ERROR: The whole domain must be a differentiable manifold!");
            break;
         }
      }   
   }
   
   if (debuggy)
   {
      cuts_pre_minimization.close();
      cuts_post_minimization.close();
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
      
      std::vector<sgnint32_t<int32_t> > dummy(2);      
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
   
   std::cout << "    done - " << tc << " seconds" << std::endl;
   
   if (conductor_on_bnd)
      std::cout << "    WARNING: conductors intersect the boundary of the mesh."
                << std::endl << "    (Ignore if this is wanted)" << std::endl;
      
   /************************ Read boundary surfaces ************************/
   linecount = 0;
   // auto num_of_hex=lines;
   lines = strtot<uint32_t>(endptr, &endptr);
   
   tc.tic();
   while (linecount < lines)
   {
      /*if (ifs.fail())
      {
         std::cout << "    Error while reading boundary surfaces" << std::endl;
         return false;
      }*/
      
      if ( (linecount%50000) == 0 )
      {
         std::cout << "    Reading facets: " << linecount;
         std::cout << "/" << lines << "\r";
         std::cout.flush();
      }
      
      auto t = parser::read_quad_line<uint32_t>(endptr, &endptr);

      uint32_t       p0(  std::get<1>(t) );
      uint32_t       p1(  std::get<2>(t) );
      uint32_t       p2(  std::get<3>(t) );
      uint32_t       p3(  std::get<4>(t) );
      uint32_t       bid( std::get<0>(t) );
      
      surface_type   tri{ p0, p1, p2, p3 };
      
      if (!physical_surfaces[bid].size())
         physical_surfaces[bid].resize(surfaces.size(),0);

      auto itor = std::lower_bound(surfaces.begin(),surfaces.end(),tri);
      physical_surfaces[bid][std::distance(surfaces.begin(),itor)]++;
      
      linecount++;
   }
   tc.toc();
   
   std::cout << "    Reading facets: " << linecount;
   std::cout << "/" << lines  << " - " << tc << " seconds"  << std::endl;

   tctot.toc();

   // std::cout << cyan << "Total time spent in reading mesh: ";
   // std::cout << tctot << " seconds" << nocolor << std::endl;
   
   // for (const auto& nn : _vtf_list)
      // assert(nn.size() == 6);
   
   // for (const auto& nn : _ftv_list)
      // assert(nn.size() == 2 || nn.size() == 1);
   
   // for (const auto& nn : _fte_list)
      // assert(nn.size() == 4);
   
   // for (const auto& nn : _etf_list)
   // {
      // assert(nn.size() >= 2 && nn.size()<=4);
   // }
   
   // for (const auto& nn : _nte_list)
      // assert(nn.size() >= 3 && nn.size() <= 6);
   // for (const auto& nn : _etn_list)
      // assert(nn.size() == 2);   
   
   return true;
}

void lean_cohomology :: MinCost(const std::vector<int32_t>& start_s, uint32_t indigen)
{
   uint32_t i, j, k;
   int e_or;
   
   k=0;
   
   // Initial feasible solutions for primal and dual
   std::vector<int32_t> fp(edges_size()), dummy_cap(edges_size(),1);
   capacities = std::move(dummy_cap);
   flow=fp;//fp=0;
   face_coeff_in_the_chain = start_s;
   // int maxflow=0;
   
   // useful_variables
   uint32_t source, sink;
   int32_t step;
   bool RMF_SUCCESS;
   pair<sgnint32_t<int32_t> ,sgnint32_t<int32_t> > adj_v;
   std::vector<int32_t> new_min_cut(edges_size()), dummy_cut(edges_size());
   
   timecounter tc;
   tc.tic();
   
   while (true)
   {
      while (k < edges_size() && ( (face_coeff_in_the_chain.at(k)==0) ||
                                             ((face_coeff_in_the_chain.at(k)>0 && flow.at(k) == capacities.at(k)) || 
                     (face_coeff_in_the_chain.at(k)<0 && flow.at(k)==-capacities.at(k))) ) )
         k++;
      
      if (k==edges_size())
      {
         // std::cout << "    Completed solution!" << std::endl;
         break;
      }
      else
      {         
         adj_v=make_pair(*(etn(k).begin()),*std::prev((etn(k).end())));
         
         if (adj_v.first < 0)
            std::swap(adj_v.first,adj_v.second);
         
         if (face_coeff_in_the_chain.at(k)<0)
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

      // Resetting should not really be needed
      k=0;
   }
      
   tc.toc();
   //~ std::cout << "    Minimizing support of generator n. " << indigen+1 << " took " << tc << " seconds" << std::endl;
   std::cout << "    Minimizing support of generators took " << tc << " seconds" << std::endl;

   std::ofstream cuts_post_minimization, h1_post_minimization;
   
   if (debuggy)
      cuts_post_minimization.open("./cuts_post_minimization.txt", std::ofstream::out | std::ofstream::app);   
   
   h1_post_minimization.open("./h1.txt", std::ofstream::out | std::ofstream::app);
   
   std::vector<int32_t> coeff_cf_array;
   std::vector<std::array<uint32_t,2> > coeff_id_array;
   uint32_t n_coeff=0;

   for ( uint32_t it=0; it<edges_size(); ++it)
   {
      if (face_coeff_in_the_chain[it] != 0)
      {
         if (!edge_in_conductor[it])
         {
            if (debuggy)
            { 
               auto dual_face_vector = print_dual_face(indigen%9+1,it,face_coeff_in_the_chain[it]>0,255,255,0);
               for (auto df : dual_face_vector)
                  cuts_post_minimization << df;
            }
            n_coeff++;
            coeff_id_array.push_back(std::array<uint32_t,2>({abs(*(_etn_list[it].begin())),
                               abs(*(_etn_list[it].rbegin()))}));
            coeff_cf_array.push_back(face_coeff_in_the_chain[it]);
         }
      }
   }

   h1_post_minimization << n_coeff << std::endl;
   for (uint32_t ii=0; ii<n_coeff; ++ii)
      h1_post_minimization << coeff_id_array[ii][0]+1 << " "
                           << coeff_id_array[ii][1]+1 << " "
                           << coeff_cf_array[ii]    << std::endl;
   
   if (debuggy)
      cuts_post_minimization.close();
   
   h1_post_minimization.close();
   
   return;
}

bool lean_cohomology :: ResMaxFlow(std::vector<int32_t>& new_min_cut, const uint32_t& k, const uint32_t& s, const uint32_t& t, const int32_t& nmc_step)
{
   // Minimum cost network flow with restriction on already achieved orthogonalities

   std::vector<int> new_fp(edges_size());
   new_min_cut.at(k)+=nmc_step;
   // std::map<uint32_t,uint32_t> parent;
   std::vector<uint32_t> parent(nodes_size(),nodes_size());
   std::vector<uint32_t> colour(nodes_size());
   std::vector<uint32_t> par_edge(nodes_size());

   // parent.at(s)=0;
   colour.at(s)=1;
   par_edge.at(s)=k;

   std::vector<uint32_t> queue;
   uint32_t i=0;
   queue.push_back(s);
   
   uint32_t qtop, curr_f, current_node;
   sgnint32_t<int32_t>  adj_v;
   
   while (i<queue.size() && colour.at(t)==0)
   {
      qtop=queue.at(i);
      // adj_f=edge_list(qtop,:);
      for (auto adj_f : nte(qtop))
      {   
         
         curr_f=abs(adj_f);
         current_node=abs(*(etn(curr_f).begin()));
         
         if (current_node==qtop)
         {
            current_node=abs(*std::prev(etn(curr_f).end()));
            adj_v=*std::prev(etn(curr_f).end());
         }
         else
            adj_v=*(etn(curr_f).begin());
         
         // std::cout << "    faccia corrente: " << curr_f << "dim_flow " << flow.size() << " " << capacities.size() << std::endl;
         
         if ((current_node != t) || (qtop != s))
         {
            // std::cout << curr_f << std::endl;
            if (adj_v<0)
            {
               if (face_coeff_in_the_chain.at(curr_f) >= 0 && flow.at(curr_f) < capacities.at(curr_f))
               {
                  if (colour.at(current_node)==0)
                  {
                     // dimQ=dimQ+1;
                     queue.push_back(current_node);
                     colour.at(current_node)+=1;
                     parent[current_node]=qtop;
                     par_edge.at(current_node)=curr_f;
                     new_fp.at(curr_f)=1;
                  }   
               }
               
               new_min_cut.at(curr_f)+=1;
            }
            else 
            {
               if (face_coeff_in_the_chain.at(curr_f) <= 0 && flow.at(curr_f) > -capacities.at(curr_f))
               {
                  if (colour.at(current_node)==0)
                  {
                     queue.push_back(current_node);
                     colour.at(current_node)+=1;
                     parent[current_node]=qtop;
                     par_edge.at(current_node)=curr_f;
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
      
      while (parent.at(adj_v) < nodes_size())
      {
         curr_f=par_edge.at(adj_v);
         flow.at(curr_f)+=new_fp.at(curr_f);
         adj_v=parent.at(adj_v);
      }
      
      curr_f=k;

      flow.at(curr_f)=flow.at(curr_f)-nmc_step;
      
      return 1;
   }
   else
      return 0;
}

void lean_cohomology :: UpdateMinCut(std::vector<int32_t>& new_min_cut)
{
   assert( face_coeff_in_the_chain.size() == new_min_cut.size() );
   
   for (uint32_t i=0; i<face_coeff_in_the_chain.size(); i++)
      face_coeff_in_the_chain.at(i)+=new_min_cut.at(i);
   
   return;
}

double lean_cohomology :: LinkingNumber(const std::vector<std::array<double,3> >& gamma0, const std::vector<std::array<double,3> >& gamma1)
{
   
   double n = 0;
   
   std::array<double,3> a,b,c,d;

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
}

std::vector<double> lean_cohomology :: LinkingNumber(const std::vector<std::array<double,3> >& gamma0)
{
   std::vector<double> ret(n_lazy,0);
   
   //~ std::array<double,3> a,b,c,d;
   for (uint32_t f=0;f<surfaces_size();++f)
   {
      if (_ftv_list[f].size()==2)
      {
         for (auto gen : HomoCoHomo.second[f])
         {
            double coeff = signbit(gen) ? -1 : 1;
            for (uint32_t i = 0; i <gamma0.size(); ++i)
            {
               std::array<double,3> a,b,c,d;
               std::array<double,3> gamma11,gamma10;
               gamma11[0] =  vol_barycenter(abs(_ftv_list[f][0]))[0];
               gamma10[0] =  vol_barycenter(abs(_ftv_list[f][1]))[0];
               gamma11[1] =  vol_barycenter(abs(_ftv_list[f][0]))[1];
               gamma10[1] =  vol_barycenter(abs(_ftv_list[f][1]))[1];
               gamma11[2] =  vol_barycenter(abs(_ftv_list[f][0]))[2];
               gamma10[2] =  vol_barycenter(abs(_ftv_list[f][1]))[2];
               
               if ( coeff < 0 )
                  std::swap(gamma11,gamma10);
                  
               if ( !(_ftv_list[f][0]<0) )
                  std::swap(gamma11,gamma10);
               
               a[0] = gamma10[0]-gamma0.at(i)[0];
               a[1] = gamma10[1]-gamma0.at(i)[1];
               a[2] = gamma10[2]-gamma0.at(i)[2];
               b[0] = gamma10[0]-gamma0.at((i+1) % gamma0.size())[0];
               b[1] = gamma10[1]-gamma0.at((i+1) % gamma0.size())[1];
               b[2] = gamma10[2]-gamma0.at((i+1) % gamma0.size())[2];
               c[0] = gamma11[0]-gamma0.at((i+1) % gamma0.size())[0];
               c[1] = gamma11[1]-gamma0.at((i+1) % gamma0.size())[1];
               c[2] = gamma11[2]-gamma0.at((i+1) % gamma0.size())[2];
               d[0] = gamma11[0]-gamma0.at(i)[0];
               d[1] = gamma11[1]-gamma0.at(i)[1];
               d[2] = gamma11[2]-gamma0.at(i)[2];
               
               ret[abs(gen)-1] += (0.25/PI)*SolidAngleQuadrilateral(a, b, c, d);
            }
         }
      }
   }

    return ret;
}


double lean_cohomology :: SolidAngleQuadrilateral(const std::array<double,3>& a, 
                                      const std::array<double,3>& b, 
                                      const std::array<double,3>& c,
                                      const std::array<double,3>& d)
{
   return SolidAngleTriangle(a, b, c) + SolidAngleTriangle(c, d, a);
}

double lean_cohomology :: SolidAngleTriangle(const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c)
{
   double determ= a[0]*(b[1]*c[2] - 
                  b[2]*c[1])+a[1]*(b[2]*c[0] - 
                  b[0]*c[2])+a[2]*(b[0]*c[1] - 
                  b[1]*c[0]);

   double al = sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2));
   double bl = sqrt(pow(b[0],2)+pow(b[1],2)+pow(b[2],2));
   double cl = sqrt(pow(c[0],2)+pow(c[1],2)+pow(c[2],2));

   double div = al * bl * cl +
                cl*(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]) + 
                bl*(c[0]*a[0]+c[1]*a[1]+c[2]*a[2]) + 
                al*(c[0]*b[0]+c[1]*b[1]+c[2]*b[2]);

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
}

const std::vector<sgnint32_t<int32_t> >& lean_cohomology :: vtf(const int32_t& v_id) const
{
   if ( _vtf_list.size() == 0 )
      throw std::invalid_argument("you have to populate the list first");
   else
      return _vtf_list[uint32_t(abs(v_id))];
}

bool lean_cohomology :: is_conductor(const uint32_t& v)
{
   return conductor_bool[domains[v]];
}

const std::vector<sgnint32_t<int32_t> >& lean_cohomology :: fte(const int32_t& f_id) const
{
   if ( _fte_list.size() == 0 )
      throw std::invalid_argument("you have to populate the list first");
   else
      return _fte_list[uint32_t(abs(f_id))];
}

const std::vector<sgnint32_t<int32_t> >& lean_cohomology :: ftv(const int32_t& f_id) const
{
   if ( _ftv_list.size() == 0 )
      throw std::invalid_argument("you have to populate the list first");
   else
      return _ftv_list[uint32_t(abs(f_id))];
}

const std::vector<sgnint32_t<int32_t> >& lean_cohomology :: etf(const int32_t& e_id) const
{
   if ( _etf_list.size() == 0 )
      throw std::invalid_argument("you have to populate the list first");
   else
      return _etf_list[uint32_t(abs(e_id))];
}

const std::vector<sgnint32_t<int32_t> >& lean_cohomology :: etn(const int32_t& e_id) const
{
   if ( _etn_list.size() == 0 )
      throw std::invalid_argument("you have to populate the list first");
   else
      return _etn_list[uint32_t(abs(e_id))];
}

const std::vector<sgnint32_t<int32_t> >& lean_cohomology :: nte(const int32_t& n_id) const
{
   if ( _nte_list.size() == 0 )
      throw std::invalid_argument("you have to populate the list first");
   else
      return _nte_list[uint32_t(abs(n_id))];
}   

std::vector<double> lean_cohomology :: face_barycenter(const uint32_t& f)
{
   std::vector<double> bc(3,0);
   
   auto p0 = std::get<0>(surfaces[f]);
   auto p1 = std::get<1>(surfaces[f]);
   auto p2 = std::get<2>(surfaces[f]);
   
   double one_third = double(1)/double(3);
   
   bc[0] = one_third*(pts[p0][0]+pts[p1][0]+pts[p2][0]);
   bc[1] = one_third*(pts[p0][1]+pts[p1][1]+pts[p2][1]);
   bc[2] = one_third*(pts[p0][2]+pts[p1][2]+pts[p2][2]);
   
   return bc;
}

std::vector<double> lean_cohomology :: vol_barycenter(const uint32_t& v)
{
   std::vector<double> bc(3,0);
   
   auto p0 = std::get<0>(volumes[v]);
   auto p1 = std::get<1>(volumes[v]);
   auto p2 = std::get<2>(volumes[v]);
   auto p3 = std::get<3>(volumes[v]);
   
   
   bc[0] = 0.25*(pts[p0][0]+pts[p1][0]+pts[p2][0]+pts[p3][0]);
   bc[1] = 0.25*(pts[p0][1]+pts[p1][1]+pts[p2][1]+pts[p3][1]);
   bc[2] = 0.25*(pts[p0][2]+pts[p1][2]+pts[p2][2]+pts[p3][2]);
   
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
   std::vector<std::array<double,3> > n;
   
   for (auto ee : fte(f))
   {
      auto n1 = std::get<0>(edges[abs(ee)]);
      auto n2 = std::get<1>(edges[abs(ee)]);
      
      nodes.insert(n1);
      nodes.insert(n2);
      // for (auto nn : etn(abs(ee)))
         // nodes.insert(abs(nn));
      
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

   std::array<double,3> v1 { x, y, z };
   
   for (auto f : _etf_list[e])
   {
      std::vector<double> fb = face_barycenter(abs(f));
      
      for (auto v : _ftv_list[abs(f)])
      {
         orient = new_or*f.Sgn()>0;
         std::vector<double> vb =  vol_barycenter(abs(v));

         std::array<double,3> v2 { eb[0]-vb[0], eb[1]-vb[1],eb[2]-vb[2] };
         std::array<double,3> v3 { fb[0]-vb[0], fb[1]-vb[1], fb[2]-vb[2] };
         
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
      auto n1 = std::get<0>(edges[abs(ee)]);
      auto n2 = std::get<1>(edges[abs(ee)]);
      
      nodes.insert(n1);
      nodes.insert(n2);
      // for (auto nn : etn(abs(ee)))
         // nodes.insert(abs(nn));
      
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
   auto n1 =  pts[nn_b];
   auto n2 =  pts[nn_e];
   
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

std::array<double,3> lean_cohomology :: stdcross(const std::array<double,3>& v1, const std::array<double,3>& v2)
{
    std::array<double,3> ret;
    
    ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
    ret[1] = v1[2]*v2[0] - v1[0]*v2[2];
    ret[2] = v1[0]*v2[1] - v1[1]*v2[0];
    
    return ret;
}

double lean_cohomology :: stddot(const std::array<double,3>& v1, const std::array<double,3>& v2)
{
    double acc=0;
    
    for (uint32_t i = 0; i < 3; i++)
        acc += v1[i] * v2[i];
        
    return acc;
}
