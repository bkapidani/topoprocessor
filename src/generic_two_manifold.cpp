//file generic_two_manifold.cpp
#include "lean_cohomology.hpp"

generic_two_manifold :: generic_two_manifold(lean_cohomology *lc)
{
   this->vol_mesh = lc;
   this->f = vol_mesh->f2d;
   this->e = vol_mesh->e2d;
   this->n = vol_mesh->n2d;
}

std::pair<h1_2d_basis,thinned_currents> generic_two_manifold :: H_to_CoH(const std::vector<uint32_t>& physical_surface, const std::vector<uint32_t>& physical_edges, const std::vector<uint32_t>& physical_nodes)
{
   // h1_2d_basis h1b;
   // thinned_currents tc;
   std::pair<h1_2d_basis,thinned_currents> p;
   std::vector<uint32_t> remaining_edges;
   std::vector<bool> cotree_edges(vol_mesh->edges_size(),true);
   uint32_t k=0;
   std::vector<uint32_t> p_queue, d_queue;
   
   std::vector<bool>              p_colour(this->vol_mesh->nodes_size(),false), d_colour(this->vol_mesh->surfaces_size(),false);
   std::map<uint32_t,uint32_t>    p_parent, p_paredge, p_distance, d_parent, d_paredge, d_distance;
   // std::map<uint32_t,uint32_t>    t_parent, t_paredge;
   // std::vector<bool>                       p_colour(this->vol_mesh->nodes_size(),false), d_colour(this->vol_mesh->surfaces_size(),false);
   // std::vector<uint32_t>                   p_parent(this->vol_mesh->nodes_size()), p_paredge(this->vol_mesh->nodes_size());
   // std::vector<uint32_t>                   p_distance(this->vol_mesh->nodes_size(),std::numeric_limits<int>::max());
   // std::vector<uint32_t>                  t_parent(this->vol_mesh->volumes_size()), t_paredge(this->vol_mesh->volumes_size());

   uint32_t j=0;
   uint32_t root_dual=0;
   timecounter t_h1, t_tc;   

   h1b.resize(vol_mesh->edges_size());
   tc.resize(vol_mesh->surfaces_size());
   // p_queue.reserve(this->n);
   // d_queue.reserve(this->f);
   
   k=0;
   t_h1.tic();
   while (p_queue.size()<this->n)
   {
      // std::cout << p_queue.size() << "------" << this->n << "-----" << vol_mesh->nodes_size() << std::endl;
      uint32_t root=0;
      while (p_colour[root] || !physical_nodes[root])
      {
         // std::cout << root << "----" << std::endl;
         root++;
      }
      // if (root < 100)
         // std::cout << "    ----" << root << "----" << std::endl;
      p_colour[root]=true;
      p_distance[root]=0;
      p_queue.push_back(root);
      uint32_t qtop;
      
      while (k < p_queue.size())
      {   
         uint32_t qtop = p_queue[k];
         // std::cout << "    ----" << qtop << "----" << std::endl;
         for (const auto& signed_ee : vol_mesh->nte(qtop))
         {
            uint32_t ee = abs(signed_ee);
            
            if (physical_edges[ee]>0)
            {
               for (const auto& signed_nn : vol_mesh->etn(ee))
               {
                  uint32_t nn = abs(signed_nn);
                  
                  if (nn != qtop)
                  {

                     if (p_colour[nn] == false )
                     {
                        p_colour[nn]=true;
                        // p_distance[nn]=p_distance[qtop]+p_length[ee];
                        p_distance[nn]=p_distance[qtop]+1;
                        p_parent[nn]=qtop;
                        p_paredge[nn]=ee;
                        cotree_edges[ee]=false;
                        p_queue.push_back(nn);
                     }
                  }
               }
            }
         }
         
         k++;
      }
      // std::cout << "    ----" << root << "----" << std::endl;
   }

   k=0;
   
    while (d_queue.size()<this->f)
   {
      root_dual=0;
      // std::cout << d_queue.size() << "------" << this->f << std::endl;
      
      while (d_colour[root_dual] || !physical_surface[root_dual])
         root_dual++;
      // std::cout << "    ----" << root_dual << "----" << std::endl;
      d_colour[root_dual]=true;
      // d_distance[root_dual]=0;
      d_queue.push_back(root_dual);
      uint32_t qtop;      
      
      
      while (k < d_queue.size())
      {   
         // std::cout << "    test_bug" << std::endl;
         qtop = d_queue[k];
         
         for (auto& signed_ee : vol_mesh->fte(qtop))
         {
            uint32_t ee = abs(signed_ee);
            if (cotree_edges[ee])
            {            
               for (const auto& signed_ff : vol_mesh->etf(ee))
               {   
                  uint32_t ff = abs(signed_ff);
                  if (ff != qtop && physical_surface[ff]>0)
                  {
                     if (d_colour[ff] == false)
                     {
                        // d_distance[ff]=d_distance[qtop]+1;
                        d_colour[ff]=true;
                        // d_parent[ff]=qtop;

                        // d_paredge[ff]=ee
                        
                        d_distance[ff]=d_distance[qtop]+1;
                        d_parent[ff]=qtop;
                        d_paredge[ff]=ee;

                        d_queue.push_back(ff);
                        cotree_edges[ee]=false;
                     }
                     else
                     {
                        cotree_edges[ee]=false;
                        remaining_edges.push_back(ee);
                     }
                  }
                  
               }
            }

         }
         
         k++;
      }
   }
   
   // for (size_t i=0; i<cotree_edges.size(); i++)
      // if (cotree_edges[i] && physical_edges[i]>0)
         // remaining_edges.push_back(i);


   t_h1.toc();
   std::cout << "    Tree primal, tree dual took: " << t_h1 << " s" << std::endl;
   // Homology generators
   int16_t n_gen = 0;

   std::vector<std::thread> h1_multithreading;

   if (remaining_edges.size()>0)
   {   
      // std::ofstream os;
      // os.open("./output/h1_lazy.txt");
         
      std::cout << "    Number of lazy generators : "
                << remaining_edges.size() << std::endl;

      // this->h1b=std::move(h1b);
      // this->tc=std::move(tc);
      this->p_colour=std::move(p_colour);
      this->p_parent=std::move(p_parent);
      this->p_paredge=std::move(p_paredge);
      this->p_distance=std::move(p_distance);
      this->d_colour=std::move(d_colour);
      this->d_parent=std::move(d_parent);
      this->d_paredge=std::move(d_paredge);
      this->d_distance=std::move(d_distance);
      this->physical_edges=std::move(physical_edges);
      this->physical_surface=std::move(physical_surface);      
      
      for (int32_t i=1; i <= remaining_edges.size(); i++)
      {
         // t_h1.tic();
         auto ee = remaining_edges[i-1];
         
         RetrieveGenAndTC(i,ee);
         // h1_multithreading.push_back(std::thread(&generic_two_manifold::RetrieveGenAndTC,this,i,ee));
      }
      // os.close();
      
   }
   
   // for (size_t i=0; i<h1_multithreading.size(); i++)
      // h1_multithreading[i].join();
   //~ if (debuggy)
   //~ {
      //~ std::ofstream os("homology_gens.txt", std::ofstream::out | std::ofstream::app);
      //~ for (uint32_t fff = 0; fff < vol_mesh->surfaces_size(); ++fff)
      //~ {
         //~ if (physical_surface[fff]>0)
            //~ os << vol_mesh->print_face(9,fff,true,160,160,160);
      //~ }
      //~ os.close();
   //~ }
   
   this->Ngen=remaining_edges.size();
   this->remaining_edges=std::move(remaining_edges);
   return std::make_pair(h1b,tc);
}

void generic_two_manifold :: RetrieveGenAndTC(const int16_t& n_gen, const uint32_t& ee)
{
   // n_gen++;
   std::ofstream os;
   auto ff = vol_mesh->etn(ee);
   uint32_t u = abs(*(ff.begin()));
   uint32_t v = abs(*(std::prev(ff.end())));
   uint32_t added_edge_u = ee;
   uint32_t added_edge_v = ee;
   uint32_t added_edge   = ee;
   uint32_t added_edge_small;
   int16_t or_edge=1;
   int16_t or_edge_v = 1;
   int16_t or_edge_u = 1;
   int16_t or_edge_small;
   
   p1_pushback_mutex.lock();
   h1b[added_edge].push_back(n_gen);
   p1_pushback_mutex.unlock();
   if (debuggy)
   {
      os.open("homology_gens.txt", std::ofstream::out | std::ofstream::app);
      os << vol_mesh->print_edge(n_gen,added_edge,true,0,255,255);
   }
   std::vector<uint32_t> HomologyEdges;
   HomologyEdges.push_back(added_edge);
   
   bool from_small_to_big=true;
   int16_t c_pr_edge, c_added_edge;
   uint32_t small,big;
   
   // if (p_distance[u]>p_distance[v])
   // {
      // std::swap(u,v);
      // from_small_to_big = !from_small_to_big;
   // }
   
   while (true)
   {
      if (p_distance[u]<p_distance[v])
      {
         small = u;
         big = v;
         added_edge = added_edge_v;
         or_edge = or_edge_v;
         added_edge_small = added_edge_u;
         or_edge_small = or_edge_u;
         from_small_to_big = !from_small_to_big;
      }
      else
      {
         small = v;
         big = u;
         added_edge = added_edge_u;
         or_edge = or_edge_u;
         added_edge_small = added_edge_v;
         or_edge_small = or_edge_v;
      }
      
      // for (auto ee : vol_mesh->nte(big))
      // {
         // if(physical_edges[abs(ee)])
         // {
            // if (abs(ee) == added_edge)
            // {
               // c_pr_edge    = ee<0 ? -1 : 1;
            // }
            // else if (abs(ee) == p_paredge[big])
            // {
               // c_added_edge = ee<0 ? -1 : 1;
            // }
         // }
      // }
      
      // std::thread pre_tr([&]{
         for (auto nn : vol_mesh->etn(added_edge))
         {
            if (abs(nn) == big)
            {
               c_pr_edge    = nn<0 ? -1 : 1;
               break;
            }
         }
      // });
      
      // std::thread succ_tr([&]{
         for (auto nn : vol_mesh->etn(p_paredge[big]))
         {
            if (abs(nn) == big)
            {
               c_added_edge = nn<0 ? -1 : 1;
               break;
            }
         }
      // });
      
      // pre_tr.join();
      // succ_tr.join();
      
      if (c_added_edge<0)
         or_edge*= c_pr_edge;
      else
         or_edge*=-c_pr_edge;

      added_edge = p_paredge[big];
      int16_t coeff = (p_parent[big]<big)? 1:-1;
      uint32_t s_edge,t_edge;
      
      if (from_small_to_big)
      {
         p1_pushback_mutex.lock();
         h1b[added_edge].push_back(coeff*n_gen);
         p1_pushback_mutex.unlock();
         
         s_edge = added_edge;
         t_edge = *HomologyEdges.begin();
         HomologyEdges.insert(HomologyEdges.begin(),added_edge);
         if (debuggy)
            os << vol_mesh->print_edge(n_gen,added_edge,coeff>0,0,255,255);
      }
      else
      {
         p1_pushback_mutex.lock();
         h1b[added_edge].push_back((-coeff)*n_gen);
         p1_pushback_mutex.unlock();
         
         s_edge = HomologyEdges.back();
         t_edge = added_edge;            
         HomologyEdges.push_back(added_edge);
         if (debuggy)
            os << vol_mesh->print_edge(n_gen,added_edge,coeff<0,0,255,255);
      }
      
      TwoSidedAlg(n_gen,s_edge,t_edge);
      
      if (p_parent[big] == small)
      {
         // std::vector<uint32_t> last_input;
         s_edge = HomologyEdges.back();
         t_edge = *HomologyEdges.begin();
         TwoSidedAlg(n_gen,s_edge,t_edge);
         
         break;
      }
      
      u=small;
      from_small_to_big = !from_small_to_big;
      added_edge_u=added_edge_small;
      or_edge_u=or_edge_small;
      v=p_parent[big];
      added_edge_v=added_edge;
      or_edge_v=or_edge;
   }
   if (debuggy)
      os.close();
}

/*
std::vector<int> generic_two_manifold :: RetrieveCoH(const std::vector<int>& gen_comb)
{
   std::vector<int> ret(vol_mesh->edges_size(),0);
   for (auto n_gen= 0; n_gen < gen_comb.size(); ++n_gen)
   {
      if (gen_comb[n_gen] != 0)
      {
         std::ofstream os("cohomology_gens.txt", std::ofstream::out | std::ofstream::app);
         int32_t gcvalue = gen_comb[n_gen]; //this is the value i have to multiply for
         
         auto ee = remaining_edges[n_gen];
         auto ff = vol_mesh->etf(ee);
         
         uint8_t nfound = 0;
         uint32_t utest, u, v;
         
         for (auto testf : ff)
         {
            if (nfound == 0)
            {
               utest = abs(testf);
               if (physical_surface[utest]>0)
               {
                  u = utest;
                  ++nfound;
               }
            }
            else if (nfound == 1)
            {
               utest = abs(testf);
               if (physical_surface[utest]>0)
               {
                  v = utest;
                  break;
               }
            }
         }
         
         uint32_t added_edge_u = ee;
         uint32_t added_edge_v = ee;
         uint32_t added_edge   = ee;
         uint32_t added_edge_small;
         int16_t or_edge=1;
         int16_t or_edge_v = 1;
         int16_t or_edge_u = 1;
         int16_t or_edge_small;
         
         p1_pushback_mutex.lock();
         h1b[added_edge].push_back(n_gen);
         p1_pushback_mutex.unlock();
         os << vol_mesh->print_edge(n_gen,added_edge,true,0,255,255);
         std::vector<uint32_t> HomologyEdges;
         HomologyEdges.push_back(added_edge);
         
         bool from_small_to_big=true;
         int16_t c_pr_edge, c_added_edge;
         uint32_t small,big;
         
         while (true)
         {
            if (d_distance[u]<d_distance[v])
            {
               small = u;
               big = v;
               added_edge = added_edge_v;
               or_edge = or_edge_v;
               added_edge_small = added_edge_u;
               or_edge_small = or_edge_u;
               from_small_to_big = !from_small_to_big;
            }
            else
            {
               small = v;
               big = u;
               added_edge = added_edge_u;
               or_edge = or_edge_u;
               added_edge_small = added_edge_v;
               or_edge_small = or_edge_v;
            }
            
            for (auto nn : vol_mesh->etf(added_edge))
            {
               if (physical_surface[nn] && abs(nn) == big)
               {
                  c_pr_edge    = nn<0 ? -1 : 1;
                  break;
               }
            }
            
            for (auto nn : vol_mesh->etf(d_paredge[big]))
            {
               if (physical_surface[nn] && abs(nn) == big)
               {
                  c_added_edge = nn<0 ? -1 : 1;
                  break;
               }
            }
            
            if (c_added_edge<0)
               or_edge*= c_pr_edge;
            else
               or_edge*=-c_pr_edge;

            added_edge = d_paredge[big];
            int16_t coeff = (d_parent[big]<big)? 1:-1;
            uint32_t s_edge,t_edge;
            
            if (from_small_to_big)
            {
               p1_pushback_mutex.lock();
               h1b[added_edge].push_back(coeff*n_gen);
               p1_pushback_mutex.unlock();
               
               os << vol_mesh->print_edge(n_gen,added_edge,coeff>0,0,255,255);
            }
            else
            {
               p1_pushback_mutex.lock();
               h1b[added_edge].push_back((-coeff)*n_gen);
               p1_pushback_mutex.unlock();
               
               os << vol_mesh->print_edge(n_gen,added_edge,coeff<0,0,255,255);
            }
            
            
            if (d_parent[big] == small)
            {
               
               break;
            }
            
            u=small;
            from_small_to_big = !from_small_to_big;
            added_edge_u=added_edge_small;
            or_edge_u=or_edge_small;
            v=d_parent[big];
            added_edge_v=added_edge;
            or_edge_v=or_edge;
         }
         os.close();
      
      }
      
   }
   return ret;
}
*/

void generic_two_manifold :: TwoSidedAlg(const int16_t& n_gen, const uint32_t& this_edge, const uint32_t& next_edge)
{
   std::ofstream os;
   timecounter t_tc;
   t_tc.tic();
   // auto this_edge = *input.begin();
   // auto next_edge = input.back();
   bool found_tgt = false;
   
   // std::cout << *n_iter << std::endl;
   std::vector<uint32_t> queue;
   // queue.reserve(vol_mesh->volumes_size());
   // std::vector<bool> colour(vol_mesh->volumes_size(),false);
   std::map<uint32_t,uint8_t>      colour; 
   std::map<uint32_t,uint32_t>    par_e, parent;
   
   uint32_t target, source;
   target = source = vol_mesh->volumes_size();
   
   // std::thread s_tr([&] {
   for (auto face : vol_mesh->etf(this_edge))
   {
      auto ff = abs(face);
      auto vols = vol_mesh->ftv(ff);
      for (auto vol : vols)
         if (vol_mesh->is_conductor(abs(vol)))
            if (abs(vol)<source)
               source=abs(vol);
   }
      
      // queue.push_back(source);
   // });
   
   // std::thread t_tr([&] {   
   for (auto face : vol_mesh->etf(next_edge))
   {
      auto ff = abs(face);
      auto vols = vol_mesh->ftv(ff);
      for (auto vol : vols)
         if (vol_mesh->is_conductor(abs(vol))) 
            if (abs(vol)<target)
               target=abs(vol);
   }
   // });
   
   // s_tr.join();
   // t_tr.join();
   t_tc.toc();
   
   // std::cout << "    Build up: " << t_tc << " s" << std::endl;
   // std::cout << "    Sorting: " << t_tc << " s" << std::endl;
   if (target != source)
   {
      t_tc.tic();
      size_t k=0;
      uint32_t qtop;
      queue.push_back(source);
      colour[source]=0;
      
      while (true)
      {
         qtop = queue[k];
         
         // timecounter t_test;
         // t_test.tic();
         for (auto ee : vol_mesh->vtf(qtop))
         {
            if (!physical_surface[abs(ee)])
            {
               for (auto ff : vol_mesh->ftv(abs(ee)))
               {
                  if (abs(ff) != qtop)
                  {
                     if (colour.find(abs(ff)) == colour.end())
                     {               
                        queue.push_back(abs(ff));
                        par_e[abs(ff)]=abs(ee);
                        parent[abs(ff)]=qtop;
                        colour[abs(ff)]=0;
                     }
                     if (abs(ff) == target)
                     {
                        found_tgt=true;
                        break;
                     }
                  }
               }
            }
            
            if (found_tgt)
               break;
         }

         if (found_tgt)
            break;
         
         k++;
         
         if (k>=queue.size())
         {
            std::ofstream os;
            os.open("./dk.txt", std::ofstream::out | std::ofstream::app);
            for (auto j=k-1;k>0;k--)
            {
               os << vol_mesh->print_dual_edge(9,parent[queue[k]],par_e[queue[k]],false,0,255,0);
               os << vol_mesh->print_dual_edge(9,queue[k],par_e[queue[k]],true,0,255,0);
            }
            
            os.close();
            // std::cout << "    boh" << std::endl;
            throw std::invalid_argument("Problem in TS algorithm");
         }
      }

      t_tc.toc();
      // std::cout << "    Breadth first search: " << t_tc << " s" << std::endl;
      
      t_tc.tic();
      auto head = target;
      //~ if (debuggy)
      //~ {
         write_thinned_current.lock();
         os.open("thinned_currents.txt", std::ofstream::out | std::ofstream::app);
         os << vol_mesh->print_dual_edge(n_gen,parent[head],par_e[head],false,255,0,0);
         os << vol_mesh->print_dual_edge(n_gen,head,par_e[head],true,255,0,0);
         os.close();
         write_thinned_current.unlock();
      //~ }
      
      while (parent[head] != source)
      {   
         for (auto vv : vol_mesh->ftv(par_e[head]))
         {
            if (abs(vv) == parent[head])
            {
               // p2_pushback_mutex.lock();
               tc[par_e[head]].push_back(n_gen*vv.Sgn());
               // p2_pushback_mutex.unlock();
            }
            else
            {
               // p2_pushback_mutex.lock();
               tc[par_e[head]].push_back(-n_gen*vv.Sgn());
               // p2_pushback_mutex.unlock();
            }
            break;
         }
      
         head = parent[head];      
         
         //~ if (debuggy)
         //~ {
            write_thinned_current.lock();
            os.open("thinned_currents.txt", std::ofstream::out | std::ofstream::app);
            os << vol_mesh->print_dual_edge(n_gen,parent[head],par_e[head],false,255,0,0);
            os << vol_mesh->print_dual_edge(n_gen,head,par_e[head],true,255,0,0);
            os.close();
            write_thinned_current.unlock();
         //~ }         
      }
      
      for (auto vv : vol_mesh->ftv(par_e[head]))
      {
         if (abs(vv) == parent[head])
         {
            // p2_pushback_mutex.lock();
            tc[par_e[head]].push_back(n_gen*vv.Sgn());
            // p2_pushback_mutex.unlock();
         }
         else
         {
            // p2_pushback_mutex.lock();
            tc[par_e[head]].push_back(-n_gen*vv.Sgn());
            // p2_pushback_mutex.unlock();
         }
         
         break;
      }
      
      t_tc.toc();
      
      // std::cout << "    Coefficient computation: " << t_tc << " s" << std::endl;
   }
   
   if (debuggy)
      os.close();
}
