#include "lean_cohomology.hpp"

// using namespace std;

int main (int argc, char** argv)
{   
    const char *splash =
"    -------------------------------------------------------------------\n"
"    |                          *** TOPOPRO ***                        |\n"
"    |   Cohomology Groups' Computation for Electromagnetic Modeling   |\n"
"    |    Bernard Kapidani (C) 2016 - 2023 bernard.kapidani@epfl.ch    |\n"
"    |                      DPIA, University of Udine                  |\n"
"    |                  ASC, Technische Universitaet Wien              |\n"
"    |            MNS, Ecole Polytechnique Federale Lausanne           |\n"
"    -------------------------------------------------------------------\n";

   std::cout << splash << std::endl;
   timecounter t_total;
   t_total.tic();
   if (argc < 2 )
   {
      std::cout << "    Correct command line argument layout is:" << std::endl;
      std::cout << "    <mesh filename>" 
                << " <mesher (netgen or gmsh, default gmsh)>" 
                << " <conductor labels filename> (default label: 1)"
                << " <insulator labels filename> (default label: 2)"
                << " <1 | 0 (lean or lazy, default is lazy)>"
                << " <1 | 0 (i.e. minimize cut or not, default is do not minimize)>"
                << std::endl;
      // std::cout << "            or: mesh_filename surface_id"             << std::endl; 
      return EXIT_FAILURE;
   }

   t_total.tic();
   if (argc > 7)
      lean_cohomology run_topopro(argv[1],argv[2],argv[3],argv[4],atoi(argv[5])>0,atoi(argv[6])>0,argv[7]);
   else if (argc > 6)
      lean_cohomology run_topopro(argv[1],argv[2],argv[3],argv[4],atoi(argv[5])>0,atoi(argv[6])>0,"no");
   else if (argc > 5)
      lean_cohomology run_topopro(argv[1],argv[2],argv[3],argv[4],atoi(argv[5])>0,false,"no");
   else if (argc > 4)
      lean_cohomology run_topopro(argv[1],argv[2],argv[3],argv[4],false,false,"no");
   else if (argc > 3)
   {
      std::ofstream i("i.txt");
      i << "2" << std::endl;
      i.close();
      lean_cohomology run_topopro(argv[1],argv[2],argv[3],"i.txt",false,false,"no");
      std::remove("i.txt");
   }
   else if (argc > 2)
   {
      std::ofstream i("i.txt");
      std::ofstream c("c.txt");
      i << "2" << std::endl;
      c << "1" << std::endl;
      i.close();
      c.close();
      lean_cohomology run_topopro(argv[1],argv[2],"c.txt","i.txt",false,false,"no");
      std::remove("i.txt");
      std::remove("c.txt");
   }
   else
   {
      std::ofstream i("i.txt");
      std::ofstream c("c.txt");
      i << "2" << std::endl;
      c << "1" << std::endl;
      i.close();
      c.close();
      lean_cohomology run_topopro(argv[1],"gmsh","c.txt","i.txt",false,false,"no");
      std::remove("i.txt");
      std::remove("c.txt");
   }

   t_total.toc();
   
   std::cout << "    Completed in: " << t_total << " seconds" << std::endl << std::endl;
   
   return 0;
}
