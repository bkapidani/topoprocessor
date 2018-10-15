#include "lean_cohomology.hpp"

// using namespace std;

int main (int argc, char** argv)
{   
    const char *splash =
"    -------------------------------------------------------------------\n"
"    |                          *** TOPOPRO ***                        |\n"
"    |    Cohomology Basis Computation for Electromagnetic Modeling    |\n"
"    |    Bernard Kapidani (C) 2017 - kapidani.bernard@spes.uniud.it   |\n"
"    |                      DPIA, University of Udine                  |\n"
"    -------------------------------------------------------------------\n";

    std::cout << splash << std::endl;
   timecounter t_total;
   t_total.tic();
   if (argc <5 )
   {
      std::cout << "    Correct use is: <mesher (netgen or gmsh)> <mesh filename> <conductors filename> <insulators filename> <1 | 0 (lean or lazy, default is lean)> <1 | 0 (i.e. minimize cut or not, default is do not minimize)>"
              << std::endl;
      // std::cout << "            or: mesh_filename surface_id"             << std::endl; 
      return EXIT_FAILURE;
   }

   t_total.tic();
   if (argc > 7)
      lean_cohomology test_lean(argv[1],argv[2],argv[3],argv[4],atoi(argv[5])>0,atoi(argv[6])>0,argv[7]);
   else if (argc > 6)
      lean_cohomology test_lean(argv[1],argv[2],argv[3],argv[4],atoi(argv[5])>0,atoi(argv[6])>0,"no");
   else if (argc > 5)
      lean_cohomology test_lean(argv[1],argv[2],argv[3],argv[4],atoi(argv[5])>0,false,"no");
   else
      lean_cohomology test_lean(argv[1],argv[2],argv[3],argv[4],true,false,"no");
   
   t_total.toc();
   
   std::cout << "    The whole mess took: " << t_total << " s" << std::endl;
   
   return 0;
}
