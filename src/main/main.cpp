/*
 * This source file is part of Topoprocessor.
 *
 * Copyright (C) 2021, Bernard Kapidani bernard(dot)kapidani(at)gmail(dot)com
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * SEE THE GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
 * YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
 * ALONG WITH THIS PROGRAM. IF NOT, SEE <https://www.gnu.org/licenses/>.
 */
//file main.cpp

#include "lean_cohomology.hpp"

// using namespace std;

int main (int argc, char** argv)
{   
    const char *splash =
"    -------------------------------------------------------------------\n"
"    |                      *** TOPOPROCESSOR ***                      |\n"
"    |   Cohomology Groups' Computation for Electromagnetic Modeling   |\n"
"    |      Bernard Kapidani (C) 2021 bernard.kapidani@gmail.com       |\n"
"    |            MNS, Ecole Polytechnique Federale Lausanne           |\n"
"    -------------------------------------------------------------------\n";

   std::cout << splash << std::endl;
   std::cout << std::endl
             << "    TOPOPROCESSOR is free software protected by the GNU General License:"
             << std::endl
             << "    check license.txt file for more info on terms of use"
             << std::endl;
             
   timecounter t_total;
   t_total.tic();
   if (argc <5 )
   {
      std::cout << "    Correct use is: " 
                << "<mesher (netgen or gmsh)> <mesh filename> <conductors filename> <insulators filename>"
                << "<1 | 0 (lean or lazy, default is lean)> <1 | 0 (i.e. minimize cut or not, default is do not minimize)>"
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
   
   std::cout << "    Completed in: " << t_total << " seconds" << std::endl;
   
   return 0;
}
