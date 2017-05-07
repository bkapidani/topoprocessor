#include "lean_cohomology.hpp"

// using namespace std;

int main (int argc, char** argv)
{	

	timecounter t_total;
	t_total.tic();
	if (argc !=4 )
	{
		std::cout << "Correct use is: mesh_filename conductor_id insulator_id" << std::endl;
		// std::cout << "            or: mesh_filename surface_id" 			   << std::endl; 
		return EXIT_FAILURE;
	}

	t_total.tic();
	lean_cohomology test_lean(argv[1],atoi(argv[2]),atoi(argv[3]));
	t_total.toc();
	
	std::cout << "The whole mess took: " << t_total << " s" << std::endl;
	
	return 0;
}
