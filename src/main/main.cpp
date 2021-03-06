#include "lean_cohomology.hpp"

// using namespace std;

int main (int argc, char** argv)
{	

	timecounter t_total;
	t_total.tic();
	if (argc <=5 )
	{
		std::cout << "Correct use is: mesher mesh_filename conductor_id insulator_id true (lean) | false (lazy)" << std::endl;
		// std::cout << "            or: mesh_filename surface_id" 			   << std::endl; 
		return EXIT_FAILURE;
	}

	t_total.tic();
	if (argc > 5)
	{
		lean_cohomology test_lean(argv[1],argv[2],atoi(argv[3]),atoi(argv[4]),atoi(argv[5])>0);
	}
	else
		lean_cohomology test_lean(argv[1],argv[2],atoi(argv[3]),atoi(argv[4]),true);
	
	t_total.toc();
	
	std::cout << "The whole mess took: " << t_total << " s" << std::endl;
	
	return 0;
}
