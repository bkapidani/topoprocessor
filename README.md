# Documentation #

This is a brief vademecum for Topoprocessor, a C++ based library which computes the cohomology generators of combinatorial 3-manifolds with boundaries. The present version of the tool was born out of a collaborative effort between Pawel Dłotko, Bernard Kapidani (main maintainer and developer of this repository), and Ruben Specogna. 

### What is this repository for? ###

If you are numerically solving electromagnetic eddy currents problems with a scalar magnetic potential formulation, or flow problems with a stream-function formulation, the problem of computing the cohomology generators of some domain of interest embedded in $R^3$ arises. The present tool supports the domain of interest to be discretized through a mesh (simplices or hexahedra) and is compatible with outputs provided by open source meshers NETGEN (visit www.ngsolve.org) in the neutral mesh format and GMSH (visit www.gmsh.info) in legacy 2.2 msh format. 
The whole computational domain is usually a topologically trivial one, $E = A \cup B$ where $A$ is the nontrivial subset for which the cohomology group needs to be generated. In accordance with the main application of the toolbox (discrete formulations for eddy current problems in electromagnetism) domain $A$ is called the insulator and domain $B$ is called the conductor.

### How do I install Topoprocessor? ###

To install the tool open a terminal in a folder owned by the current user run the following:

    TOPODIR=./topoprocessor
    mkdir $TOPODIR/topo-build
    mkdir $TOPODIR/topo-install
    cd $TOPODIR
    git clone https://github.com/bkapidani/topoprocessor.git topo-src
    cd topo-build
    cmake -DCMAKE_INSTALL_PREFIX=../topo-install ../topo-src/src
    make -j
    make install

where the TOPODIR variable can be set to any specific folder.
If you have admin privileges and want to install Topoprocessor for all users, just omit the CMAKE_INSTALL_PREFIX flag and run
    
    sudo make install

instead of the last instruction above.

The tool is in principle cross-platform, although developed on UNIX based architectures and originally adapted to compile it on Microsoft Windows through the POSIX compatibility layer Cygwin.
Recently, the introduction of WSL on Windows architectures should have made the compilation of the library seamless on Microsoft Windows too, but this has not been officially tested by the developers. Feel free to contact them if you have done so yourself and want to share your successful build-chain.

### How do I run the executable? ###

Topoprocessor takes a variable number of command line inputs, the minimal input consists in specifying the mesher, the mesh filename and the files containing material labels for the conductor and insulator domains. These last two files are very easy to produce and are infact just text files containing an integer on the first line, indicating the number of labels, followed by all the labels that have to be regarded as conductors or insulators, respectively. An example of calling Topoprocessor is then as follows:

    topoprocessor gmsh mymesh.msh conductors.txt insulators.txt

where we note that the mesh format compatible for gmsh meshes is the 2.2 version, while for NETGEN the neutral mesh format is recognized. In general the syntax for calling Topoprocessor can be recalled by running the executable without input arguments, which will yield the following message:

    topoprocessor  <mesher ("netgen" or "gmsh")>  <mesh filename>  <conductors filename> <insulators filename> <1 | 0 (lean or lazy, default is lean)> <1 | 0 (i.e. minimize cut or not, default is do not minimize)>

where the optional arguments are used to choose the lazy version of the algorithm (which doubles the number of generators, see references below), and to choose if to run a post-processing step which minimizes the support of the generators, based on the solution of the discrete Plateau problem (see again references below).

### What do I get as an output? ###

Topoprocessor returns a text file, called h1.txt which looks as follows: The first line is the number of ($H^1$) cohomology generators found. Then, for each generator you will have a line indicating the size N_i of the support of the generator (the number of nonzero coefficients) and then you will have exactly N_i lines with 3 integers, the first two are the mesh node indices which are the end-points of a given edge (the labels match the indices in the mesh file you provided as an input) and the third number is the coefficient (a signed integer) associated the edge for the given generator.

### New: structured grids ###

The algorithms now work also for structured grids, the support is still experimental, and has been used recently in high-order splines methods, with MATLAB scripts to provide the results as inputs for the the GeoPDEs library (https://github.com/rafavzqz/geopdes).
To access the version for structured grids, which can be used to combine the algorithms with spline bases in GeoPDEs, one needs to check-out to the separate branch. The plan is in the future whole choice to be parametrised on the main branch. Since this version is still experimental consider contacting the developer before use.

### How do I cite the proper sources? ###

If you find Topoprocessor useful and you use it in your research please cite the relevant papers. It will motivate the developers, who are researchers themselves, to improve the algorithms and maintain the code as bug-free as possible. For the lazy and lean version of the main code cite:

[1] P. Dłotko and R. Specogna, ‘Lazy Cohomology Generators: A Breakthrough in (Co)homology Computations for CEM’, IEEE Transactions on Magnetics, vol. 50, no. 2, pp. 577–580, Feb. 2014, doi: 10.1109/TMAG.2013.2281076.

[2] P. Dłotko, B. Kapidani, and R. Specogna, ‘Topoprocessor: An Efficient Computational Topology Toolbox for h-Oriented Eddy Current Formulations’, IEEE Transactions on Magnetics, vol. 53, no. 6, pp. 1–4, Jun. 2017, doi: 10.1109/TMAG.2017.2661480.

[3] P. Dłotko, B. Kapidani, and R. Specogna, ‘Lean Cohomology Computation for Electromagnetic Modeling’, IEEE Transactions on Magnetics, vol. 54, no. 3, pp. 1–4, Mar. 2018, doi: 10.1109/TMAG.2017.2749618.

If you use the support minimizing procedure to reduce the fill-in of the final system please additionally cite:

[4] P. Dłotko, B. Kapidani, and R. Specogna, ‘Fast Computation of Cuts With Reduced Support by Solving Maximum Circulation Problems’, IEEE Transactions on Magnetics, vol. 51, no. 3, pp. 1–4, Mar. 2015, doi: 10.1109/TMAG.2014.2359976.

If you use the structured version, please cite:

[5] B. Kapidani, M. Merkel, S. Schöps, and R. Vázquez, ‘Tree–cotree decomposition of isogeometric mortared spaces in H(curl) on multi-patch domains’, Computer Methods in Applied Mechanics and Engineering, vol. 395, p. 114949, May 2022, doi: 10.1016/j.cma.2022.114949.


### All this looks useful, but scares me, who do I talk to? ###

If you need further infos, or you want to report bugs, the main developer/mantainer of the Topoprocessor tool can be reached at:

bernard(dot)kapidani(at)gmail(dot)com
