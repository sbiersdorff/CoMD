# Major Update (November 6 2012)

The new version includes various updates. Major additions are the ability to initialize a FCC lattice with a deformation gradient, Chebyshev evaluation of EAM tables, computation of virial stress and variable box factor. See below for more detail.

# CoMD: A Molecular Dynamics Proxy Applications Suite

CoMD is one of many _proxy applications_ that support of the [ExMatEx](http://exmatex.lanl.gov) Co-Design Center--one of three Exascale Co-Design Centers. ExMatEx is funded by the DoE office of [Advanced Scientific Computing Research](http://science.energy.gov/ascr) (ASCR). Dr. Karen Pao is the program manager and Dr. William Harrod is the director of the ASCR Research Division. Dr. Tim Germann (Los Alamos) is the ExMatEx Center Director.

## Description of the Code

CoMD is designed to be an extensible molecular dynamics proxy applications suite. The idea is that we will capture the essence of different features of molecular dynamics in this software and extend it to multiple platforms. The code will serve as a vehicle for co-design by allowing others to extend it as needed to get better performance on different architectures.

CoMD allows the user to use either the Lennard-Jones potential or the Embedded Atom Method potential. Both these are discussed in reference (2) below. All we do in the main subroutine in `src-flat/pmd.c` is:

 *  read a data file
 *  write out a file that can be read by `clsman` (http://www.t12.lanl.gov/home/afv/)
 *  call `do_compute_work()` to compute on this data
 *  write out a final configuration that is `clsman` compatible
 *  free data

All other variants such as `src-viz` and `src-ocl` build off this.

The default is now to generate an FCC lattice, of size 20x20x20 lattice constants, for a total of 32,000 atoms. The lattice is automatically scaled according to the potential chosen. The atoms are assigned to cubic link cells, the size of which can be controlled by the '-b' box-factor parameter. Increasing the box factor forms larger link cells with more partcles, which may improve performance on some architectures.

When EAM potentials are used, the code generates a set of corresponding Chebyshev coefficients from the table data. The choice of using the EAM table lookup or the Chebyshev reconstruction is controlled by the `USE_CHEBY` flag in the cmake configuration step.

## Quickstart

Building:

Create a 'build' directory at the top level.
   
    mkdir build
    cd build
    ccmake ..
    <enter c c g >
    make

Several options are available in ccmake, including:

   * toggle single/double precision
   * toggle OpenCL build
   * toggle in-situ viz
   * toggle VTK use for viz
   * toggle viz interop (for OpenCL version)
   * toggle approximate centrosymmetry computation in viz
   * toggle the use of Chebyshev coefficients to compute the EAM values

Running (LJ):

    ./CoMD 

Running (EAM):

    ./CoMD -e

Running OpenCL

    ./CoMDOCL (all options available in base version)

Command line parameters are:

 * `-h` : print command line parameters
 * `-f <filename>` : the name of the input file
 * `-e` : use EAM potentials
 * `-p <potname>` : name of the EAM potential
 * `-z` : disable periodic boundary conditions
 * `-x <Nx> -y <Ny> -z <Nz>` : create an FCC lattice with the given dimensions (in lattice cell dimensions).
 * `-s <defgrad>` : perform a 1D stretch (in x) of the domain with the specified deformation gradient.
 * `-b <boxfactor>` : set the minimum box factor for the link cells. Must be > 1.0

Additional command line parameters for the OpenCL version:

 * `-g` : use the GPU if available

## Directories

`src-flat`:
This is the base code with the data structure flattened out.  We assume that all particles will be mapped to a grid that has cells that are at least one cuttoff in each dimension.  Interactions between particles are handled as interactions between particles in pairs of mesh cells.  `Main` is in `pmd.c`.

`src-ocl`:
This contains the OpenCL implementation. The OpenCL version uses the IO routines in `src-flat` and makes copies of that data (with different data layout, as appropriate). The OpenCL code runs the base code after the OpenCL runs, to test for correctness. 

`data`:
This contains a few input files for testing.  Note that the code in the `src-flat` directory has been verified against `clsman` and we have confidence that it gives correct results for EAM potential.

_Note: the LJ potential has not been verified on `src-flat`_.

`pots`: 
This directory contains potentials for the EAM potential. These potentials are described in a few different papers, beginning with _Ni_ and _Al_ in Voter and Chen, 1987 [1]. The most complete description of the fitting procedure is in a chapter in a book on intermetallics [2].  This complete description can also be found in a Los Alamos technical report, along with the parameters for all seven _fcc_ potentials [3].  This technical report and the original Voter-Chen paper are now included with this potentials package as pdf files. For the _cu_ potential, please cite Ref. [4], which gives the parameters and cites the earlier papers (or, cite both Ref. [4] and Ref. [2]). For the silver potential (_ag_), please cite Refs. [5] and [2]. For _Pd_, _Pt_, and _Au_, simply cite reference [3]. The _cu1_ and _ag1_ potentials are designed to work with each other for _Cu-Ag_ alloy simulations.  The cross potential is simply fit to the heat of solution at infinite dilution for _Cu_ in _Ag_ and _Ag_ in _Cu_. As pure elements (_cu1_, _ag1_), they work exactly like the _cu_ and _ag_ pots.

[1] A.F. Voter and S.P. Chen, Mat. Res. Soc. Symp. Proc. Vol. 82, 175 (1987). _"Accurate Interatomic Potentials for Ni, Al, and Ni3Al"_

[2] A.F. Voter, in Intermetallic Compounds: Principles and Practice, edited by J.H. Westbrook and R.L. Fleischer, (John Wiley and Sons, Ltd, 1995), Vol. 1, p. 77. _"The Embedded Atom Method"_

[3] A.F. Voter, Los Alamos Unclassified Technical Report LA-UR-93-3901. _"Embedded Atom Method Potentials for Seven FCC Metals: Ni, Pd, Pt, Cu, Ag, Au, and Al"_

[4] A.F. Voter, Phys. Rev. B, 57, 13985 (1998). _"Parallel Replica Method for Dynamics of Infrequent Events"_

[5] A.F. Voter, in Modeling of Optical Thin Films, M.R. Jacobson, Ed., Proc. SPIE Vol. 821, 214 (1987). _"Simulation of the Layer-Growth Dynamics in Silver Films: Dynamics of Adatom Clusters and Vacancy Clusters on Ag(100)"_

## Additional Information

A good introduction to molecular dynamics can be found in _"Computer Simulation of Liquids"_ by M.P. Allen and D.J. Tildesley.

A description of how to run MD at large scales can be found in _"369 Tflop/s Molecular Dynamics Simulations on the Roadrunner General-Purpose Heterogeneous Supercomputer"_ by S. Swaminarayan, K. Kadau, T.C. Germann, and G. Fossum published in Proceedings of SC'08.

## Requirements

There are a number of features that can be configured to be on or off in the Makefiles.  CMake is of course needed to build anything, but you can build the non-OpenCL code without OpenCL, the non-viz code without GLUT, and the viz code with or without VTK.  To build everything, you would want CMake, OpenCL, GLUT, and VTK.

Most features have been tested on Mac and on Linux (Fedora).

## Authors

The initial authors of CoMD (all at Los Alamos National Laboratory) are:

 * Sriram Swainarayan
 * Jamaludin Mohd-Yusof
 * Christopher Sewell

CoMD is released as LA-CC-11-119.
