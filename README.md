# CoMD: A Molecular Dynamics Proxy Applications Suite

CoMD is one of many _proxy applications_ that support of the [ExMatEx](http://exmatex.lanl.gov) Co-Design Center--one of three Exascale Co-Design Centers. ExMatEx is funded by the DoE office of [Advanced Scientific Computing Research](http://science.energy.gov/ascr) (ASCR). Dr. Karen Pao is the program manager and Dr. William Harrod is the director of the ASCR Research Division. Dr. Tim Germann (Los Alamos) is the ExMatEx Center Director.

## Description of the Code

CoMD is designed to be an extensible molecular dynamics proxy applications suite. The idea is that we will capture the essence of different features of molecular dynamics in this software and extend it to multiple platforms. The code will serve as a vehicle for co-design by allowing others to extend it as needed to get better performance on different architectures.

CoMD allows the user to use either the Lennard-Jones potential or the Embedded Atom Method potential. Both these are discussed in reference (2) below. All we do in the main subroutine in `src-flat/pmd.c` is:

 *  read in a data file
 *  write out a file that can be read by `clsman` (http://www.t12.lanl.gov/home/afv/)
 *  call `do_compute_work()` to compute on this data
 *  write out a final configuration that is `clsman` compatible
 *  free data

All other variants such as `src-viz` and `src-ocl` build off this.


## Quickstart

Building:

    ccmake .
    <enter c c g >
    make

Running (EAM):

    ./comd -p ag -e -f data/8k.inp.gz

Running (LJ):

    ./comd -f data/8k.inp.gz


Command line parameters are:

 * `-h` : print command line parameters
 * `-f <filename>` : the name of the input file
 * `-e` : use EAM potentials
 * `-p <potname>` : name of the EAM potential
 * `-d <potdir>` : directory where EAM potential files reside
 * `-z` : disable periodic boundary conditions

## Directories

`src-flat`:
This is the base code with the data structure flattened out.  We assume that all particles will be mapped to a grid that has cells that are at least one cuttoff in each dimension.  Interactions between particles are handled as interactions between particles in pairs of mesh cells.  `Main` is in `pmd.c`.


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