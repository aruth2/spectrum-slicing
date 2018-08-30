# spectrum-slicing
Easy to use and well-documented (hopefully soon ;) Spectrum Slicing.

The Shift-and-Invert Parallel Eigenproblem Solver (SIPS) also known as the spectrum-slicing algorithm is a numerical method for solving large eigenvalue problems. The crux of the algorithm is to divide the eigenspectrum into slices which can be solved independently. Because each slice can be solved independently, interprocess communication is kept to a minimum, and efficient performance is maintained up to very large matrices with a very large number of processes in use. An implementation of the spectrum-slicing algorithm demonstrated good performance up to a N=500,000 matrix with n=200,000 processes.[ref] Like many high-performance products, the algorithm features numerous controls, optimization strategies, and difficulties that may arise based on the specifics of the problem.

# Goals of this library:
1. Provide easy-to-use interfaces which work well in most situations.
  Currently, 3 interfaces are provided of gradually increasing difficulty of use (dsygvs, dsygvsx, and SIPSsolve). 
2. Explain the basic performance considerations of the algorithm, limitations, effective strategies, and available controls.

Outline:
*  [Is the Spectrum-Slicing algorithm the right solution to my problem?](#is-the-spectrum-slicing-algorithm-the-right-solution-to-my-problem)
*  [How do I use the Spectrum-Slicing algorithm?](#how-do-i-use-the-spectrum-slicing-algorithm)
  * [dsygvs - Solving a square matrix](#dsygvs---solving-a-square-matrix)
  * [dsygvsx - Expert Driver for square matrices](#dsygvsx---expert-driver-for-square-matrices)
  * [SIPSolve - Solver for PETSc matrices](#sipsolve---solver-for-petsc-matrices)
*  Performance characteristics and theory.


# Is the Spectrum-Slicing algorithm the right solution to my problem?

The spectrum slicing algorithm may be appropriate for your problem if *all* the following conditions are met:

1. [The matrix to be diagonalized is sparse and symmetric](#time-to-solution-versus-sparsity)
1. The matrix is real (At present only routines for real matrices are presented. Routines for complex matrices will be available in a future release, however they will be restricted to full multicommunicator mode which may limit performance in some situations). 

Additionally *at least one* of the following must be true for reasonable performance:
1. [The eigenspectrum is nearly linearly distributed.](*poor-distribution-of-eigenspectrum)
1. A rough idea of the probability distribution of the eigenspectrum is known a priory.
1. The eigenspectrum will be solved repeatedly. 
1. The problem will be solved with few processes.

*Optionally*, the algorithm can make use of a speedup from computing only a portion of the eigenspectrum.

# How do I use the Spectrum-Slicing algorithm?
## dsygvs - Solving a square matrix
## dsygvsx - Expert Driver for square matrices
## SIPSolve - Solver for PETSc matrices

# What are the important knobs and performance characteristics of the solver?
## Time to solution versus sparsity
## Time to solution versus number of processes
## Poor distribution of eigenspectrum
