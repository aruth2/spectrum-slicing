# spectrum-slicing
Easy to use and well-documented Spectrum Slicing.

The Shift-and-Invert Parallel Eigenproblem Solver (SIPS) also known as the spectrum-slicing algorithm is a numerical method for solving large eigenvalue problems. The crux of the algorithm is to divide the eigenspectrum into slices which can be solved independently. Because each slice can be solved independently, interprocess communication is kept to a minimum, and efficient performance is maintained up to very large matrices with a very large number of processes in use. An implementation of the spectrum-slicing algorithm demonstrated good performance up to a N=500,000 matrix with n=200,000 processes.[ref] Like many high-performance products, the algorithm features numerous controls, optimization strategies, and difficulties that may arise based on the specifics of the problem.

# Goals of this library:
1. Provide easy-to-use interfaces which work well in most situations.
  Currently, 3 interfaces are provided of gradually increasing difficulty of use (dsygvs, dsygvsx, and SIPSsolve). 
2. Explain the basic performance considerations of the algorithm, limitations, effective strategies, and available controls.

Outline:
*  [Is the Spectrum-Slicing algorithm the right solution to my problem?](#is-the-spectrum-slicing-algorithm-the-right-solution-to-my-problem)
*  How do I use it?
*  Performance characteristics and theory.


# Is the Spectrum-Slicing algorithm the right solution to my problem?

# How do I use the Spectrum-Slicing algorithm?

# What are the important knobs and performance characteristics of the solver?
