#ifndef _SCALAPACK_INTERFACE_H_
#define _SCALAPACK_INTERFACE_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi/mpi.h>

#pragma GCC visibility push(default)
int localToGlobal1DZeroStart(int, int, int, int);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
int localToGlobal1DOneStart(int, int, int, int);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
void globalToLocal1DZeroStart(int, int, int, int *, int *);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
void globalToLocal1DOneStart(int, int , int , int *, int *);  
#pragma GCC visibility pop

#pragma GCC visibility push(default)
void transferMatrixToScalapack(double *, double *, int, int, int, int);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
void transferScalapackToMatrix(double *, double *, int, int, int, int);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
void solveEigenvalueScalapack(double *, double *, int);
#pragma GCC visibility pop

#endif
