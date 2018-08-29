#ifndef _SIPS_SQUARE_H_
#define _SIPS_SQUARE_H_

#include <stdio.h>
#include <stdlib.h>
#include "sips.h"

#define SIPSReal

#pragma GCC visibility push(default)
int dsyevs(int *, double  *, double *);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
int dsygvs(int *, double  *, double *, double *);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
int dsygvsx( int 	*, double  *, double  *, double 	*, PetscInt, double *, int, char 	*, char *);
#pragma GCC visibility pop


#pragma GCC visibility push(default)
int squareFromEPS(EPS *,PetscScalar *,double *,int);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
int eigenvaluesFromEPS(EPS *,double *,int);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
int eigenvectorsFromEPS(EPS *,PetscScalar *,int);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
Mat squareToPetsc(PetscScalar *, int);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
void printmatrix(double *, int, int,int,int,int);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
void printdatatypes();
#pragma GCC visibility pop

#pragma GCC visibility push(default)
void UTtoSym(PetscScalar *, int, int);
#pragma GCC visibility pop

#endif
