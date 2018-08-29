#ifndef _SIPS_H_
#define _SIPS_H_

#include <petsctime.h>
#include <slepceps.h>

#define verbose 0
#define printf if(verbose) printf
#define PetscPrintf if(verbose) PetscPrintf
#define EPSView if(verbose) EPSView
#define MatView if(verbose) MatView
#define VecView if(verbose) VecView

//#define CalcDenMat


extern PetscErrorCode SIPSShowInertias(EPS);
extern PetscErrorCode SIPSShowVersions();
extern PetscErrorCode SIPSShowErrors(EPS);

#pragma GCC visibility push(default)
extern PetscErrorCode SIPSolve(Mat *,Mat *,PetscReal *,PetscInt , PetscReal *, PetscInt, EPS *);
#pragma GCC visibility pop

#pragma GCC visibility push(default)
extern PetscErrorCode SIPSCreateDensityMat(EPS *,PetscReal*, PetscInt, Mat*);
#pragma GCC visibility pop

#endif
