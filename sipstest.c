#include "sips.h"

static char help[] = "Generalized Symmetric eigenproblem: Ax = lambda Bx \n\n"
  "The command line options are:\n"
  "  -fA <binary data file for A>\n"
  "  -fB <binary data file for B> \n\n";
/* 
   Example:
   mpiexec -n 8 ./dftb2 -fA $D/graphene_xxs_A -fB $D/graphene_xxs_B -nintervals 4 (-eps_krylovschur_partitions)
   mpiexec -n 2 ./dftb2 -fA $D/small -show_errors -interval 0.0,4.01
*/

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
 
  SlepcInitialize(&argc,&argv,"dftbOption",help); 
  /** Variables used only by main **/
  PetscBool fileAFlag, fileBFlag, commandLineFlag;
  PetscViewer    fd;              
  char           file[2][PETSC_MAX_PATH_LEN];     /* input file name */  
  
  
  /**Variables passed by main to SIPSolve **/
  Mat            A,B;        
  PetscBool      show_inertias=PETSC_FALSE,show_errors=PETSC_FALSE,show_versions=PETSC_FALSE;
  PetscInt       nintervals,nmax,maxits; //nmax appears to be related to nintervals, but nmax is never used in the solver
  PetscReal      interval[2];
  
  
  /**Variables resolved or declared in main and SIPSolve**/
  MPI_Comm       comm=MPI_COMM_WORLD;
  PetscErrorCode ierr;
 
  PetscMPIInt    rank,size;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  
  PetscInt       n,m;
  //ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr); 
  
  PetscBool      headProcess = (rank == 0);
  
//  PetscLogDouble time0,time1;

 
  
  /**Default values set by main**/
  nintervals  = size;
  interval[0] = -0.8; 
  interval[1] = 0.2;
  nmax        = 2;
  maxits = 2;
  
  
  ierr = PetscOptionsGetBool(NULL,NULL,"-show_versions",&show_versions,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-show_inertias",&show_inertias,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-show_errors",&show_errors,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-nintervals",&nintervals,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetRealArray(NULL,NULL,"-interval",interval,&nmax,&commandLineFlag);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-sips_maxits",&maxits,NULL);CHKERRQ(ierr);


  /* Load matrices A and B over MPI_COMM_WORLD */
  /* ----------------------------------------- */

//  ierr = PetscLogStagePush(stage[0]);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-fA",file[0],PETSC_MAX_PATH_LEN,&fileAFlag);CHKERRQ(ierr);
  if (!fileAFlag) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Must indicate binary file with the -fA options");
  ierr = PetscViewerBinaryOpen(comm,file[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = MatCreate(comm,&A);CHKERRQ(ierr);
  ierr = MatSetType(A,MATSBAIJ);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatLoad(A,fd);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr); 
  if (headProcess) 
    ierr = PetscPrintf(PETSC_COMM_SELF,"\nA: %d %d\n",m,n);CHKERRQ(ierr); 
  

  ierr = PetscOptionsGetString(NULL,NULL,"-fB",file[1],PETSC_MAX_PATH_LEN,&fileBFlag);CHKERRQ(ierr); 
  if (fileBFlag) {
    ierr = PetscViewerBinaryOpen(comm,file[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = MatCreate(comm,&B);CHKERRQ(ierr);
    ierr = MatSetType(B,MATSBAIJ);CHKERRQ(ierr);
    ierr = MatSetFromOptions(B);CHKERRQ(ierr);
    ierr = MatLoad(B,fd);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
    ierr = MatGetSize(B,&m,&n);CHKERRQ(ierr);
    if (headProcess) printf("B: %d %d\n",m,n);
  } else {
    B = NULL;
    if (headProcess) printf("Standard eigenvalue problem \n");
  }
  
  //ierr = PetscLogStagePop();CHKERRQ(ierr);



  EPS eps; 
  SIPSolve(&A, &B, interval,&eps);
   
  /* Free objects */
//  ierr = PetscLogStagePush(stage[5]);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
 
  //ierr = PetscLogStagePop();CHKERRQ(ierr);
  
  ierr = SlepcFinalize();
  return 0;
}
