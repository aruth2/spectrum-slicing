#include "sips.h"

#undef __FUNCT__
#define __FUNCT__ "SIPSShowInertias"
PetscErrorCode SIPSShowInertias(EPS eps)
{	
	PetscReal      *shifts;
	PetscInt       *inertias,nShifts;
	PetscInt       i;
    PetscErrorCode ierr;
    	
    ierr = EPSKrylovSchurGetInertias(eps,&nShifts,&shifts,&inertias);  
    PetscPrintf(PETSC_COMM_WORLD,"Subintervals (after setup):\n");
    for (i=0;i<nShifts;i++) 
    {
        PetscPrintf(PETSC_COMM_WORLD," %d - Shift %g  Inertia %d, ",i,shifts[i],inertias[i]);
        if (i)
            PetscPrintf(PETSC_COMM_WORLD," [%g %g] nevals %d ",shifts[i-1],shifts[i],inertias[i]-inertias[i-1]);
        PetscPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
    
    ierr = PetscFree(shifts);CHKERRQ(ierr);
    ierr = PetscFree(inertias);CHKERRQ(ierr);
      
    return(0);
}

#undef __FUNCT__
#define __FUNCT__ "SIPSShowVersions"
PetscErrorCode SIPSShowVersions()
{
    PetscErrorCode ierr;
	  /* Print git versions */
	  /* ------------------- */
	char           petsc_version[PETSC_MAX_PATH_LEN],slepc_version[PETSC_MAX_PATH_LEN];
    ierr = PetscGetVersion(petsc_version, PETSC_MAX_PATH_LEN);CHKERRQ(ierr);
    ierr = SlepcGetVersion(slepc_version, PETSC_MAX_PATH_LEN);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_SELF,"PETSc version:  %s \n", petsc_version); 
    PetscPrintf(PETSC_COMM_SELF,"SLEPc version:  %s \n", slepc_version);

    return(0);
}

#undef __FUNCT__
#define __FUNCT__ "SIPSShowErrors"
PetscErrorCode SIPSShowErrors(EPS eps)
{
    PetscErrorCode ierr;

    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	return(0);
}

#undef __FUNCT__
#define __FUNCT__ "SIPSIntervalOptimization"
PetscErrorCode SIPSIntervalOptimization(EPS eps,PetscReal *interval)
{
        /* This section of the code uses matrix-inertia calculations to optimize the interval. The downside is that at least 3 inertia calculation must be performed before the run.
      *  My tests indicate that the time required to do this was small, a bad interval can significantly slow down the calculation, and so this is an effective strategy when
      * the limits of the eigenspectrum are not known a priori. Better strategies exist when a good guess of the eigenspectrum is known which can also incorporate load balancing.
      * */
    PetscInt ierr;
    PetscReal      *shifts;
    PetscInt        nshifts=1,m,n;
    PetscInt       *inertias;
    //inertias[0]=inertias[1]=0;
    
    Mat A,B;
    ierr = EPSGetOperators(eps,&A,&B);CHKERRQ(ierr);
    ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr); 
    

    double intervalgrowthrate = 0.1,growth;
    do
    { 
    ierr = EPSSetInterval(eps,interval[0],interval[1]);CHKERRQ(ierr);
    ierr = EPSSetUp(eps);CHKERRQ(ierr);
    ierr = EPSKrylovSchurGetInertias(eps,&nshifts,&shifts,&inertias);
//    PetscPrintf(PETSC_COMM_WORLD,"interval %g %g inertias are %d %d\n",interval[0],interval[1],inertias[0],inertias[nshifts-1]);
    growth = (interval[1]-interval[0])*intervalgrowthrate;
    interval[0] -= growth;
    interval[1] += growth;
    }
    while(inertias[nshifts-1]-inertias[0]<m);
    do
    {
        growth = (interval[1]-interval[0])*intervalgrowthrate;
        interval[1] -= growth;
        ierr = EPSSetInterval(eps,interval[0],interval[1]);CHKERRQ(ierr);
        ierr = EPSSetUp(eps);CHKERRQ(ierr);
        ierr = EPSKrylovSchurGetInertias(eps,&nshifts,&shifts,&inertias);
//      PetscPrintf(PETSC_COMM_WORLD,"interval %g %g inertias are %d %d\n",interval[0],interval[1],inertias[0],inertias[nshifts-1]);
    }
    while(inertias[nshifts-1]-inertias[0]>=m);//Go in reverse subtracting from top until top eigenvalue pops out
    interval[1] += growth; //Add back last bit that push it over
    do
    {
        growth = (interval[1]-interval[0])*intervalgrowthrate;
        interval[0] += growth;
        ierr = EPSSetInterval(eps,interval[0],interval[1]);CHKERRQ(ierr);
        ierr = EPSSetUp(eps);CHKERRQ(ierr);
        ierr = EPSKrylovSchurGetInertias(eps,&nshifts,&shifts,&inertias);
//      PetscPrintf(PETSC_COMM_WORLD,"interval %g %g inertias are %d %d\n",interval[0],interval[1],inertias[0],inertias[nshifts-1]);
    }
    while(inertias[nshifts-1]-inertias[0]>=m);//Repeat, this time adding to bottom side
    interval[0] -= growth; //Add back last bit that push it over
    ierr = EPSSetInterval(eps,interval[0],interval[1]);CHKERRQ(ierr);
    
    return 0;
}



#undef __FUNCT__
#define __FUNCT__ "SIPSSubIntervalFromEigenvalues"
PetscErrorCode SIPSSubIntervalFromEigenvalues(EPS eps, double *interval, double *eig_prev, PetscInt nEigPrev)
{
    double epsilon = 1e-4;
    double offset = epsilon * (interval[1]-interval[0]);
    int index;
    MPI_Comm       comm=MPI_COMM_WORLD;
    PetscErrorCode ierr;
    PetscReal      *subint;
    PetscMPIInt    rank,size;
    PetscInt       nintervals,n,m;
    Mat A,B;
    
    ierr = EPSGetOperators(eps,&A,&B);CHKERRQ(ierr);
    ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr); 
    if(nEigPrev == 0)
    nEigPrev = m;
    
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
    nintervals = size;
    ierr = PetscMalloc1(nintervals+1,&subint);CHKERRQ(ierr);
    //printf("number of previous eigenvalues %d\n",nEigPrev);
    
    int i;
    for(i=0;i<=size;i++)
    {
    index =  i * nEigPrev / size;
    //printf("Picking divider %d from previous eigenvalue %d\n",i,index);
    if(i == 0)
    subint[i] = eig_prev[index]-offset;
    else if(i == size)
    subint[i] = eig_prev[index-1] + offset;
    else
    subint[i] = (eig_prev[index-1] + eig_prev[index])/2.0;
    }
    PetscPrintf(PETSC_COMM_WORLD,"Subinterval is \n"); printmatrix(subint,1,size+1,1,1,0);
//    uniformSubIntervals(interval[0],interval[1],subint,rank,size);
    interval[0] = subint[0];
    interval[1] = subint[size];  
    ierr = EPSSetInterval(eps,interval[0],interval[1]);CHKERRQ(ierr);
    ierr = EPSKrylovSchurSetSubintervals(eps,subint); CHKERRQ(ierr);
    ierr = PetscFree(subint);CHKERRQ(ierr);
    return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SIPSSolve"
PetscErrorCode  SIPSolve(Mat *A_pttr, Mat *B_pttr, PetscReal *interval, PetscInt interval_optimization, PetscReal *eig_prev, PetscInt nEigPrev, EPS  *returnedeps)
{
    
    Mat A = *A_pttr;
    Mat B;
    if(B_pttr != NULL)
    B = *B_pttr;	
    else
    B = NULL;
    if (interval == NULL)
    {
    interval = malloc(2*sizeof(double));
    interval[0] = -1.000001;
    interval[1] = 1;
    }
    
    /** Variables used only during solving**/
//    PetscLogStage  stage[4];
    ST             st; 
    KSP            ksp;
    PC             pc;

    PetscInt       nconv,nintervals;
    EPS            eps;
            
    /**Variables resolved or declared in main and SIPSolve**/
    MPI_Comm       comm=MPI_COMM_WORLD;
    PetscErrorCode ierr;
 
    PetscMPIInt    rank,size;
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  
    nintervals = size;
  

    
//  if (show_versions) SIPSShowVersions();
  
    PetscPrintf(PETSC_COMM_SELF,"solving the eigen problem with %d nintervals\n",nintervals);

//    ierr = PetscLogStageRegister("EPSSetUp", &stage[0]);CHKERRQ(ierr);
//    ierr = PetscLogStageRegister("EPSSolve1", &stage[1]);CHKERRQ(ierr);
//    ierr = PetscLogStageRegister("EPSSolve2-end", &stage[2]);CHKERRQ(ierr);
//   ierr = PetscLogStageRegister("Display", &stage[3]);CHKERRQ(ierr);
  

//    ierr = PetscLogStagePush(stage[0]);CHKERRQ(ierr);
  
    /* Setup eigensolver */
    /* Get global interval */
    /*---------------------*/
    ierr = EPSCreate(comm,&eps);CHKERRQ(ierr);
    ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRQ(ierr);
    ierr = EPSKrylovSchurSetPartitions(eps,nintervals);CHKERRQ(ierr);

	//Regular Eigenproblem or generalized eigen problem
    if (B) {
    ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_SELF,"Generalized Eigenvalue Problem\n");}	
    else {
    ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_SELF," Eigenvalue Problem\n");}
    

    /* Set shift-and-invert with Cholesky; select MUMPS if available */
    ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
    ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
    ierr = STSetMatStructure(st,SAME_NONZERO_PATTERN);CHKERRQ(ierr); /* A and B have same sparsity pattern */
    ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
    ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCCHOLESKY);CHKERRQ(ierr);
    ierr = EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE);CHKERRQ(ierr);  
//    ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);CHKERRQ(ierr); 
  

    ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);CHKERRQ(ierr);
    ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr); /* CC: Set operators in global EPS entails a MatRedundant operation */
    ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  

//	ierr = PetscLogStagePop();CHKERRQ(ierr);
//	ierr = PetscLogStagePush(stage[1]);CHKERRQ(ierr);
  

    PetscPrintf(PETSC_COMM_WORLD,"\n  EPSSolve ...\n ------------------- \n"); 
    if(eig_prev == NULL)
    {
        if(!interval_optimization)
        {
        PetscPrintf(PETSC_COMM_WORLD,"Using given interval %g %g\n",interval[0],interval[1]);     
        ierr = EPSSetInterval(eps,interval[0],interval[1]);CHKERRQ(ierr);
        ierr = EPSSetUp(eps);CHKERRQ(ierr);
//    ierr = EPSKrylovSchurGetInertias(eps,&nshifts,&shifts,&inertias);
        }
        else
        {
            PetscPrintf(PETSC_COMM_WORLD,"Optimizing interval\n");     
            SIPSIntervalOptimization(eps,interval);
        }
    }
    else //eig_prev exists, use it to set interval
    {
        PetscPrintf(PETSC_COMM_WORLD,"Retrieving subintervals from previous eigenvalues\n"); 
        SIPSSubIntervalFromEigenvalues(eps, interval, eig_prev,nEigPrev);
    
    }
    PetscPrintf(PETSC_COMM_WORLD,"Interval is %g %g\n",interval[0],interval[1]);
    
    ierr = EPSSolve(eps);CHKERRQ(ierr);
//	ierr = EPSGetConvergedReason(eps,&convergedreason);CHKERRQ(ierr);	
	
    SIPSShowInertias(eps);
          	
	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD," %D eigenvalues converged in [%g, %g]; nconv %D\n",nconv,interval[0],interval[1],nconv);
	
//  ierr = PetscLogStagePop();CHKERRQ(ierr);
//	ierr = PetscLogStagePush(stage[3]);CHKERRQ(ierr);
   
    *returnedeps = eps;
//    ierr = PetscLogStagePop();CHKERRQ(ierr);
   
    EPSView(eps,PETSC_VIEWER_STDOUT_WORLD);

#ifdef CalcDenMat
    Mat D_Mat;
    double weight[nconv];
    for(i=0;i<nconv;i++)
    if(i<nconv/2)
    weight[i]=2;
    else
    weight[i]=0;

    SIPSCreateDensityMat(&eps,weight,nconv,&D_Mat);
    MatView(D_Mat,PETSC_VIEWER_STDOUT_WORLD);  
   
#endif   
   
    return(0);
}

#undef __FUNCT__
#define __FUNCT__ "SIPSCreateDensityMat"
int SIPSCreateDensityMat(EPS *epsaddr, double *weights, int N, Mat *P)
{	
	PetscPrintf(PETSC_COMM_WORLD,"Trying to extract density Mat from converged solution\n");
	PetscInt     nconv,nShifts,*inertias;
    PetscReal    *shifts;

    Vec xr,xi;
	Vec xr_scat;
	VecScatter xctx;
	PetscScalar kr,ki;
	const PetscScalar *xr_arr;
	PetscErrorCode ierr;
	EPS eps = *epsaddr;
    Mat A_Mat, B_Mat, D_Mat;
    PetscInt first_row,last_row;
//  PetscInt m,n;
    PetscInt iVec,col,lrow;
    PetscScalar value;
    
    ierr = EPSGetOperators(eps,&A_Mat,&B_Mat);
    ierr = MatDuplicate(A_Mat,MAT_DO_NOT_COPY_VALUES,&D_Mat);CHKERRQ(ierr);

	ierr = MatCreateVecs(A_Mat,NULL,&xr);CHKERRQ(ierr);
	ierr = MatCreateVecs(A_Mat,NULL,&xi);CHKERRQ(ierr);
	
    ierr = MatGetOwnershipRange(A_Mat,&first_row,&last_row);CHKERRQ(ierr);
	
//	m = last_row-first_row;
//	n = N;
    
	ierr = EPSKrylovSchurGetInertias(eps,&nShifts,&shifts,&inertias);CHKERRQ(ierr);
	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
	
    int iCol, nCols;
    const int *cols;
    int a=0,b=0,c=0;
    ierr = VecScatterCreateToAll(xr,&xctx,&xr_scat);CHKERRQ(ierr);		
	for(iVec=0;iVec<nconv;iVec++)
	{
        a++;
		ierr = EPSGetEigenpair(eps, iVec, &kr, &ki, xr, xi);CHKERRQ(ierr);        
        //ierr = VecView(xr,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

		ierr = VecScatterBegin(xctx,xr,xr_scat,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
		ierr = VecScatterEnd(xctx,xr,xr_scat,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
        //ierr = VecView(xr_scat,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
		ierr = VecGetArrayRead(xr_scat,&xr_arr);CHKERRQ(ierr);		
		
        for(lrow = first_row;lrow<last_row;lrow++)
        {
            b++;
            ierr = MatGetRow(A_Mat,lrow,&nCols,&cols,NULL);CHKERRQ(ierr);		
            for(iCol=0;iCol<nCols;iCol++)
            {
                c++;
                col = cols[iCol];
                value = weights[iVec]*xr_arr[lrow]*xr_arr[col]; 
                ierr = MatSetValues(D_Mat,1,&lrow,1,&col,&value,ADD_VALUES);CHKERRQ(ierr);
            }
            ierr = MatRestoreRow(A_Mat,lrow,&nCols,&cols,NULL);
           
        }
		ierr = VecRestoreArrayRead(xr_scat,&xr_arr);CHKERRQ(ierr);	
	}
    printf("Loop counters %d %d %d\n",a,b,c);
    ierr = MatAssemblyBegin(D_Mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(D_Mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    *P = D_Mat;
	return 0;
}
