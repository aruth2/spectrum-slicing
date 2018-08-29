#include "sips_square.h"

void printdatatypes()
{
	printf("Sizes: PetscScalar: %ld, PetscReal: %ld, Float(s): %ld, Double(d): %ld, Complex(c): %ld Double Complex(z): %ld\n",sizeof(PetscScalar),sizeof(PetscReal),sizeof(float),sizeof(double),sizeof(float complex),sizeof(double complex)); 
}

void printmatrix(double *matrix, int rows, int cols, int stride, int rowmajor, int brief)
{
	printf("Printing a matrix with dimensions of %d X %d using a stride of %d brevity %d\n",rows,cols,stride,brief);
	/****************************************************/
	/* This function outputs the contents of a matrix in*/
	/* nice neat format.								*/
	/****************************************************/
	int row,col;
	double value;
    int nlimit = 4,lastprint=1;
	for(row=0;row<rows;row++)
    for(col=0;col<cols;col++)
    {
        if(brief)
        {
            if((col >= nlimit) && (col < cols-nlimit))
            {
            if(lastprint)
            {
                printf("...\t");
                lastprint = 0;
            }
            continue;
            }
            if((row >= nlimit) && (row < rows-nlimit))
            {
                if(lastprint)
                {
                    printf("...\n");
                    lastprint = 0;
                }
                continue;    
            }
        }
        if(rowmajor)
            value = *(matrix+(row*cols+col)*stride);
        else
            value = *(matrix+(col*rows+row)*stride);
        printf("%.03f \t",value);
			
        if(col == cols-1)
            printf("\n");

        lastprint = 1;
	}
}

double listmax(double *list, int num)
{
	double max = *list;
	int index;
	for(index = 1;index<num;index++)
        if(*(list+index) > max)
            max = *(list+index);
	return max;
}

double listmin(double *list, int num)
{
	double min = *list;
	int index;
	for(index = 1;index<num;index++)
        if(*(list+index) < min)
            min = *(list+index);
	return min;
}

void UTtoSym(PetscScalar *matrix, int rows, int columns)
{
	/* Given an upper triangular matrix, this fills in the lower triangle to make it symmetric
	 * */
	if(matrix == NULL)
        return;
    
    int row,column;
	
	for(row=1;row<rows;row++)
        for(column=0;column<row;column++)
            *(matrix + row*columns+column) = *(matrix + column*rows+row);
}

void LTtoSym(PetscScalar *matrix, int rows, int columns)
{
	/* Given a lower triangular matrix, this fills in the lower triangle to make it symmetric
	 * */
	if(matrix == NULL)
        return;
    int row,column;
	
	for(row=1;row<rows;row++)
	for(column=0;column<row;column++)
        *(matrix + column*rows+row) = *(matrix + row*columns+column);
}

Mat squareToPetsc(PetscScalar *A, int N)
{
    if (A == NULL)
        return NULL;
	printf("Making a square matrix of dimensions %d\n",N);
	Mat A_Mat;
	
	MPI_Comm       comm=MPI_COMM_WORLD;
	PetscErrorCode ierr;

	PetscInt       n,m,first_row,last_row, i,j;
    
     
	PetscMPIInt    rank,size;
	ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
	ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

	
	ierr = MatCreate(comm,&A_Mat);CHKERRQ(ierr);
	
	ierr = MatSetType(A_Mat,MATAIJ);CHKERRQ(ierr);
	ierr = MatSetSizes(A_Mat,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A_Mat);CHKERRQ(ierr);    
	ierr = MatSetUp(A_Mat);CHKERRQ(ierr);
    //Implicit assumption that matrix is symmetric
	ierr = MatSetOption(A_Mat,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);

	ierr = MatGetOwnershipRange(A_Mat,&first_row,&last_row);CHKERRQ(ierr);
	
	m = last_row-first_row;
	n = N;
	printf("[%d] M %d N %d m %d n %d firstrow %d lastrow %d\n",rank,N,N,m,n,first_row,last_row);
    int nonzerocounter=0;
    
	for(i=first_row;i<last_row;i++)
        for(j=0;j<n;j++)
            if (isnormal(*(A+N*i+j)))
            {
                //printf("Attempting to insert values i: %d %lx j: %d %lx matrix: %lx\n",i,&i,j,&j,A+i*N+j);
                ierr = MatSetValues(A_Mat,1,&i,1,&j,A+i*N+j,INSERT_VALUES);CHKERRQ(ierr);
                nonzerocounter++;
            }

	printf("%d values in matrix %g per row\n",nonzerocounter,(double)nonzerocounter/N);
    
	ierr = MatAssemblyBegin(A_Mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A_Mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	//MatView(A_Mat,PETSC_VIEWER_STDOUT_(comm));

	return A_Mat;
}

int squareFromEPS(EPS *epsaddr,PetscScalar *A,double *W,int N)
{	
    EPS eps = *epsaddr;
	PetscInt     nconv;
	//Eigenvectors x, and eigenvalues k on exit from sipssolve
	//These eigenvectors/values are distributed accross the nodes
    Vec xr,xi;
	Vec xr_scat;
	PetscScalar kr,ki;
	const PetscScalar *xr_arr;
	PetscErrorCode ierr;
	//Mat A_Mat = *Aaddr;
    Mat A_Mat, B_Mat;
    ierr = EPSGetOperators(eps,&A_Mat,&B_Mat);

	ierr = MatCreateVecs(A_Mat,NULL,&xr);CHKERRQ(ierr);
	ierr = MatCreateVecs(A_Mat,NULL,&xi);CHKERRQ(ierr);
	
	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
	
	printf("Copying %d eigenvalues to square matrix\n",nconv);
	VecScatter ctx;
	PetscInt row,col;
	
	for(row=0;row<N;row++)
	for(col=0;col<N;col++)
	*(A+row*N+col) = 0; //Set matrix to zero before copying over eigenvectors
	
		
	for(row=0;row<nconv;row++)
	{
		ierr = EPSGetEigenpair(eps, row, &kr, &ki, xr, xi);CHKERRQ(ierr);
		W[row] = kr;
		//ierr = VecView(xr,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		ierr = VecScatterCreateToAll(xr,&ctx,&xr_scat);CHKERRQ(ierr);
		ierr = VecScatterBegin(ctx,xr,xr_scat,ADD_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
		ierr = VecScatterEnd(ctx,xr,xr_scat,ADD_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
		//ierr = VecView(xr_scat,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

		//Loading from the scattered vector
		ierr = VecGetArrayRead(xr_scat,&xr_arr);CHKERRQ(ierr);		
		for(col=0;col<N;col++)
		{
			*(A+row*N+col) = xr_arr[col];
		}
		ierr = VecRestoreArrayRead(xr_scat,&xr_arr);CHKERRQ(ierr);	
	}
	return 0;
}


int eigenvectorsFromEPS(EPS *epsaddr,PetscScalar *A,int N)
{	
	printf("Trying to extract eigenvectors from converged solution\n");
	PetscInt     nconv;
	//Eigenvectors x, and eigenvalues k on exit from sipssolve
	//These eigenvectors/values are distributed accross the nodes
    Vec xr,xi;
	Vec xr_scat;
	const PetscScalar *xr_arr;
	PetscErrorCode ierr;
	EPS eps = *epsaddr;
	//Mat A_Mat = *Aaddr;
    Mat A_Mat, B_Mat;
    ierr = EPSGetOperators(eps,&A_Mat,&B_Mat);

	ierr = MatCreateVecs(A_Mat,NULL,&xr);CHKERRQ(ierr);
	ierr = MatCreateVecs(A_Mat,NULL,&xi);CHKERRQ(ierr);
	
	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
    //ierr = SIPSShowInertias(eps);CHKERRQ(ierr);
	
	printf("Copying %d eigenvectors to square matrix\n",nconv);
	VecScatter ctx;
	PetscInt row,col;
	
	for(row=0;row<N;row++)
        for(col=0;col<N;col++)
            *(A+row*N+col) = 0; //Set matrix to zero before copying over eigenvectors
	
		
	for(row=0;row<nconv;row++)
	{
		ierr = EPSGetEigenvector(eps, row, xr, xi);CHKERRQ(ierr);
		//ierr = VecView(xr,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        		
		ierr = VecScatterCreateToAll(xr,&ctx,&xr_scat);CHKERRQ(ierr);
		ierr = VecScatterBegin(ctx,xr,xr_scat,ADD_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
		ierr = VecScatterEnd(ctx,xr,xr_scat,ADD_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
		//ierr = VecView(xr_scat,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
		ierr = VecGetArrayRead(xr_scat,&xr_arr);CHKERRQ(ierr);		
		for(col=0;col<N;col++)
		{
			*(A+row*N+col) = xr_arr[col];
		}
		ierr = VecRestoreArrayRead(xr_scat,&xr_arr);CHKERRQ(ierr);	
	}
	return 0;
}


int eigenvaluesFromEPS(EPS *epsaddr,double *W,int N)
{	
	printf("Trying to extract eigenvalues from converged solution\n");
	PetscInt     nconv;
	PetscScalar kr,ki;
	PetscErrorCode ierr;
	EPS eps = *epsaddr;
	
	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
    
	PetscInt row;		
		
	for(row=0;row<nconv;row++)
	{
		ierr = EPSGetEigenvalue(eps, row, &kr, &ki);CHKERRQ(ierr);
		W[row] = kr;
	}
    return 0;
}

int dsyevs(int 	*N, double  *A, double 	*W)
{
    //SIPSolve accepts B_Mat = NULL, the simplest call is to 
    return dsygvs(N,A,NULL,W);
}

int dsygvs(int 	*N, double  *A, double  *B, double 	*W)
{
	//Double symmetric generalized EigenValue/Vector Sips
	//This function converts the passed matrices A and B into Petsc matrices and then calls SIPSolve to solve the system
    //Variables similar to LAPACK have retained their names
    //JOBZ - 'V' or 'v' compute vectors, 'N' or 'n' do not compute vectors
    //N - the dimension of the matrix
    //A - assumed to be a square symmetric matrix in C-style row-major order, or Fortran-style column-major order.
    // Assuming the symmetry, and that both upper and lower triangles are filled, loading the Petsc matrix will proceed assuming occupied rows 
    //On exit, if JOBZ == 'V' or 'v' A contains the eigenvectors
    //LDA - leading dimensions of matrix A
    //W - On exit, the eigenvalues of A
    //interval - the interval overwhich the problem is solved interval[0] - lower bound, interval[1] upperbound
    //nintervals - the number of subintervals over which the problem will be divided.
    PetscBool isInitialized;
  	EPS eps;
	PetscErrorCode ierr;
    
    SlepcInitialized(&isInitialized);
    if(!isInitialized)
    SlepcInitialize(NULL,NULL,NULL,NULL); 	
    
    //printf("Double precision generalized eigenvalue problem dimensions NXLDAxLDB: %d %d %d\n Matrix A: %lx matrix B: %lx\n",*N,*LDA,*LDB,*Aprime,*Bprime);
   

    Mat A_Mat = squareToPetsc(A,*N);
    Mat B_Mat = squareToPetsc(B,*N);
	
	SIPSolve(&A_Mat, &B_Mat, NULL,1,NULL,NULL,&eps);
	squareFromEPS(&eps,A,W,*N);	
    
    //This code is for testing the efficiency of using the eigenvalues to set the interval. The code should be timed in parallel with and without this second diagonalization call 
    /*
    EPSDestroy(&eps);
    SIPSolve(&A_Mat, &B_Mat, NULL,1,W,&eps);//Repeat a second time, this time with the previous eigenvalues
    squareFromEPS(&eps,A,W,*N);	
*/
    
    
    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
    ierr = MatDestroy(&A_Mat);CHKERRQ(ierr);
    ierr = MatDestroy(&B_Mat);CHKERRQ(ierr);
    	
	return 0;
}

int dsygvsx( int 	*N, double  *A, double  *B, double 	*W, PetscInt interval_optimization, double *interval, int prev_eigenvalues, char 	*JOBZ, char *UPLO)
{
	//Expert Driver for Double symmetric generalized EigenValue/Vector Sips
	//This function converts the passed matrices A and B into Petsc matrices and then calls SIPSolve to solve the system
    //Variables similar to LAPACK have retained their names
    //N - On input, the dimension of the matrix, on exit the number of converged eigenvalues if intervalOptimization was not used.
    //A - assumed to be a square symmetric matrix in C-style row-major order, or Fortran-style column-major order.
    // Assuming the symmetry, and that both upper and lower triangles are filled, loading the Petsc matrix will proceed assuming occupied rows 
    //On exit, if JOBZ == 'V' or 'v' of NULL A contains the eigenvectors
    //W - On exit, the eigenvalues of A in ascending order. If prev_eigenvalues, then W is passed to form the initial subintervals for spectrum slicing assuming ascending order.
    //In this case prev_eigenvalues should also be the number of eigenvalues.
    //JOBZ - 'V' or 'v' copy vectors, 'N' or 'n' do not copy vectors
    //UPLO - 'U' or 'u' - upper triangular, 'L' or 'l', lower triangular, NULL symmetric
    //prev_eigenvalues - see definition of W
    //interval - on input the interval over which the problem is solved interval[0] - lower bound, interval[1] upperbound.
    //prev_eigenvalues overrides this selection. On output - the interval that was actually used
    //interval_optimization - If used, the interval will be expanded and shrunk to fit the entire eigenspectrum within the interval. 
    //Error Codes 1 - Invalide JOBZ, 2 - Invalide UPLO, 
   
    PetscBool isInitialized;
    PetscInt ierr;
  	EPS eps;
	
    SlepcInitialized(&isInitialized);
    if(!isInitialized)
    SlepcInitialize(NULL,NULL,NULL,NULL); 	
    
    //printf("Double precision generalized eigenvalue problem dimensions NXLDAxLDB: %d %d %d\n Matrix A: %lx matrix B: %lx\n",*N,*LDA,*LDB,*Aprime,*Bprime);
   
    
    if(UPLO != NULL)
    {
        if(*UPLO == 'U' || *UPLO == 'u')
        { 
            UTtoSym(A,*N,*N);
            UTtoSym(B,*N,*N);
        }
        else if(*UPLO == 'L' || *UPLO == 'l')
        {
            LTtoSym(A,*N,*N);
            LTtoSym(B,*N,*N);
        }
        else
        return 2;
    }
    
    Mat A_Mat = squareToPetsc(A,*N);
    Mat B_Mat = squareToPetsc(B,*N);
	
    if(prev_eigenvalues)
        SIPSolve(&A_Mat, &B_Mat, interval, interval_optimization ,W,prev_eigenvalues,&eps);
    else
        SIPSolve(&A_Mat, &B_Mat, interval, interval_optimization ,NULL,NULL,&eps);
    
    ierr = EPSGetConverged(eps,N);CHKERRQ(ierr);
    
    if(JOBZ == NULL)
        squareFromEPS(&eps,A,W,*N);	
    else if(*JOBZ == 'V' || *JOBZ == 'v')
        squareFromEPS(&eps,A,W,*N);	
    else if(*JOBZ == 'N' || *JOBZ == 'N')
        eigenvaluesFromEPS(&eps,W,*N);
    else
        return 1;
        
    
    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
    ierr = MatDestroy(&A_Mat);CHKERRQ(ierr);
    ierr = MatDestroy(&B_Mat);CHKERRQ(ierr);
    	
	return 0;
}

