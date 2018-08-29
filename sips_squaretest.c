#include "sips_square.h"

//#define USE_SCALAPACK

#ifdef USE_SCALAPACK
#include "scalapack_interface.h"
#endif

int solver = 0; //0 SIPS, 1 LAPACK, 2 SCALAPACK
int rows = 10000;
int nonzerodiagonals = 2;
int doublediag = 0;

int doubletest(int argc, char **argv)
{
	printf("Attempting to diagonalize a double matrix\n");	
	int i,j;
	int cols = rows;
	double *Hamiltonian = malloc(rows*rows*sizeof(double));
	//double value;
    int nconv=rows;
    for(i=0;i<rows;i++)
	{
		//nearest neighbor hamiltonian of 1-D chain
		 for(j=0;j<nonzerodiagonals/2;j++)
        {
		if(i-1-j>=0)
		Hamiltonian[i*cols+i-1-j] = 1;
		if(i+1+j<rows)
		Hamiltonian[i*cols+i+1+j] = 1;
        }
        
		//Identity Matrix
		Hamiltonian[i*cols+i] = 0.001;
	}
	char jobz = 'v';
	char uplo = 'u';
	int lwork = rows * rows * 4;
	double *work = malloc(lwork*sizeof(double));
	int liwork = 6 * rows + 3;
	double *iwork = malloc(liwork*sizeof(double)); 
	int info = 0;//After dsyevd info = 0 if the calculation went correctly	
	double *eigenvalues = malloc(rows*sizeof(double));

	if(verbose)
	{
	printf("Hamiltonian before:\n");
	printmatrix(Hamiltonian,rows,cols,1,1,1);

	}
	
	switch(solver)
    {
    case 0:
    if(doublediag)
    {
        //First diagonalization uses interval optimization and does not copy eigenvectors so matrix can be reused, second diagonalization uses previous eigenvalues to set subinterval
        const char *JOBZ = "N";
        dsygvsx( &nconv, Hamiltonian, NULL, eigenvalues, 1, NULL, NULL, (char *)JOBZ, NULL);
        dsygvsx( &nconv, Hamiltonian, NULL, eigenvalues, NULL, NULL, nconv, NULL, NULL);
    }
    else
    dsyevs( &rows, Hamiltonian, eigenvalues);
    break;
    case 1:
	//Diagonalize with Lapack
	dsyev_( &jobz, &uplo, &rows, Hamiltonian, &cols, eigenvalues, work, &lwork, iwork,  &liwork, &info );
	break;
    case 2:
    //Solve with SCALAPACK
#ifdef USE_SCALAPACK    
    solveEigenvalueScalapack(Hamiltonian, eigenvalues, rows);
#endif
    break;
    
    default:
    printf("No solver found\n");
    break;
    }
	if(verbose)
	{
	printf("Eigenvectors:\n");
	printmatrix(Hamiltonian,nconv,cols,1,1,1);
	
	printf("Eigenvalues:\n");
	printmatrix(eigenvalues,nconv,1,1,1,0);
	}
	//memcpy(eigenvectors,Hamiltonian,rows*cols*sizeof(double));
	//double *leftEigenVectors = matrix;
	//free(work);
	//free(iwork);
	return 0;
}

int generalizeddoubletest(int argc, char **argv)
{
	printf("Attempting to diagonalize a generalized double matrix\n");	
	int i,j;
	int cols = rows;
	double *Hamiltonian = malloc(rows*rows*sizeof(double));
	double *B = malloc(rows*rows*sizeof(double));
	for(i=0;i<rows;i++)
	{
		//nearest neighbor hamiltonian of 1-D chain
		 for(j=0;j<nonzerodiagonals/2;j++)
        {
		if(i-1-j>=0)
		Hamiltonian[i*cols+i-1-j] = 1;
		if(i+1+j<rows)
		Hamiltonian[i*cols+i+1+j] = 1;
        }
        
		//Identity Matrix
		B[i*cols+i] = 1;
		//Hamiltonian[i*cols+i] = 0.001;
	}
	char jobz = 'v';
	char uplo = 'u';
	int itype = 1;
	int lwork = rows * rows * 4;
	double *work = malloc(lwork*sizeof(double));
	int liwork = 6 * rows + 3;
	double *iwork = malloc(liwork*sizeof(double)); 
	int info = 0;//After dsyevd info = 0 if the calculation went correctly	
	double *eigenvalues = malloc(rows*sizeof(double));

	if(verbose)
	{
	printf("Hamiltonian before:\n");
	printmatrix(Hamiltonian,rows,cols,1,1,1);

	printf("B before:\n");
	printmatrix(B,rows,cols,1,1,1);
	}
	
	switch(solver)
    {
	case 0:
    dsygvs( &rows, Hamiltonian, B , eigenvalues);
    break;
    
    case 1:
    //Diagonalize with Lapack
	dsygv_( &itype, &jobz, &uplo, &rows, Hamiltonian, &cols, B, &cols, eigenvalues, work, &lwork, iwork,  &liwork, &info );
	break;
    
        case 2:
    //Solve with SCALAPACK
#ifdef USE_SCALAPACK    
#endif
    break;
    
    default:
    printf("No solver found\n");
    break;
    
	}
	
	if(verbose)
	{
	printf("Eigenvectors:\n");
	printmatrix(Hamiltonian,rows,cols,1,1,1);
	
	printf("Eigenvalues:\n");
	printmatrix(eigenvalues,rows,1,1,1,0);
	}

	return 0;
}





int main(int argc,char **argv)
{
	printdatatypes();
	int i;
	for(i=0;i<argc;i++)
	{
        if(!strcmp(argv[i],"-lapack"))
            solver = 1;
        if(!strcmp(argv[i],"-scalapack"))
            solver = 2;
        //	if(!strcmp(argv[i],"-noverbose"))
        //	verbose = 0;
        if(!strcmp(argv[i],"-rows"))
            rows = atoi(argv[i+1]);
        if(!strcmp(argv[i],"-nonzerodiagonals"))
            nonzerodiagonals = atoi(argv[i+1]);
        if(!strcmp(argv[i],"-doublediag"))
            doublediag = 1;
	}
	
    SlepcInitialize(NULL,NULL,NULL,NULL);
	doubletest(argc,argv);
	//generalizeddoubletest(argc,argv);
	
    SlepcFinalize();
	return 0;
}
