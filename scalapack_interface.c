/*
 * scalapack_interface.c
 */

#include "scalapack_interface.h"
#define chkerr if(info != 0)printf("Error code is %d\n",info);
enum {DTYPE_MAT,CTXT_MAT,M_MAT,N_MAT,MB_MAT,NB_MAT,RSRC_MAT,CSRC_MAT,LLD_MAT};
int one = 1, zero=0,minusone=-1;


#define DEBUG
#ifndef DEBUG
#define printf //printf 
#endif
int localToGlobal1DZeroStart(int NB, int NP, int p, int localindex)  
{
	/* Based on http://netlib.org/scalapack/slug/node76.html
	 * NB (NB) - the block size
	 * NP (P) - the number of processes the array is distributed over
	 * p (p) - the id of the process the globalindex belongs to
	 * localindex (l*NB+x) - the index that the data should be located at in the local array 
	 * returns globalindex (I) - the index of the full-size global array
	 * 
	 * Fortran versions need to be corrected for arrays starting at 1 
	 * globalindex, localindex, and x should all be 1 greater, whilst p and l remain unchanged.
	 * 
	 * int x = localindex % NB;
	 * int l = (localindex-1) / NB;
	 * int globalindex = (l*NP+p)*NB+x;
	 * */
	 
	 
	 int x = localindex % NB;// (x) - the index within the block 
	 int l = localindex / NB;// (l) - the block number of the data
	 
	 int globalindex = (l*NP+p)*NB+x;
	 return globalindex;
}

int localToGlobal1DOneStart(int NB, int NP, int p, int localindex)  
{
	/* Based on http://netlib.org/scalapack/slug/node76.html
	 * NB (NB) - the block size
	 * NP (P) - the number of processes the array is distributed over
	 * p (p) - the id of the process the globalindex belongs to
	 * localindex (l*NB+x) - the index that the data should be located at in the local array 
	 * returns globalindex (I) - the index of the full-size global array
	 * 
	 * */
	 
	 int x = (localindex-1) % NB+1;
	 int l = (localindex-1) / NB;
	 int globalindex = (l*NP+p)*NB+x;
	 
	 return globalindex;
}

void globalToLocal1DZeroStart(int globalindex, int NB, int NP, int *p, int *localindex)  
{
	/* Based on http://netlib.org/scalapack/slug/node76.html
	 * globalindex - the index of the full-size global array
	 * NB - the block size
	 * NP - the number of processes the array is distributed over
	 * returns p - the id of the process the globalindex belongs to
	 * returns localindex - the index that the data should be located at in the local array 
	 * 
	 * Fortran versions need to be corrected for arrays starting at 1 
	 * globalindex, localindex, and x should all be 1 greater, whilst p and l remain unchanged.
	 * 
	 * int x = ((globalindex-1) % NB) + 1;
	 * int l = (globalindex-1) / (NP*NB);
	 * *p = ((globalindex-1)/(NB)) % NP;
	 * *localindex = l*NB+x;
	 * */
	 
	 int x = globalindex % NB;// (x) - the index within the block //Add 1 for fortran
	 int l = globalindex / (NP*NB);// (l) - the block number of the data //Add 1 for fortran
	 *p = (globalindex/(NB)) % NP;
	 *localindex = l*NB+x;
}

void globalToLocal1DOneStart(int globalindex, int NB, int NP, int *p, int *localindex)  
{
	/* Based on http://netlib.org/scalapack/slug/node76.html
	 * globalindex - the index of the full-size global array
	 * NB - the block size
	 * NP - the number of processes the array is distributed over
	 * returns p - the id of the process the globalindex belongs to
	 * returns localindex - the index that the data should be located at in the local array 
	 * */
	 
	 int x = ((globalindex-1) % NB) + 1;
	 int l = (globalindex-1) / (NP*NB);
	 *p = ((globalindex-1)/(NB)) % NP;
	 *localindex = l*NB+x;
}

void transferMatrixToScalapack(double *A, double *Alocal, int M, int MB, int N, int NB)
{
	int rank, size, NPROW, NPCOL, MYROW, MYCOL, mlocal, nlocal;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
	printf("MPI rank,size %d %d\n",rank,size);
    NPROW = sqrt(size); NPCOL = size/NPROW;     
	MYROW = rank/NPCOL; MYCOL = rank % NPCOL;
	mlocal    = numroc_( &M, &MB, &MYROW, &zero, &NPROW );
    nlocal    = numroc_( &N, &NB, &MYCOL, &zero, &NPCOL );
	
	int lrow,lcolumn,grow,gcolumn;
	
	for(lrow=0;lrow<mlocal;lrow++)
	{
		grow = localToGlobal1DZeroStart(MB,NPROW,MYROW,lrow);
		//printf("Local row %d maps to global row %d\n",lrow,grow);
		for(lcolumn=0;lcolumn<nlocal;lcolumn++)
		{
			gcolumn = localToGlobal1DZeroStart(NB,NPCOL,MYCOL,lcolumn);
			//printf("Local column %d maps to global column %d\n",lcolumn,gcolumn);
			Alocal[lcolumn*mlocal+lrow] = A[gcolumn*M+grow];		
		}
	}
}

void transferScalapackToMatrix(double *A, double *Alocal, int M, int MB, int N, int NB)
{
	// Here we copy the matrix over to every single process. We do this by telling Scalapack that there is a (M*NPROW x N*NPCOL) matrix
	// and each (MxN) submatrix is in fact the matrix we want. Because of Scalapacks block distribution, this will result in a correct local (sub)matrix
	int rank, size, NPROW, NPCOL, MYROW, MYCOL, mlocal;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
	printf("MPI rank,size %d %d\n",rank,size);
    NPROW = sqrt(size); NPCOL = size/NPROW;     
	MYROW = rank/NPCOL; MYCOL = rank % NPCOL;
	mlocal    = numroc_( &M, &MB, &MYROW, &zero, &NPROW );
//    nlocal    = numroc_( &N, &NB, &MYCOL, &zero, &NPCOL );

	printf("Global process grid %d x %d %d x %d\n",NPROW,NPCOL,MYROW,MYCOL);
  	int descgather[9], desca[9];
    int bigM = M*NPROW,bigN = N*NPCOL;
    int info;
	int ictxt;
	sl_init_( &ictxt, &NPROW, &NPCOL);
	descinit_(descgather, &bigM, &bigN, &M, &N, &zero, &zero, &ictxt, &M, &info); chkerr
	descinit_(desca, &M, &N, &MB, &NB, &zero, &zero, &ictxt, &mlocal, &info); chkerr
	int prow,pcol;
	int coloffset,rowoffset;
	for(prow=0;prow<NPROW;prow++)
	{
		rowoffset = prow*M+1;
        for(pcol=0;pcol<NPCOL;pcol++)
        {
			coloffset = pcol*N+1;
			pdgemr2d_(&M, &N, Alocal, &one, &one, desca, A, &rowoffset, &coloffset, descgather, &ictxt);
		}
	}
}

void solveEigenvalueScalapack(double *A, double *lambda, int M)
{
    printf("Start of solve eigenvalue function\n");
    int MB = 64;
	int NB=MB;
	//Definition of the matrices
	 
	printf("Matrix A:\n");printmatrix(A,M,M,1,0,1);

	//Variables associated with SCALAPACK and MPI
    int rank,size;
	int NPROW, NPCOL;
	int MYROW, MYCOL;
	int desca[9],descev[9];
	int ictxt;
	int mlocal,nlocal;
	double *Alocal,*evlocal;
	//Variables used for LAPACK
	int info;
	//int ipiv[M];
	const char *JOBZ = "V";
	const char *UPLO = "U";
	double doubleLWORK, doubleLRWORK, doubleLIWORK;//The double versions are used for the workspace queries in scalapack     
	int LWORK, LRWORK, LIWORK;
	 
	//Initialize MPI     
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
	printf("MPI rank: %d, size: %d\n",rank,size);
    NPROW = sqrt(size); NPCOL = size/NPROW; 
     
    //Initials the scalapack process grid
	sl_init_( &ictxt, &NPROW, &NPCOL);
    blacs_gridinfo_( &ictxt, &NPROW, &NPCOL, &MYROW, &MYCOL);
    mlocal    = numroc_( &M, &MB, &MYROW, &zero, &NPROW );
    nlocal    = numroc_( &M, &NB, &MYCOL, &zero, &NPCOL );
         
    //Create descriptors for the matrices and initalize the local storage of the matrices
    printf("NumProcesses: %d x %d, MyProcess: %d rank %d x %d, Global matrix dimensions: %d x %d, local matrix dimensions: %d x %d\n",
    NPROW,NPCOL,rank, MYROW,MYCOL,M,M,mlocal,nlocal);
	descinit_( desca, &M, &M, &MB, &NB, &zero, &zero, &ictxt, &mlocal, &info); chkerr
	descinit_( descev, &M, &M, &MB, &NB, &zero, &zero, &ictxt, &mlocal, &info); chkerr
	Alocal = malloc(mlocal*nlocal*sizeof(double));
	evlocal = malloc(mlocal*nlocal*sizeof(double));
	
	transferMatrixToScalapack(A, Alocal, M, MB, M, MB);
	
	printf("Local Matrix A:\n");printmatrix(Alocal,mlocal,nlocal,1,0,1);
	  	
	pdsyevd_(JOBZ, UPLO, &M, Alocal, &one, &one, desca, lambda, evlocal, &one, &one, descev,
	 &doubleLWORK, &minusone, &doubleLRWORK, &minusone, &doubleLIWORK,&minusone,&info);
	
	LWORK = doubleLWORK; LRWORK=doubleLRWORK; LIWORK=doubleLIWORK;
	LRWORK = fmax(LRWORK,2*LWORK);//I keep getting LRWORK=0 from the workspace query
	LIWORK = 7*M+ 8*NPCOL+2;//Straight from manual
	printf("LWORK %d LRWORK %d LIWORK %d\n",LWORK,LRWORK,LIWORK);
	double *WORK = malloc(LWORK*sizeof(double));
	double *RWORK = malloc(LRWORK*sizeof(double));
	double *IWORK = malloc(LIWORK*sizeof(double)); 
	
	pdsyevd_(JOBZ, UPLO, &M, Alocal, &one, &one, desca, lambda, evlocal, &one, &one, descev,
	 WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &info); 
	
	free(WORK); free(RWORK); free(IWORK);
	if(info != 0)
	printf("Error code %d \n",info);
	 
	printf("Local eigenvectors\n"); printmatrix(evlocal,mlocal,nlocal,1,0,1);
	
	transferScalapackToMatrix(A, evlocal, M, MB, M, MB); 
	
    //Transfer the eigenvalues from the rank 0 process to all processes
    MPI_Bcast(lambda, M, MPI_DOUBLE, 0,MPI_COMM_WORLD);
     
	printf("Eigenvectors of A:\n");printmatrix(A,M,M,1,0,1);
	printf("Eigenvalues lambda_i:\n");printmatrix(lambda,M,1,1,0,1);
	 
	free(Alocal); free(evlocal);
} 

