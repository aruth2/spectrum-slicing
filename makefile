#This is the makefile for spectrum-slicing interface. 

#suppresses a warning from CHKERRQ
CFLAGS     = -Wno-int-conversion -Wno-implicit-function-declaration

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

DEPS = sips.h sips_square.h
#USE_SCALAPACK must also be defined inside sips_squaretest.c
USE_SCALAPACK = 0

clean::
	${RM} libsips.so libsips_square.so libscalapack_interface.so sips sips_square

#Library for sips calls
sips: sips.o chkopts
	-${CLINKER} -shared -o libsips.so sips.o ${SLEPC_LIB}
	${RM} sips.o

#Executable, accepts two petsc matrices as command like arguments and diagonalizes using -fA A_matrix_file -fB B_matrix_file
sipstest: sips.o sipstest.o chkopts
	-${CLINKER} -o sips sips.o sipstest.o ${SLEPC_LIB}
	${RM} sips.o sipstest.o

#A dynamically linked executable for sipstest.
sipstest_shared: sips sipstest.o chkopts
	-${CLINKER} -o sips sipstest.o ${SLEPC_LIB} -L. -lsips
	${RM} sipstest.o

#Library for calls between square double matrices and PETSc matrices.
sips_square: sips sips_square.o chkopts
	-${CLINKER} -shared -o libsips_square.so sips_square.o ${SLEPC_LIB} -L. -lsips
	${RM} sips_square.o	

#Library for calls to scalapack using square double matrix on every process
scalapack_interface: scalapack_interface.o chkopts
	-${CLINKER} -shared -o libscalapack_interface.so scalapack_interface.o ${SLEPC_LIB}
	${RM} scalapack_interface.o	


#Test program for initializing square double matrix(s) and diagonalizing. 
ifeq ($(USE_SCALAPACK),1)
sips_squaretest: sips.o sips_square.o sips_squaretest.o scalapack_interface.o chkopts
	-${CLINKER} -o sips_square sips.o sips_square.o sips_squaretest.o scalapack_interface.o ${SLEPC_LIB} 
	${RM} sips.o sips_square.o sips_squaretest.o scalapack_interface.o
else
sips_squaretest: sips.o sips_square.o sips_squaretest.o chkopts
	-${CLINKER} -o sips_square sips.o sips_square.o sips_squaretest.o ${SLEPC_LIB}
	${RM} sips.o sips_square.o sips_squaretest.o
endif	


#Dynamically linked executable for sips_squaretest
ifeq ($(USE_SCALAPACK),1)
sips_squaretest_shared: sips sips_square sips_squaretest.o chkopts scalapack_interface
	-${CLINKER} -o sips_square sips_squaretest.o ${SLEPC_LIB} -L. -lsips -lsips_square -lscalapack_interface 
else
sips_squaretest_shared: sips sips_square sips_squaretest.o chkopts
	-${CLINKER} -o sips_square sips_squaretest.o ${SLEPC_LIB} -L. -lsips -lsips_square
endif	
	${RM} sips_squaretest.o 

.DEFAULT_GOAL := sips_squaretest_shared
 
#------------------------------------------------------------------------------------

.PHONY: test

test:: runsips_square_mpi1
	

runsips_square_verbose: 
	${MPIEXEC} -n 2 ./sips_square -rows 5 

runsips_square_v:
	valgrind ./sips_square -log_view

runsips_square_l:
	time ./sips_square -lapack -rows 4000

runsips_square: 
	time ./sips_square -rows 10000

runsips_square_mpi1: 
	time ${MPIEXEC} -n 1 ./sips_square -rows 4000
	
runsips_square_mpi1_4columns: 
	time ${MPIEXEC} -n 1 ./sips_square -rows 4000 -nonzerodiagonals 4	
	
runsips_square_mpi2: 
	time ${MPIEXEC} -n 2 ./sips_square -rows 4000

NUMBERS = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
runsips_square_time: 
		$(foreach var,$(NUMBERS),time ${MPIEXEC} -n $(var) ./sips_square -rows 4000 -nonzerodiagonals 100;)

runsips_square_time_double: 
		$(foreach var,$(NUMBERS),time ${MPIEXEC} -n $(var) ./sips_square -rows 4000 -nonzerodiagonals 100 -doublediag;)

runsscalapack_square_time: 
		$(foreach var,$(NUMBERS),time ${MPIEXEC} -n $(var) ./sips_square -rows 4000 -nonzerodiagonals 100 -scalapack;)
