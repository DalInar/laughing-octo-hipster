#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[]) {
	int rank, size;
	MPI_Init(&argc, &argv); //Start up MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get current processor number
	MPI_Comm_size(MPI_COMM_WORLD, &size); //Get number of processors in communicator
	std::cout<<"Hello worlds from process "<<rank<<" of "<<size<<std::endl;
	MPI_Finalize();

  return 0;
}
