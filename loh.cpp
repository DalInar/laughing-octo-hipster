#include <iostream>
#include <math.h>
#include <mpi.h>

enum partition {squares,rows,cols};

double square_partition(int n, int rank, int size, int num_iter) {
	double start_t, end_t;
	int test=0;

	if(rank==0){
		std::cout<<"Partitioning with squares"<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0) {
		start_t=MPI_Wtime();
	}

	for(int i=0; i<1000000000; i++) {
		test+=1;
	}




	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0) {
		end_t=MPI_Wtime();
	}

	return end_t-start_t;
}

double col_partition(int n, int rank, int size, int num_iter) {
	if(rank==0){
			std::cout<<"Partitioning with columns"<<std::endl;
	}

	double start_t, end_t;
	int size_x = ceil((1.0*n)/size);;
	int size_y = n;

	if(rank == size-1) {
		size_x = n - (size-1)*size_x;
	}

	std::cout<<"Rank: "<<rank<<", size_x: "<<size_x<<", size_y: "<<size_y<<std::endl;

	double * A,B;
	int length = size_y*(size_x+2);		//Space for this procs grid, plus ghost (halo) cells
	A = new double [length];	//Need two grid so that updates do not overwrite
	B = new double [length];	//data that is still needed

	int row,col;
	for(int i=0; i<length; i++) {
		row = i%size_y;
		col = (i-i%size_y)/size_y;
		if(row == 0) {
			A[i] = 0.0;
		}
		else if(row == size_y-1) {
			A[i] = 5*sin(M_PI * ((1.0*col)/n) * ((1.0*col)/n));
		}
		else {
			A[i] = 0.5;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0) {
		start_t=MPI_Wtime();
	}

	//nonblocking send first col
	//nonblocking send last col
	//nonblocking rec first ghost col
	//nonblocking rec last ghost col

	//update B (inner cols first, then check to see if ghost data has arrived to update end cols)
	//swap A,B pointers

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0) {
		end_t=MPI_Wtime();
	}

	return end_t-start_t;
}

int main(int argc, char *argv[]) {
	int n=10;	//Grid is nxn
	int num_iter=500;
	partition part = cols;
	double time;

	int rank, size;
	MPI_Init(&argc, &argv); //Start up MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get current processor number
	MPI_Comm_size(MPI_COMM_WORLD, &size); //Get number of processors in communicator
	std::cout<<"Process "<<rank<<" of "<<size<<" online"<<std::endl;

	MPI_Barrier(MPI_COMM_WORLD);
	if(part==cols) {
		time = col_partition(n,rank,size,num_iter);
	}
	else if(part==squares) {
		time = square_partition(n,rank,size,num_iter);
	}
	else {
		std::cout<<"Error! Invalid partition!"<<std::endl;
	}

	if(rank==0) {
		std::cout<<"Elapsed time: "<<time<<" seconds"<<std::endl;
	}

	MPI_Finalize();

  return 0;
}
