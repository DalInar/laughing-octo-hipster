#include <iostream>
#include <math.h>
#include <mpi.h>

enum partition {squares,rows,cols};

void print_subgrid(double * A, int rank, int size_x, int size_y, partition part) {
	if (part == cols) {
		std::cout<<"(Rank: "<<rank<<", Col: GL) ";
		for(int row=0; row < size_y; row++) {
			std::cout<<A[row]<<" ";
		}
		std::cout<<std::endl;

		for(int col=0; col<size_x; col++) {
			std::cout<<"(Rank: "<<rank<<", Col: "<<col<<") ";
			for(int row=0; row < size_y; row++) {
				std::cout<<A[size_y+col*size_y+row]<<" ";
			}
			std::cout<<std::endl;
		}

		std::cout<<"(Rank: "<<rank<<", Col: GR) ";
		for(int row=0; row < size_y; row++) {
			std::cout<<A[row+(size_x+1)*size_y]<<" ";
		}
		std::cout<<std::endl<<std::endl;
	}
	else {
		std::cout<<"Error! Invalid partition in print_subgrid function."<<std::endl;
	}
}

void print_grid(double * A, int rank, int comm_size, int size_x, int size_y, partition part) {
	int flag;
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		flag=1;
		print_subgrid(A,rank,size_x,size_y,part);
		if(comm_size>1) {
			MPI_Send(&flag,1,MPI_INT,1,0,MPI_COMM_WORLD);
		}
	}
	else {
		MPI_Recv(&flag,1,MPI_INT,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		print_subgrid(A,rank,size_x,size_y,part);
		if(comm_size > rank+1) {
			MPI_Send(&flag,1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
}

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
	int col_offset = rank*size_x;	//Each processor needs to know which columns it actually has from A

	if(rank == size-1) {
		size_x = n - (size-1)*size_x;
	}

	std::cout<<"Rank: "<<rank<<", size_x: "<<size_x<<", size_y: "<<size_y<<std::endl;

	double * A;
	double * B;
	double * C;
	int length = size_y*(size_x+2);		//Space for this procs grid, plus ghost (halo) cells
	A = new double [length];	//Need two grid so that updates do not overwrite
	B = new double [length];	//data that is still needed

	int row,col;
	for(int i=size_y; i<length-size_y; i++) {
		row = i%size_y;
		col = (i-i%size_y)/size_y + col_offset;
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
	print_grid(A, rank, size, size_x, size_y, cols);

	int sendL=(rank-1+size)%size;
	int recvL=sendL;
	int sendR=(rank+1)%size;
	int recvR=sendR;

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0) {
		start_t=MPI_Wtime();
	}

	//nonblocking send first col
	MPI_Send(&A[size_y],size_y,MPI_DOUBLE,sendL,0,MPI_COMM_WORLD);
	//nonblocking send last col
	MPI_Send(&A[size_x*size_y],size_y,MPI_DOUBLE,sendR,0,MPI_COMM_WORLD);
	//nonblocking rec first ghost col
	MPI_Request requestL;
	MPI_Irecv(&A[0],size_y,MPI_DOUBLE,recvL,0,MPI_COMM_WORLD,&requestL);
	//nonblocking rec last ghost col
	MPI_Request requestR;
	MPI_Irecv(&A[(size_x+1)*size_y],size_y,MPI_DOUBLE,recvR,0,MPI_COMM_WORLD,&requestR);

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
