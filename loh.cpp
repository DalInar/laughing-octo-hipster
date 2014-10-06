#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
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

double col_partition(int n, int rank, int size, int num_iter, std::ofstream & output) {
	if(rank==0){
			std::cout<<"Partitioning with columns on "<<n<<" by "<<n<<std::endl;
	}

	double start_t, end_t, local_verify, global_verify;
	int size_x = ceil((1.0*n)/size);
	int size_y = n;
	int col_offset = rank*size_x;	//Each processor needs to know which columns it actually has from A

	if(rank == size-1) {
		size_x = n - (size-1)*size_x;
	}

	double * A;
	double * B;
	double * C;
	int length = size_y*(size_x+2);		//Space for this procs grid, plus ghost (halo) cells
	A = new double [length];	//Need two grid so that updates do not overwrite
	B = new double [length];	//data that is still needed

	int row,col;
	for(int i=size_y; i<length-size_y; i++) {
		row = i%size_y;
		col = (i-i%size_y)/size_y - 1 + col_offset;
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
	//print_grid(A, rank, size, size_x, size_y, cols);

	int sendL=(rank-1+size)%size;
	int recvL=sendL;
	int sendR=(rank+1)%size;
	int recvR=sendR;

	std::cout<<"Rank: "<<rank<<", size_x: "<<size_x<<", size_y: "<<size_y<<std::endl;

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0) {
		start_t=MPI_Wtime();
	}

	for(int k=0; k<500; k++) {
		//nonblocking send first col
		MPI_Request sendrequestL;
		MPI_Isend(&A[size_y],size_y,MPI_DOUBLE,sendL,0,MPI_COMM_WORLD, &sendrequestL);

		//nonblocking send last col
		MPI_Request sendrequestR;
		MPI_Isend(&A[size_x*size_y],size_y,MPI_DOUBLE,sendR,1,MPI_COMM_WORLD, &sendrequestR);	//Use tags in case two procs communicate on both left and right
		//nonblocking rec first ghost col

		MPI_Request requestL;
		MPI_Irecv(&A[0],size_y,MPI_DOUBLE,recvL,1,MPI_COMM_WORLD,&requestL);
		//nonblocking rec last ghost col
		MPI_Request requestR;
		MPI_Irecv(&A[(size_x+1)*size_y],size_y,MPI_DOUBLE,recvR,0,MPI_COMM_WORLD,&requestR);

		//update B (inner cols first, then check to see if ghost data has arrived to update end cols)
//		for(int i=2*size_y; i<length-2*size_y; i++) {
//			row = i%size_y;
//			if(row!=0 && row!=n-1) {
//				col = (i-i%size_y)/size_y - 1;
//				B[i]=(A[i]+
//						A[size_y+row+(col-1)*size_y]+A[size_y+row+(col+1)*size_y]+
//						A[size_y+row-1+(col)*size_y]+A[size_y+row+1+(col)*size_y]+
//						A[size_y+row-1+(col-1)*size_y]+A[size_y+row-1+(col+1)*size_y]+
//						A[size_y+row+1+(col-1)*size_y]+A[size_y+row+1+(col+1)*size_y])/9;
//			}
//			else {
//				B[i]=A[i];
//			}
//		}

		//different indexing scheme
		int i;
		for(col=1; col<size_x-1; col++) {
			i = size_y + col*size_y;
			B[i]=A[i];
			B[i+size_y-1]=A[i+size_y-1];
			for(row=1; row<size_y-1; row++) {
				i+=1;
				B[i]=A[i-1]+A[i]+A[i+1];
				B[i]+=A[i-size_y-1]+A[i-size_y]+A[i-size_y+1];
				B[i]+=A[i+size_y-1]+A[i+size_y]+A[i+size_y+1];
				B[i]=B[i]/9;
			}
		}

		//Wait until ghost cell data arrives
		MPI_Wait(&requestL,MPI_STATUS_IGNORE);
		MPI_Wait(&requestR,MPI_STATUS_IGNORE);

		//Update boundary columns
//		B[size_y]=A[size_y];
//		B[2*size_y-1]=A[2*size_y-1];
//		B[(size_x)*size_y]=A[(size_x)*size_y];
//		B[(size_x+1)*size_y-1]=A[(size_x+1)*size_y-1];
//		for(int i=1; i<n-1; i++) {
//			row=i;
//			col = 0;
//			B[size_y+i]=(A[size_y+i]+
//					A[size_y+row+(col-1)*size_y]+A[size_y+row+(col+1)*size_y]+
//					A[size_y+row-1+(col)*size_y]+A[size_y+row+1+(col)*size_y]+
//					A[size_y+row-1+(col-1)*size_y]+A[size_y+row-1+(col+1)*size_y]+
//					A[size_y+row+1+(col-1)*size_y]+A[size_y+row+1+(col+1)*size_y])/9;
//			col = size_x-1;
//			B[(size_x)*size_y+i]=(A[(size_x)*size_y+i]+
//					A[size_y+row+(col-1)*size_y]+A[size_y+row+(col+1)*size_y]+
//					A[size_y+row-1+(col)*size_y]+A[size_y+row+1+(col)*size_y]+
//					A[size_y+row-1+(col-1)*size_y]+A[size_y+row-1+(col+1)*size_y]+
//					A[size_y+row+1+(col-1)*size_y]+A[size_y+row+1+(col+1)*size_y])/9;
//		}

		for(col=0; col<size_x; col=col+size_x-1) {
			i = size_y + col*size_y;
			B[i]=A[i];
			B[i+size_y-1]=A[i+size_y-1];
			for(row=1; row<size_y-1; row++) {
				i+=1;
				B[i]=A[i-1]+A[i]+A[i+1];
				B[i]+=A[i-size_y-1]+A[i-size_y]+A[i-size_y+1];
				B[i]+=A[i+size_y-1]+A[i+size_y]+A[i+size_y+1];
				B[i]=B[i]/9;
			}
		}

		//swap A,B pointers
		C=B;
		B=A;
		A=C;
	}

	local_verify=0.0;
	for(col=0; col<size_x; col++) {
		row=col+col_offset;
		local_verify+=A[size_y+col*size_y+row];
	}

	MPI_Reduce(&local_verify, &global_verify, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0) {
		end_t=MPI_Wtime();
		std::cout<<"Global Verification = "<<global_verify<<std::endl;
	}

	if(rank==0) {
		output<<"Time\t"<<end_t-start_t<<std::endl;
		output<<"Verify\t"<<global_verify<<std::endl;
	}

	delete []A;
	delete []B;
	//print_grid(A, rank, size, size_x, size_y, cols);

	return end_t-start_t;
}

int main(int argc, char *argv[]) {
	int n=atoi(argv[1]);	//Grid is nxn
	int num_iter=500;
	partition part = cols;
	double time,verify;

	int rank, size;
	MPI_Init(&argc, &argv); //Start up MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get current processor number
	MPI_Comm_size(MPI_COMM_WORLD, &size); //Get number of processors in communicator
	std::cout<<"Process "<<rank<<" of "<<size<<" online"<<std::endl;

	std::ofstream output;
	std::stringstream filename;
	if(rank==0) {
		filename << "loh_n"<<n<<"_p"<<size;
		output.open(filename.str().c_str());
		output<<"Size\t"<<n<<std::endl;
		output<<"Procs\t"<<size<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(part==cols) {
		time = col_partition(n,rank,size,num_iter,output);
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

	output.close();

	MPI_Finalize();

  return 0;
}
