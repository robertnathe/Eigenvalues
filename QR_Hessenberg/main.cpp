include <stdio.h>
#include <mpi.h>
#include <iostream>
#include <chrono>
int main(int argc, char** argv) {

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

	// Using time point and system_clock
    std::chrono::time_point<std::chrono::system_clock> start, end;
 
    start = std::chrono::system_clock::now();

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    end = std::chrono::system_clock::now();
 
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
 
//	std::cout << "finished computation at " << std::ctime(&end_time)
//              << "elapsed time: " << elapsed_seconds.count() << "s\n";
              
    // Finalize the MPI environment.  
    MPI_Finalize();
}
