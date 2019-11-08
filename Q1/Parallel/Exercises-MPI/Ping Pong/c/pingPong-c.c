#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Maximum array size 2^10 = 1024 elements
#define MAX_ARRAY_SIZE (1<<10)

int main(int argc, char **argv)
{
    // Variables for the process rank and number of processes
    int myRank, numProcs;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    // Find out MPI communicator size and process rank
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if (numProcs <= 1)
    {
        printf("Too few Processors");
        exit(1);
    }

    // PART B
    srandom(MPI_Wtime()*100000 + myRank*137);
    int numberOfElementsToSend = random() % 100;
    // Allocate an array big enough to hold event the largest message
    int *myArray = (int *)malloc(sizeof(int)*MAX_ARRAY_SIZE);
    if (myArray == NULL)
    {
        printf("Not enough memory\n");
        exit(1);
    }
    int numberOfElementsReceived;

    // Have only the first process execute the following code
    if (myRank == 0)
    {
        printf("We have %i Processseses and i am %i\n", numProcs, myRank);
        for(int i = 1; i <= (numProcs-1); i++){
            printf("Sending %d elements to %d \n", numberOfElementsToSend, i);
            MPI_Send(&numberOfElementsToSend, 1, MPI_INT, i,
                        0, MPI_COMM_WORLD);
            MPI_Send(myArray, numberOfElementsToSend, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    else // myRank == 1
    {
        MPI_Recv(&numberOfElementsReceived, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Recv(myArray, numberOfElementsReceived, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Received %i elements on %i\n", numberOfElementsReceived, myRank);

        printf("Sending back %i elements from %i \n", numberOfElementsToSend, myRank);
        MPI_Send(&numberOfElementsToSend, 1, MPI_INT, 0,
                    0, MPI_COMM_WORLD);
        MPI_Send(myArray, numberOfElementsToSend, MPI_INT, 0, 0, MPI_COMM_WORLD);
        printf("DONE on Process %i \n", myRank);
    }
    if (myRank == 0){
                for(int i = 1; i <= (numProcs-1); i++){
        

            MPI_Recv(&numberOfElementsReceived, 1, MPI_INT, i,
                        0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Recv(myArray, numberOfElementsReceived, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            printf("Received %d elements from %i \n", numberOfElementsReceived, i);
        }
        printf("DONE on Process %i \n", myRank);
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
