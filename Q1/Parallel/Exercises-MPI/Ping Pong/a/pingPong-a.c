#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    int myRank, numProcs;
    int number = 1234;

/* Initialize the infrastructure necessary for communication */
    MPI_Init(&argc, &argv);

    /* Identify this process */
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    /* Find out how many total processes are active */
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);


    if (myRank == 0)
    {
        MPI_Send(&number, sizeof(int), MPI_INT, 1,
                    0, MPI_COMM_WORLD);

        MPI_Recv(&number, sizeof(int), MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Received %i on %d\n", number, myRank);
    }
    // TODO: Do proper receive and send in any other process
    else
    {
        MPI_Recv(&number, sizeof(int), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Received Ping %i on %d\n", number, myRank);

        number = -1*number;
        MPI_Send(&number, sizeof(int), MPI_INT, 0,
                    0, MPI_COMM_WORLD);

    }

    MPI_Finalize();

    return 0;
}
