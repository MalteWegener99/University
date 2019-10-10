#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    int myRank, numProcs;
    int pingCount = 0;
    int pongCount = 0;

/* Initialize the infrastructure necessary for communication */
    MPI_Init(&argc, &argv);

    /* Identify this process */
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    /* Find out how many total processes are active */
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);


    if (myRank == 0)
    {
        printf("Sending Ping (# %i)\n", pingCount);
        MPI_Send(&pingCount, sizeof(int), MPI_INT, 1,
                    0, MPI_COMM_WORLD);

        MPI_Recv(&pongCount, sizeof(int), MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Received Pong (# %i)\n", pongCount);
    }
    // TODO: Do proper receive and send in any other process
    else
    {
        MPI_Recv(&pingCount, sizeof(int), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Received Ping (# %i)\n", pingCount);

        pongCount = -1*pingCount;
        MPI_Send(&pongCount, sizeof(int), MPI_INT, 0,
                    0, MPI_COMM_WORLD);
        printf("Sending Pong (# %i)\n", pongCount);

    }

    MPI_Finalize();

    return 0;
}
