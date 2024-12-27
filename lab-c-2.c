#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    int rank, size;
    int data;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if ((size & (size - 1)) != 0) {
        if (rank == 0) {
            printf("Number of processors must be a power of 2\n");
        }
        MPI_Finalize();
        return 0;
    }

    // 初始化各进程的值为rank值
    data = rank;

    // 对各个进程进行二叉树求和
    int step;
    for (step = 1; step < size; step <<= 1) {
        if (rank % (2 * step) == 0) {
            int partner = rank + step;
            int received_data;
            if (partner < size) {
                MPI_Recv(&received_data, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                data += received_data;
            }
        } else if ((rank % step) == 0) {
            int partner = rank - step;
            MPI_Send(&data, 1, MPI_INT, partner, 0, MPI_COMM_WORLD);
            break;
        }
    }

    // 各进程求和值进行广播
    for (step = step >> 1; step > 0; step >>= 1) {
        if (rank % (2 * step) == 0) {
            int partner = rank + step;
            if (partner < size) {
                MPI_Send(&data, 1, MPI_INT, partner, 0, MPI_COMM_WORLD);
            }
        } else if ((rank % step) == 0) {
            int partner = rank - step;
            MPI_Recv(&data, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    
    MPI_Bcast(&data, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("Processor %d has the total sum: %d\n", rank, data);

    MPI_Finalize();
    return 0;
}
