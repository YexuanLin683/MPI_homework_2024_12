#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    int rank, size;
    int data;
    float start, end;
    float start1, end1;
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
    
    start = MPI_Wtime();
    // 初始化每个进程的数据为rank值
    data = rank;

    // 蝶式规约求和
    int step;
    for (step = 1; step < size; step <<= 1) {
        int partner = rank ^ step;
        int received_data;

        // 交换各部分的数据
        MPI_Sendrecv(&data, 1, MPI_INT, partner, 0,
                     &received_data, 1, MPI_INT, partner, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // 累加求和
        data += received_data;
    }

    // 最后 打印所有进程的获得值
    printf("Processor %d has the total sum: %d\n", rank, data);
    
    end = MPI_Wtime();
    if(rank == 0){
    printf("MPI sum time: %f seconds\n", end-start);
    }

    //int i; 
    //int sum_serial = 0;
    //start1 = MPI_Wtime();
    //for(i=0;i<size;i++){
    // sum_serial = sum_serial + i;
    //}
    //end1 = MPI_Wtime();
    //if(rank == 0){
    //printf("The serial sum is %d\n",sum_serial);
    //printf("The serial sum time : %f\n",end1 -start1);
    //}
    MPI_Finalize();
    return 0;
}

