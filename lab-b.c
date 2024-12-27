#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define DATA_SIZE 10  // 每个进程要发送的数据大小

int main(int argc, char** argv) {
    int rank, size;
    int *send_data, *recv_data;
    int i, j;

    double start_time1, end_time1;
    double start_time2, end_time2;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    send_data = (int*) malloc(DATA_SIZE * size * sizeof(int));
    recv_data = (int*) malloc(DATA_SIZE * size * sizeof(int));

    // 初始化发送数据
    for (i = 0; i < DATA_SIZE*size; i++) {
        send_data[i] = rank * DATA_SIZE * size  + i;  // 使得每个进程的数据不同
    }

    // 模拟 MPI_Alltoall，使用 MPI_Send 和 MPI_Recv
    // 对于每个目标进程，发送数据并接收数据

    start_time1 = MPI_Wtime();
    for (i = 0; i < size; i++) {
       if (1) {
        MPI_Send(send_data+i*DATA_SIZE, DATA_SIZE, MPI_INT, i, 0, MPI_COMM_WORLD);
    }}

    for (i = 0; i < size; i++) {
       if (1) {
        MPI_Recv(&recv_data[i * DATA_SIZE], DATA_SIZE, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }}
    end_time1 = MPI_Wtime();

    // 每个进程输出接收到的数据
    printf("Here are the MPI_Send & MPI_Recv:\n");
    printf("Process %d received data: ", rank);
    for (i = 0; i < size * DATA_SIZE; i++) {
        printf("%d ", recv_data[i]);
    }
    printf("\n");
    printf("Time taken: %f seconds\n", end_time1 - start_time1);
    printf("***************************");
    printf("\n");

    // 使用 MPI_Alltoall 执行全体对全体的通信
    start_time2 = MPI_Wtime();
    MPI_Alltoall(send_data, DATA_SIZE, MPI_INT, recv_data, DATA_SIZE, MPI_INT, MPI_COMM_WORLD);
    end_time2 = MPI_Wtime();
    
    // 每个进程输出接收到的数据
    printf("Here are the MPI_Alltoall:\n");
    printf("Process %d received data: ", rank);
    for (i = 0; i < size * DATA_SIZE; i++) {
        printf("%d ", recv_data[i]);
    }
    printf("\n");
    printf("Time taken: %f seconds\n", end_time2 - start_time2);
    printf("***************************");
    printf("\n");

    free(send_data);
    free(recv_data);
    MPI_Finalize();
    return 0;
}

