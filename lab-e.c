#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
// 参数服务器数
#define N 2
int main()
{
    // n个进程，2个参数服务器进程，n-2个工作进程
    int comm_sz, rank, ser_num;
    int a[6] = {0, 1, 2, 3, 4, 5};
    float value, buff;
    MPI_Group world_group, ser_group;
    MPI_Comm work1, ser;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // 创建服务器子通信域
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_incl(world_group, N, a, &ser_group);
    MPI_Comm_create(MPI_COMM_WORLD, ser_group, &ser);
    if (rank / N == 0)
        MPI_Comm_size(ser, &ser_num);
    // 创建服务器与工作的通行域
    int color;
    color = rank % N;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &work1);
    int i, rank1;
    MPI_Comm_rank(work1, &rank1);
    srand((unsigned int)rank);
    // 简化模拟3轮
    for (i = 1; i < 4; i++)
    {
        if (rank1 == 0)
            value = 0.0;
        else
            value = rand() % 10;
        // 服务器进程回收工作进程数据

        MPI_Reduce(&value, &buff, 1, MPI_FLOAT, MPI_SUM, 0, work1);
        MPI_Barrier(MPI_COMM_WORLD);
        // 服务器之间通信
        if (rank1 == 0)
        {

            MPI_Allreduce(&buff, &value, 1, MPI_FLOAT, MPI_SUM, ser);
            value = value / (comm_sz - ser_num);
        }
        // 服务器与各自工作进程通信
        MPI_Bcast(&value, 1, MPI_FLOAT, 0, work1);
        printf("i = %d , my_id = %d , value = %f\n", i, rank, value);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
