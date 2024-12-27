#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    int rank, size, node_rank, num_nodes;
    int color, key;
    int root_value = 100; // root进程发送的值

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    num_nodes = 2; 
    color = rank % num_nodes; // 每个节点为一个颜色
    key = rank / num_nodes; // 进程的原始排名作为 key
    MPI_Comm new_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &new_comm);
    int new_rank;
    MPI_Comm_rank(new_comm, &new_rank);
    int new_size;
    MPI_Comm_size(new_comm, &new_size);

    printf("Original Rank: %d, Node Rank: %d, New Group Rank: %d, Group Size: %d\n", rank, color, new_rank, new_size);
    MPI_Barrier(MPI_COMM_WORLD);
    int root = 0 , value;
    if (new_rank == 0){
       value = 210;
    }
    else{
       value = 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (new_rank == 0){
    printf("Node %d's root is sending the value %d\n", color, value);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Original Rank %d but Rank %d in Node %d has the orginal value %d\n",rank, key , color , value);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&value, 1, MPI_INT, root, new_comm);
    printf("Original Rank %d but Rank %d in Node %d has received the value %d\n",rank, key , color , value);
    MPI_Comm_free(&new_comm);
    MPI_Finalize();
    return 0;
}

