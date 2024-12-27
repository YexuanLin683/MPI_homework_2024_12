#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int rank, size, node_rank, num_nodes;
    int color, key;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    num_nodes = 2; 
    color = rank % num_nodes; // new group's rank : color 
    key = rank / num_nodes; //  new group's id : key
    MPI_Comm new_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &new_comm);
    int new_rank;
    MPI_Comm_rank(new_comm, &new_rank);
    int new_size;
    MPI_Comm_size(new_comm, &new_size);

    // 输出节点分组情况
    printf("Original Rank: %d, Node Rank: %d, New Group Rank: %d, Group Size: %d\n", rank, color, new_rank, new_size);

   
    MPI_Comm_free(&new_comm);
    MPI_Finalize();
    return 0;
}

