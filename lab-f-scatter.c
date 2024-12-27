#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

double **A, **B, **C;
double *a, *b;
double *tmp_1, *tmp_2;
int ql;
int rank;
int p, sp, dg, dl, dl2; // dg为手动输入的矩阵大小
int my_col, my_row;
MPI_Status status;
double f(int x, int y)
{
    if (x == -1 || x == dl)
        return tmp_1[y];
    else if (y == -1 || y == dl)
        return tmp_2[x];
    else
        return a[x * dl + y];
}
void comp_B(int x)
{
    int i, j;
    for (i = 0; i < dl; i++)
    {
        if (my_row == 0 && i == 0 || my_row == sp - 1 && i == dl - 1)
        {
            for (j = 0; j < dl; j++)
                b[i * dl + j] += 0.0;
            continue;
        }
        for (j = 0; j < dl; j++)
        {
            if (my_col == 0 && j == 0 || my_col == sp - 1 && j == dl - 1)
                b[i * dl + j] += 0.0;
            else
                b[i * dl + j] += (f(i + x, j) + f(i, j + x)) / 4.0;
        }
    }
}
void init_A()
{
    int i, j;
    srand((unsigned int)time(NULL));
    for (i = 0; i < dg; i++)
    {
        for (j = 0; j < dg; j++)
        {
            A[i][j] = rand() % 100;
            B[i][j] = 0.0;
        }
    }
    for (i = 0; i < dl2; i++)
        b[i] = 0.0;
}
void scatter_A()
{
    if (rank == 0)
    {
        int m, n;
        int i, j, k;
        int i_min, j_min, i_max, j_max;
        for (i = 1; i < p; i++)
        {
            i_min = i / sp * dl;
            i_max = i_min + dl;
            j_min = i % sp * dl;
            j_max = j_min + dl;
            m = 0;
            for (j = i_min; j < i_max; j++)
            {
                n = 0;
                for (k = j_min; k < j_max; k++)
                {
                    a[m * dl + n] = A[j][k];
                    n++;
                }
                m++;
            }
            MPI_Send(a, dl2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        for (i = 0; i < dl; i++)
        {
            for (j = 0; j < dl; j++)
            {
                a[i * dl + j] = A[i][j];
            }
        }
    }
    else
    {

        MPI_Recv(a, dl2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }
}
void send_A()
{
    int i;
    double *ql;
    ql = (double *)malloc(dl * sizeof(double));
    for (i = 0; i < dl; i++)
    {
        tmp_1[i] = a[i];      // 上
        tmp_2[i] = a[i * dl]; // 左
    }
    memcpy(ql, tmp_2, dl * sizeof(double));
    MPI_Sendrecv(ql, dl, MPI_DOUBLE, my_row * sp + (my_col - 1 + sp) % sp, 0, tmp_2, dl, MPI_DOUBLE, my_row * sp + (my_col + 1) % sp, 0, MPI_COMM_WORLD, &status);
    memcpy(ql, tmp_1, dl * sizeof(double));
    MPI_Sendrecv(ql, dl, MPI_DOUBLE, (my_row - 1 + sp) % sp * sp + my_col, 0, tmp_1, dl, MPI_DOUBLE, (my_row + 1) % sp * sp + my_col, 0, MPI_COMM_WORLD, &status);
    comp_B(1);
    for (i = 0; i < dl; i++)
    {
        tmp_1[i] = a[dl2 - dl + i];     // 下
        tmp_2[i] = a[(i + 1) * dl - 1]; // 右
    }
    memcpy(ql, tmp_2, dl * sizeof(double));
    MPI_Sendrecv(ql, dl, MPI_DOUBLE, my_row * sp + (my_col + 1) % sp, 0, tmp_2, dl, MPI_DOUBLE, my_row * sp + (my_col - 1 + sp) % sp, 0, MPI_COMM_WORLD, &status);
    memcpy(ql, tmp_1, dl * sizeof(double));
    MPI_Sendrecv(ql, dl, MPI_DOUBLE, (my_row + 1) % sp * sp + my_col, 0, tmp_1, dl, MPI_DOUBLE, my_col + (my_row - 1 + sp) % sp * sp, 0, MPI_COMM_WORLD, &status);
    comp_B(-1);
}
void gather_B()
{
    int m, n;
    if (rank == 0)
    {
        int i, j, k;
        int i_min, j_min, i_max, j_max;
        for (i = 0; i < dl; i++)
        {
            for (j = 0; j < dl; j++)
            {
                B[i][j] = b[i * dl + j];
            }
        }
        for (i = 1; i < p; i++)
        {

            MPI_Recv(b, dl2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
            i_min = i / sp * dl;
            i_max = i_min + dl;
            j_min = i % sp * dl;
            j_max = j_min + dl;
            m = 0;
            for (j = i_min; j < i_max; j++)
            {
                n = 0;
                for (k = j_min; k < j_max; k++)
                {
                    B[j][k] = b[m * dl + n];
                    n++;
                }
                m++;
            }
        }
    }
    else
    {
        MPI_Send(b, dl2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
}
void print(double **a, char *str)
{
    printf("%s", str);
    int i, j;
    for (i = 0; i < dg; i++)
    {
        for (j = 0; j < dg; j++)
        {
            printf("%3.0f ", a[i][j]);
        }
        printf("\n");
    }
}
int main(int argc, char *argv[])
{
    int i, j;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    dg = atoi(argv[1]);
    // 块划分
    sp = sqrt(p);
    double start_time1 = MPI_Wtime();
    if (sp * sp != p)
    {
        if (rank == 0)
            printf("p error!\n");
        MPI_Finalize();
        exit(1);
    }
    if (argc != 2)
    {
        if (rank == 0)
            printf("dg number error!\n");
        MPI_Finalize();
        exit(1);
    }
    if (dg % sp != 0)
    {
        if (rank == 0)
            printf("dg / sp error!\n");
        MPI_Finalize();
        exit(1);
    }
    // 获得初始分块矩阵
    dl = dg / sp;
    dl2 = dl * dl;
    my_col = rank % sp;
    my_row = rank / sp;
    a = (double *)malloc(dl2 * sizeof(double));
    b = (double *)malloc(dl2 * sizeof(double));
    tmp_1 = (double *)malloc(dl * sizeof(double));
    tmp_2 = (double *)malloc(dl * sizeof(double));
    // 初始化矩阵A,B
    if (rank == 0)
    {
        A = (double **)malloc(dg * sizeof(double *));
        B = (double **)malloc(dg * sizeof(double *));
        C = (double **)malloc(dg * sizeof(double *));
        for (i = 0; i < dg; i++)
        {
            A[i] = (double *)malloc(dg * sizeof(double));
            B[i] = (double *)malloc(dg * sizeof(double));
            C[i] = (double *)malloc(dg * sizeof(double));
        }
        init_A();
    }
    // 分发初始化矩阵
    scatter_A();
    send_A();
    // 收集结果分块
    gather_B();
    double end_time1 = MPI_Wtime();
    if (rank == 0)
    {
        //print(B, "B\n");
        printf("Time for Parallel taken: %f seconds\n", end_time1 - start_time1);
    }
    
    // 串行计算
    if (rank == 0)
    {

    double start_time2 = MPI_Wtime();
    for (int i = 0; i < dg; i++) {
      for (int j = 0; j < dg; j++) {
            C[i][j] = 0;
        }
    }
        for (i = 1; i < dg - 1; i++)
        {
            for (j = 1; j < dg - 1; j++)
            {
                C[i][j] = (A[i - 1][j] + A[i + 1][j] + A[i][j - 1] + A[i][j + 1]) / 4.0;
            }
        }
        //print(A, "A\n");
        // print(C, "串行B\n");
    double end_time2 = MPI_Wtime();

    printf("Time for Serial taken: %f seconds\n", end_time2 - start_time2);
    }

    MPI_Finalize();
    return 0;
}
