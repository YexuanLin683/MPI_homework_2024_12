#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
// 全局变量
float **A, **B, **C;
float *a, *b, *c, *tmp_a, *tmp_b;
int dg, dl, dl2, p, sp;   //dg是你手动输入的参数 控制矩阵的大小
int my_rank, my_row, my_col;
MPI_Comm world1;
int rank1;
MPI_Status status;
// 初始化矩阵
void random_A_B()
{
    int i, j;
    srand((unsigned int)time(NULL));
    for (i = 0; i < dg; i++)
    {
        for (j = 0; j < dg; j++)
        {
            A[i][j] = rand() % 100;
            B[i][j] = rand() % 100;
            C[i][j] = 0.0;
        }
    }
}
// 矩阵初始分块
void scatter_A_B()
{
    int i, j, k, l;
    int p_imin, p_imax, p_jmin, p_jmax;
    for (k = 0; k < p; k++)
    {
        p_jmin = (k % sp) * dl;
        p_jmax = (k % sp + 1) * dl - 1;
        p_imin = k / sp * dl;
        p_imax = (k / sp + 1) * dl - 1;
        l = 0;
        for (i = p_imin; i <= p_imax; i++)
        {
            for (j = p_jmin; j <= p_jmax; j++)
            {
                tmp_a[l] = A[i][j];
                tmp_b[l] = B[i][j];
                l++;
            }
        }
        if (k == 0)
        {
            memcpy(a, tmp_a, dl2 * sizeof(float));
            memcpy(b, tmp_b, dl2 * sizeof(float));
        }
        else
        {

            MPI_Send(tmp_a, dl2, MPI_FLOAT, k, 1, MPI_COMM_WORLD);

            MPI_Send(tmp_b, dl2, MPI_FLOAT, k, 2, MPI_COMM_WORLD);
        }
    }
}
// 广播A块、计算、平移B块,k用来控制取的列数
void comp_ABC(int k)
{
    int i, j, m;
    if ((my_row + k) % sp == rank1)
    {
        memcpy(tmp_a, a, dl2 * sizeof(float));
    }
    MPI_Bcast(tmp_a, dl2, MPI_FLOAT, (my_row + k) % sp, world1);
    for (i = 0; i < dl; i++)
    {
        for (j = 0; j < dl; j++)
        {
            for (m = 0; m < dl; m++)
            {
                c[i * dl + j] += tmp_a[i * dl + m] * b[m * dl + j];
            }
        }
    }
    memcpy(tmp_b, b, dl2 * sizeof(float));
    MPI_Sendrecv(tmp_b, dl2, MPI_FLOAT, (my_row - 1 + sp) % sp * sp + my_col, 0, b, dl2, MPI_FLOAT, (my_row + 1) % sp * sp + my_col, 0, MPI_COMM_WORLD, &status);
    return;
}
// 回收结果
void collect_C()
{
    int i, j, k, i2, j2;
    int p_imin, p_imax, p_jmin, p_jmax;
    for (i = 0; i < dl; i++)
    {
        for (j = 0; j < dl; j++)
        {
            C[i][j] = c[i * dl + j];
        }
    }
    for (k = 1; k < p; k++)
    {

        MPI_Recv(c, dl2, MPI_FLOAT, k, 1, MPI_COMM_WORLD, &status);
        p_jmin = (k % sp) * dl;
        p_jmax = (k % sp + 1) * dl - 1;
        p_imin = k / sp * dl;
        p_imax = (k / sp + 1) * dl - 1;
        i2 = 0;
        for (i = p_imin; i <= p_imax; i++)
        {
            j2 = 0;
            for (j = p_jmin; j <= p_jmax; j++)
            {
                C[i][j] = c[i2 * dl + j2];
                j2++;
            }
            i2++;
        }
    }
}
// 打印矩阵
void print(float **a, char *str)
{
    printf("%s", str);
    int i, j;
    for (i = 0; i < dg; i++)
    {
        for (j = 0; j < dg; j++)
        {
            printf("%4.0f ", a[i][j]);
        }
        printf("\n");
    }
}
int main(int argc, char *argv[])
{
    int i, color;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    sp = sqrt(p);
    // 判断进程数是否为完全平方数
    if (sp * sp != p)
    {
        if (my_rank == 0)
            printf("p error\n");
        MPI_Finalize();
        exit(1);
    }
    if (argc != 2)
    {
        if (my_rank == 0)
            printf("dg error\n");
        MPI_Finalize();
        exit(1);
    }
    double start_time1 = MPI_Wtime();
    // 按行分组
    color = my_rank / sp;
    MPI_Comm_split(MPI_COMM_WORLD, color, my_rank, &world1);
    MPI_Comm_rank(world1, &rank1);
    dg = atoi(argv[1]);
    if (dg % sp != 0)
    {
        if (my_rank == 0)
            printf("dg / sp error!\n");
        MPI_Finalize();
        exit(1);
    }
    dl = dg / sp;
    dl2 = dl * dl;
    my_col = my_rank % sp;
    my_row = my_rank / sp;
    // 获得初始分块矩阵
    a = (float *)malloc(dl2 * sizeof(float));
    b = (float *)malloc(dl2 * sizeof(float));
    c = (float *)malloc(dl2 * sizeof(float));
    for (i = 0; i < dl2; i++)
        c[i] = 0.0;
    tmp_a = (float *)malloc(dl2 * sizeof(float));
    tmp_b = (float *)malloc(dl2 * sizeof(float));
    if (my_rank == 0)
    {
        A = (float **)malloc(dg * sizeof(float *));
        B = (float **)malloc(dg * sizeof(float *));
        C = (float **)malloc(dg * sizeof(float *));
        for (i = 0; i < dg; i++)
        {
            A[i] = (float *)malloc(dg * sizeof(float));
            B[i] = (float *)malloc(dg * sizeof(float));
            C[i] = (float *)malloc(dg * sizeof(float));
        }
        random_A_B();
        scatter_A_B();
    }
    else
    {

        MPI_Recv(a, dl2, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);

        MPI_Recv(b, dl2, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status);
    }
    // 分别计算、按规则广播，平移B
    for (i = 0; i < sp; i++)
    {
        comp_ABC(i);
    }
    if (my_rank == 0)
    {
        collect_C();
        //print(A, "A矩阵\n");
        //print(B, "B矩阵\n");
        //print(C, "FOX算法_C\n");
    }
    else{
        MPI_Send(c, dl2, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
    }
    double end_time1 = MPI_Wtime();
    if (my_rank == 0)
    {
        int j, k;
        double start_time2 = MPI_Wtime();
        for (i = 0; i < dg; i++)
        {
            for (j = 0; j < dg; j++)
            {
                C[i][j] = 0.0;
                for (k = 0; k < dg; k++)
                {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        double end_time2 = MPI_Wtime();
        //print(C, "串行结果_C\n");
        printf("Time for Parallel taken: %f seconds\n", end_time1 - start_time1);
        printf("Time for Serial taken: %f seconds\n", end_time2 - start_time2);
    }
    MPI_Finalize();
    return 0;
}
