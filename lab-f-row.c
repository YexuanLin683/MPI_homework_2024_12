#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
double **A, **B, **C;
double *a, *b;
double *tmp1, *tmp2;
int rank;
int p, dg, dl, dl2;  //dg为手动输入的矩阵大小
MPI_Status status;
double f(int x, int y)
{
  double val;
  if (x == -1)
    val = tmp1[y];
  else if (x == dl)
    val = tmp2[y];
  else
    val = a[x * dg + y];
  return val;
}
// 随机初始化矩阵A
void Rand_A()
{
  int i, j;
  srand((unsigned int)time(NULL));
  for (i = 0; i < dg; i++)
  {
    for (j = 0; j < dg; j++)
    {
      A[i][j] = rand() % 100;
    }
  }
}
// 按行块分发给各进程
void Scatter_A()
{
  int i, j, k, l, p_imin, p_imax;
  for (k = 1; k < p; k++)
  {
    p_imin = k * dl;
    p_imax = p_imin + dl;
    l = 0;
    for (i = p_imin; i < p_imax; i++)
    {
      for (j = 0; j < dg; j++)
      {
        a[l] = A[i][j];
        l++;
      }
    }
    MPI_Send(a, dl2, MPI_DOUBLE, k, 1, MPI_COMM_WORLD);
  }
  for (i = 0; i < dl; i++)
  {
    for (j = 0; j < dg; j++)
    {
      a[i * dg + j] = A[i][j];
    }
  }
}
// 各进程计算值
void compu_AB()
{
  int i, j;
  memcpy(tmp1, a, dg * sizeof(double));
  memcpy(tmp2, a + (dl - 1) * dg, dg * sizeof(double));
  if (rank < p - 1)
  {
    MPI_Send(tmp2, dg, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
    MPI_Recv(tmp2, dg, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
  }
  if (rank > 0)
  {
    MPI_Send(tmp1, dg, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
    MPI_Recv(tmp1, dg, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
  }
  for (i = 0; i < dl; i++)
  {
    if ((rank == 0 && i == 0) || (rank == p - 1 && i == dl - 1))
    {
      for (j = 0; j < dg; j++)
      {
        b[i * dg + j] = 0.0;
        continue;
      }
      continue;
    }
    for (j = 0; j < dg; j++)
    {
      if (j == 0 || j == dg - 1)
      {
        b[i * dg + j] = 0.0;
        continue;
      }
      b[i * dg + j] = (f(i - 1, j) + f(i + 1, j) + f(i, j - 1) + f(i, j + 1)) / 4.0;
    }
  }
}
// 主进程收集信息
void collect_C()
{
  int i, j, k, l, p_imin, p_imax, i2, j2;
  for (i = 0; i < dl; i++)
  {
    for (j = 0; j < dg; j++)
    {
      B[i][j] = b[i * dg + j];
    }
  }
  for (k = 1; k < p; k++)
  {

    MPI_Recv(b, dl2, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &status);
    p_imin = k * dl;
    p_imax = p_imin + dl;
    i2 = 0;
    for (i = p_imin; i < p_imax; i++)
    {
      for (j = 0; j < dg; j++)
      {
        B[i][j] = b[i2 * dg + j];
      }
      i2++;
    }
  }
}
// 打印
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
  double start_time1 = MPI_Wtime();
  // 行分块
  if (dg % p != 0)
  {
    if (rank == 0)
      printf("无法完整划分\n");
    MPI_Finalize();
    exit(1);
  }
  if (argc != 2)
  {
    if (rank == 0)
      printf("error\n");
    MPI_Finalize();
    exit(1);
  }
  dl = dg / p;
  dl2 = dl * dg;
  // 获得初始分块矩阵
  a = (double *)malloc(dl2 * sizeof(double));
  b = (double *)malloc(dl2 * sizeof(double));
  tmp1 = (double *)malloc(dg * sizeof(double));
  tmp2 = (double *)malloc(dg * sizeof(double));
  for (i = 0; i < dg; i++)
    b[i] = 0.0;
  if (rank == 0)
  {
    A = (double **)malloc(dg * sizeof(double));
    B = (double **)malloc(dg * sizeof(double));
    C = (double **)malloc(dg * sizeof(double));
    for (i = 0; i < dg; i++)
    {
      A[i] = (double *)malloc(dg * sizeof(double));
      B[i] = (double *)malloc(dg * sizeof(double));
      C[i] = (double *)malloc(dg * sizeof(double));
    }
    Rand_A();
    Scatter_A();
  }
  else
    MPI_Recv(a, dl2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
  // 分别计算分块B值
  compu_AB();
  if (rank == 0)
  {
    collect_C();
    //print(A, "A\n");
    //print(B, "B\n");
  }
  else{
    MPI_Send(b, dl2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  }
  double end_time1 = MPI_Wtime();
  printf("Time for Parallel taken: %f seconds\n", end_time1 - start_time1);
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
    double end_time2 = MPI_Wtime();
    //print(C, "串行B\n");
    
    printf("Time for Serial taken: %f seconds\n", end_time2 - start_time2);
  }
  MPI_Finalize();
  return 0;
}
