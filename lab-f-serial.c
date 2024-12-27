#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 8  // 矩阵大小

int* compute(int A[N][N]) {
    int* B = (int*)malloc(N * N * sizeof(int));
    if (B == NULL) {
        // 处理内存分配失败
        return NULL;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[i*N + j] = 0;
        }
    }
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            B[i*N + j] = (A[i-1][j] + A[i][j+1] + A[i+1][j] + A[i][j-1]) / 4.0;
        }
    }
    return B;
}

int main() {
    int Aa[N][N];
    int* B_result;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Aa[i][j] = rand() % 100 ;  // 初始化0-100随机数
        }
    }
      for (int i = 0; i < N ; i++) {
        for (int j = 0; j < N ; j++) {
            printf("%d ", Aa[i][j]);
          }
      printf("\n");
      }
    

    // 记录开始时间
    clock_t start_time = clock();

    B_result = compute(Aa);

    // 记录结束时间
    clock_t end_time = clock();
    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Total computation time (serial): %f seconds\n", total_time);
      for (int i = 0; i < N ; i++) {
          for (int j = 0; j < N; j++) {
              printf("%d ", B_result[i*N+j]);
          }
      printf("\n");
      }
    free(B_result);
    return 0;
}

