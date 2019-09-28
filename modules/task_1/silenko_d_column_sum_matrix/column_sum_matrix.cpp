// Copyright 2019 SIlenko Dmitrii

#include <mpi.h>
#include <iostream>
#include <random>
#include <ctime>
#include <numeric>
#include <vector>
#include <stdexcept>
#include "../../../modules/task_1/silenko_d_column_sum_matrix/column_sum_matrix.h"

std::vector<std::vector <int>> getRandomMatrixE(int n, int m) {
  std::vector <std::vector <int>> Matrix(n, std::vector <int>(m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      Matrix[i][j] = 1;
    }
  }
  return Matrix;
}

std::vector<std::vector <int>> getRandomMatrixO(int n, int m) {
  std::vector <std::vector <int>> Matrix(n, std::vector <int>(m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      Matrix[i][j] = j + 1;
    }
  }
  return Matrix;
}

std::vector <std::vector <int>> TransposedMatrix(const std::vector <std::vector <int>> &a, int n, int m) {
  std::vector <std::vector <int>> transposed(m, std::vector <int>(n));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      transposed[j][i] = a[i][j];
  return transposed;
}

std::vector <int> ColumnSumMatrix(const std::vector <std::vector <int>> &a, int n, int m) {
  std::vector <std::vector <int>> transposed(m, std::vector <int>(n));
  std::vector <int> ansv(m);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  transposed = TransposedMatrix(a, n, m);

  if (rank == 0) {
    for (int i = 0; i < m; i++) {
      if (i % size) {
        MPI_Send(&transposed[i][0], n, MPI_INT, i % size, 5, MPI_COMM_WORLD);
      }
    }
  }
  MPI_Status status;
  std::vector <int> b(n);
  std::vector <int> ans(m);
  if (rank == 0) {
    for (int i = 0; i < m; i += size) {
      for (int j = 0; j < n; j++) {
        ans[i] += transposed[i][j];
      }
    }
  } else {
    for (int i = rank; i < m; i += size) {
      MPI_Recv(&b[0], n, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);
      for (int j = 0; j < n; j++)
        ans[i] += b[j];
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    for (int i = 0; i < m; i++) {
      if (i % size == 0)
        ansv[i] = ans[i];
      else
        MPI_Recv(&ansv[i], 1, MPI_INT, i % size, 9, MPI_COMM_WORLD, &status);
    }
  } else {
    for (int i = rank; i < m; i += size) {
      MPI_Send(&ans[i], 1, MPI_INT, 0, 9, MPI_COMM_WORLD);
    }
  }
  return ansv;
}
