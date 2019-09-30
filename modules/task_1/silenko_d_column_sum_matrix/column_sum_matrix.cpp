// Copyright 2019 SIlenko Dmitrii

#include <mpi.h>
#include <iostream>
#include <random>
#include <ctime>
#include <numeric>
#include <vector>
#include <stdexcept>
#include "../../../modules/task_1/silenko_d_column_sum_matrix/column_sum_matrix.h"

std::vector<std::vector <int>> getRandomMatrixE(const int n, const int m) {
  if (n <= 0) {
    throw "Wrong rows";
  }
  else if (m <= 0) {
    throw "wrong columns";
  }
  std::vector <std::vector <int>> Matrix(n, std::vector <int>(m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      Matrix[i][j] = 1;
    }
  }
  return Matrix;
}

std::vector<std::vector <int>> getRandomMatrixO(const int n, const int m) {
  if (n <= 0) {
    throw "Wrong rows";
  }
  else if (m <= 0) {
    throw "wrong columns";
  }
  std::vector <std::vector <int>> Matrix(n, std::vector <int>(m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      Matrix[i][j] = j + 1;
    }
  }
  return Matrix;
}

std::vector <std::vector <int>> TransposedMatrix(const std::vector <std::vector <int>> &a, const int n, const int m) {
  std::vector <std::vector <int>> transposed(m, std::vector <int>(n));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      transposed[j][i] = a[i][j];
  return transposed;
}

std::vector <int> ColumnSumMatrix(const std::vector <std::vector <int>> &a, const int n, const int m) {

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector <int> ans(m);
  MPI_Status status;
  int error;

  if (rank == 0) {
    if (a.size() != n || a[0].size() != m) {
      error = -1;
    }
    else if (n <= 0) {
      error = -2;
    }
    else if (m <= 0) {
      error = -3;
    }
    else {
      error = 0;
    }
    for (int i = 1; i < size; ++i)
      MPI_Send(&error, 1, MPI_INT, i, 8, MPI_COMM_WORLD);
  }
  else {
    MPI_Recv(&error, 1, MPI_INT, 0, 8, MPI_COMM_WORLD, &status);
  }

  switch (error) {
  case 0:
    break;
  case -1:
    throw std::runtime_error("Size doesn't match description");
  case -2:
    throw std::runtime_error("Number of rows 0");
  case -3:
    throw std::runtime_error("Nubmer of columns 0");
  }

  if (rank == 0) {
    std::vector <std::vector <int>> transposed(m, std::vector <int>(n));
    transposed = TransposedMatrix(a, n, m);
    for (int i = 0; i < m; i++) {
      if (i % size) {
        MPI_Send(&transposed[i][0], n, MPI_INT, i % size, 5, MPI_COMM_WORLD);
      }
    }
    for (int i = 0; i < m; i += size) {
      for (int j = 0; j < n; j++) {
        ans[i] += transposed[i][j];
      }
    }
  }
  std::vector <int> b(n);
  if (rank != 0) {
    for (int i = rank; i < m; i += size) {
      MPI_Recv(&b[0], n, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);
      for (int j = 0; j < n; j++)
        ans[i] += b[j];
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    for (int i = 0; i < m; i++) {
      if (i % size != 0)
        MPI_Recv(&ans[i], 1, MPI_INT, i % size, 9, MPI_COMM_WORLD, &status);
    }
  }
  else {
    for (int i = rank; i < m; i += size) {
      MPI_Send(&ans[i], 1, MPI_INT, 0, 9, MPI_COMM_WORLD);
    }
  }
  return ans;
}
