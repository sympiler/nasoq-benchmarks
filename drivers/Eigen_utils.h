//
// Created by kazem on 4/21/19.
//

#ifndef PROJECT_EIGEN_UTILS_H
#define PROJECT_EIGEN_UTILS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
/*
 * Converting matrix market format to Eigen sparse matrix
 */
int convert_mtx_to_eigen(int n_row, int n_col,
                         int *Ap, int *Ai, double *Ax,
                         Eigen::SparseMatrix<double> &A){
 typedef Triplet<double> T;
 std::vector<T> tripletList;
 for (int i = 0; i < n_col; ++i) {
  for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
   tripletList.push_back(T(Ai[j],i,Ax[j]));
  }
 }
 A.setFromTriplets(tripletList.begin(),tripletList.end());
 return 1;
}

#endif //PROJECT_EIGEN_UTILS_H
