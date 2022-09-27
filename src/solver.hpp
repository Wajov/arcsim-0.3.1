#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>

#include "sparse.hpp"

Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b);

std::vector<Vec3> linear_solve(const SpMat<Mat3x3>& A, const std::vector<Vec3>& b);

std::vector<double> linear_solve(const SpMat<double>& A, const std::vector<double>& b);

#endif