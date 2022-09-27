#include "solver.hpp"

Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b) {
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> cholesky;
    cholesky.compute(A);
    return cholesky.solve(b);
}

std::vector<Vec3> linear_solve(const SpMat<Mat3x3>& A, const std::vector<Vec3>& b) {
    int n = b.size();
    Eigen::SparseMatrix<double> At(3 * n, 3 * n);
    for (int i = 0; i < A.rows.size(); i++)
        for (int j = 0; j < A.rows[i].indices.size(); j++)
            for (int k = 0; k < 3; k++)
                for (int h = 0; h < 3; h++)
                    At.coeffRef(3 * i + k, 3 * A.rows[i].indices[j] + h) = A.rows[i].entries[j](k, h);
    
    Eigen::VectorXd bt(3 * n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < 3; j++)
            bt(3 * i + j) = b[i][j];

    Eigen::VectorXd x = solve(At, bt);

    std::vector<Vec3> ans(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < 3; j++)
            ans[i][j] = x(3 * i + j);
    
    return ans;
}

std::vector<double> linear_solve(const SpMat<double>& A, const std::vector<double>& b) {
    int n = b.size();
    Eigen::SparseMatrix<double> At(n, n);
    for (int i = 0; i < A.rows.size(); i++)
        for (int j = 0; j < A.rows[i].indices.size(); j++)
            At.coeffRef(i, A.rows[i].indices[j]) = A.rows[i].entries[j];
    
    Eigen::VectorXd bt(n);
    for (int i = 0; i < n; i++)
        bt(i) = b[i];
    
    Eigen::VectorXd x = solve(At, bt);

    std::vector<double> ans(n);
    for (int i = 0; i < n; i++)
        ans[i] = x(i);
    
    return ans;
}