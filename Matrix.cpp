#include "Matrix.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <iterator>
#include <iterator>


Matrix::Matrix(int rows, int cols) : rows(rows), cols(cols) {
    mat.resize(rows, std::list<double>(cols));
}

void Matrix::input() {
    std::cout << "Enter elements of " << rows << "x" << cols << " matrix row-wise:" << std::endl;
    for (auto& row : mat) {
        for (auto& elem : row) {
            std::cin >> elem;
        }
    }
}

void Matrix::display() const {
    for (const auto& row : mat) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

Matrix Matrix::transpose() const {
    Matrix result(cols, rows);
    auto resultIt = result.mat.begin();
    for (const auto& row : mat) {
        auto resultRowIt = resultIt->begin();
        for (auto colIt = row.begin(); colIt != row.end(); ++colIt, ++resultRowIt) {
            *resultRowIt = *colIt;
        }
        ++resultIt;
    }
    return result;
}

Matrix Matrix::pseudoInverse() const {
    Matrix transposed = transpose();
    Matrix product = transposed.multiply(*this);

    if (product.determinant(product.mat, product.rows) != 0) {
        Matrix inverseProduct = product.inverseGauss();
        Matrix pseudoInverse = inverseProduct.multiply(transposed);
        return pseudoInverse;
    }
}

Matrix Matrix::inverseGauss() const {
    if (!isSquare()) {
        throw std::logic_error("Non-square matrices cannot be inverted using Gauss.");
    }

    int n = rows;
    std::list<std::list<double>> a = mat;
    std::list<std::list<double>> inv(n, std::list<double>(n, 0));

    // Initialize the identity matrix
    auto invRowIt = inv.begin();
    for (int i = 0; i < n; ++i, ++invRowIt) {
        auto invColIt = invRowIt->begin();
        std::advance(invColIt, i); // Move to the i-th element
        *invColIt = 1; // Set it to 1
    }

    // Perform Gauss-Jordan elimination
    auto rowIt = a.begin();
    for (int i = 0; i < n; ++i, ++rowIt) {
        double diagElement = *std::next(rowIt->begin(), i);
        if (diagElement == 0) {
            throw std::logic_error("Matrix is singular and cannot be inverted.");
        }

        auto colIt = rowIt->begin();
        auto invColIt = std::next(inv.begin(), i)->begin();
        for (int j = 0; j < n; ++j, ++colIt, ++invColIt) {
            *colIt /= diagElement;
            *invColIt /= diagElement;
        }

        auto kRowIt = a.begin();
        auto invKRowIt = inv.begin();
        for (int k = 0; k < n; ++k, ++kRowIt, ++invKRowIt) {
            if (k != i) {
                double factor = *std::next(kRowIt->begin(), i);
                auto kColIt = kRowIt->begin();
                auto iColIt = rowIt->begin();
                auto invKColIt = invKRowIt->begin();
                auto invIColIt = std::next(inv.begin(), i)->begin();

                for (int j = 0; j < n; ++j, ++kColIt, ++iColIt, ++invKColIt, ++invIColIt) {
                    *kColIt -= factor * (*iColIt);
                    *invKColIt -= factor * (*invIColIt);
                }
            }
        }
    }

    Matrix inverseMatrix(n, n);
    inverseMatrix.mat = inv;
    return inverseMatrix;
}

double Matrix::determinant(const std::list<std::list<double>>& m, int n) const {
    if (n == 1) {
        return m.front().front();
    }
    if (n == 2) {
        return m.front().front() * m.back().back() - m.front().back() * m.back().front();
    }

    double det = 0;
    auto mIt = m.begin();
    for (int x = 0; x < n; ++x, ++mIt) {
        std::list<std::list<double>> submatrix = getCofactor(m, 0, x, n);
        det += ((x % 2 == 0 ? 1 : -1) * m.front().front() * determinant(submatrix, n - 1));
    }

    return det;
}

std::list<std::list<double>> Matrix::adjoint() const {
    if (!isSquare()) {
        throw std::logic_error("Non-square matrices do not have an adjoint.");
    }

    int n = rows;
    if (n == 1) {
        return { {1} };
    }

    std::list<std::list<double>> adj(n, std::list<double>(n));
    auto adjIt = adj.begin();
    for (int i = 0; i < n; ++i, ++adjIt) {
        auto adjRowIt = adjIt->begin();
        for (int j = 0; j < n; ++j, ++adjRowIt) {
            std::list<std::list<double>> cofactor = getCofactor(mat, i, j, n);
            *adjRowIt = determinant(cofactor, n - 1) * ((i + j) % 2 == 0 ? 1 : -1);
        }
    }

    return adj;
}

std::list<std::list<double>> Matrix::getCofactor(const std::list<std::list<double>>& m, int p, int q, int n) const {
    std::list<std::list<double>> temp(n - 1, std::list<double>(n - 1));
    auto tempIt = temp.begin();
    for (auto rowIt = m.begin(); rowIt != m.end(); ++rowIt) {
        if (std::distance(m.begin(), rowIt) != p) {
            auto tempRowIt = tempIt->begin();
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
                if (std::distance(rowIt->begin(), colIt) != q) {
                    *tempRowIt = *colIt;
                    ++tempRowIt;
                }
            }
            ++tempIt;
        }
    }
    return temp;
}

Matrix Matrix::inverseCramer() const {
    if (!isSquare()) {
        throw std::logic_error("Non-square matrices cannot be inverted using Cramer.");
    }

    int n = rows;
    double det = determinant(mat, n);
    if (det == 0) {
        throw std::logic_error("Matrix is singular and cannot be inverted.");
    }

    std::list<std::list<double>> adj = adjoint();
    std::list<std::list<double>> inv(n, std::list<double>(n));
    auto invIt = inv.begin();
    for (auto& adjRow : adj) {
        auto invRowIt = invIt->begin();
        for (auto& adjElem : adjRow) {
            *invRowIt = adjElem / det;
            ++invRowIt;
        }
        ++invIt;
    }

    Matrix inverseMatrix(n, n);
    inverseMatrix.mat = inv;
    return inverseMatrix;
}

Matrix Matrix::multiply(const Matrix& other) const {
    if (cols != other.rows) {
        throw std::logic_error("Matrices dimensions do not match for multiplication.");
    }

    Matrix result(rows, other.cols);
    auto resultIt = result.mat.begin();
    for (const auto& row : mat) {
        auto resultRowIt = resultIt->begin();
        for (int j = 0; j < other.cols; ++j, ++resultRowIt) {
            *resultRowIt = 0;
            auto otherRowIt = other.mat.begin();
            auto rowIt = row.begin();
            for (int k = 0; k < cols; ++k, ++rowIt, ++otherRowIt) {
                auto otherColIt = otherRowIt->begin();
                std::advance(otherColIt, j);
                *resultRowIt += *rowIt * *otherColIt;
            }
        }
        ++resultIt;
    }

    return result;
}

void Matrix::roundSmallValues(double threshold) {
    for (auto& row : mat) {
        for (auto& elem : row) {
            if (std::abs(elem) < threshold) {
                elem = 0;
            }
        }
    }
}

bool Matrix::isIdentity() const {
    if (!isSquare()) {
        return false;
    }

    auto rowIt = mat.begin();
    for (int i = 0; i < rows; ++i, ++rowIt) {
        auto colIt = rowIt->begin();
        for (int j = 0; j < cols; ++j, ++colIt) {
            if ((i == j && std::abs(*colIt - 1) >= 1e-8 && std::abs(*colIt - 1) > 0) || (i != j && std::abs(*colIt) >= 1e-8 && std::abs(*colIt - 1) > 0)) {
                return false;
            }
        }
    }

    return true;
}