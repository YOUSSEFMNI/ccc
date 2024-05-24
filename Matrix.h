#ifndef MATRIX_H
#define MATRIX_H

#include <list>
#include <iostream>

class Matrix {
private:
    std::list<std::list<double>> mat;
    int rows;
    int cols;

public:
    Matrix(int rows, int cols);
    void input();
    void display() const;
    Matrix inverseGauss() const;
    Matrix inverseCramer() const;
    double determinant(const std::list<std::list<double>>& m, int n) const;
    std::list<std::list<double>> adjoint() const;
    std::list<std::list<double>> getCofactor(const std::list<std::list<double>>& m, int p, int q, int n) const;
    Matrix multiply(const Matrix& other) const;
    bool isIdentity() const;
    bool isSquare() const { return rows == cols; }
    void roundSmallValues(double threshold = 1e-8);
    Matrix transpose() const;
    Matrix pseudoInverse() const;
    int getRows() const { return rows; }
    int getCols() const { return cols; }
};

#endif // MATRIX_H