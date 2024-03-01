#include "DLU.h"
#include "Eigen/Dense"
#include <vector>
#include <iostream>
#include <tuple>
#include <cstdio>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;

bool sinCerosDiagonal(MatrixXd &A) {
    for (int i = 0; i < A.cols(); i++) {
        if (A(i, i) == 0) {
            return false;
        }
    }
    return true;
}

tuple<MatrixXd, MatrixXd, MatrixXd> eigenDLU(MatrixXd &A) {
    vector<int> diagonal;
    int n = A.cols();   // Matrices cuadradas --> A.rows() = A.cols()

    MatrixXd D(n, n);
    MatrixXd L(n, n);
    MatrixXd U(n, n);

    D.setZero();
    L.setZero();
    U.setZero();


    for (int i = 0; i < n; i++) {
        D(i,i) = A(i,i);
        for (int j = 0; j < n; j++) {
            if (j > i) {
                if (A(i,j) == 0) {
                    U(i,j) = 0;         // Esto lo hago porque sino el programa coloca -0 en la matriz
                } else {                // No sé si después sigue operando de forma normal o no.
                    U(i, j) = -A(i, j);
                }
            } else if (j < i) {
                if (A(i,j) == 0) {      // Lo hago por la misma razón de antes. Mejor prevenir que curar
                    L(i,j) = 0;
                } else {
                    L(i,j) = -A(i,j);
                }
            }
        }
    }


    tuple<MatrixXd, MatrixXd, MatrixXd> dlu;
    get<0>(dlu) = D;
    get<1>(dlu) = L;
    get<2>(dlu) = U;
    return dlu;
}





