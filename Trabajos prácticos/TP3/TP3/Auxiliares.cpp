#include "Eigen/Dense"
#include <vector>
#include <fstream>
#include <iostream>
#include <tuple>
#include <cstdio>
#include "Auxiliares.h"


using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;

MatrixXd read_txt(string matriz){
    ifstream fin(matriz);
    // Leemos la matriz
    int nrows, ncols;
    fin >> nrows >> ncols;

    MatrixXd A(nrows, ncols);
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            fin >> A(i, j);
        }
    }

    fin.close();

    return A;
}

double metodoPotencia(MatrixXd &A, int N, double e) {

    int size = A.cols();
    VectorXd x(size);
    x[0] = 1; // para asegurarnos de que no es nulo
    for (int i = 1; i < size; ++i ) {
        x[i] = rand();
    }

    VectorXd res (size);
    double lambda1 = x.transpose() * A * x;
    double lambda2 = x.transpose() * A * x;
    int i = 0;
    do {
        lambda1 = lambda2;
        res = A * x;
        res.normalize();
        x = res;
        lambda2 = x.transpose() * A * x;
        i++;
    }
    while (i < N && abs(lambda1-lambda2) >= e);
    double lambda = x.transpose() * A * x;

    return lambda;
}

float ECM(VectorXd x, VectorXd b){
    float res = 0.0;
    for(int i  = 0; i < x.size(); i++){
        res += (x(i) - b(i)) * (x(i) - b(i));
    }

    res = res / x.size();

    return res;
}

MatrixXd toDouble(MatrixXi &entry) {
	int size = entry.cols();
	MatrixXd A(size, size);
	for (int i = 0; i < size; i++ ) {
		for (int j = 0; j < size; j++) {
			A(i, j) = double(entry(i, j));
		}
	}
	return A;
}