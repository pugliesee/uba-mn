#include "Eigen/Dense"
#include <vector>
#include <iostream>
#include <tuple>
#include <cstdio>
#include "JACOBI.h"
#include "Auxiliares.h"


using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;


tuple<MatrixXd, MatrixXd> matrizIteracionJacobi(tuple<MatrixXd, MatrixXd, MatrixXd> &DLU) {
	MatrixXd D = get<0>(DLU);
	MatrixXd L = get<1>(DLU);
	MatrixXd U = get<2>(DLU);

	if (!sinCerosDiagonal(D)) {
		printf("No se puede aplicar el método de Jacobi, D no es inversible");
        tuple<MatrixXd, MatrixXd> res;
        get<0>(res) = D;
        get<1>(res) = L + U;
		return res;
	}

	D = D.inverse();
	MatrixXd LU = L + U;
	MatrixXd T_jacobi = D * LU;

	tuple<MatrixXd, MatrixXd> res;
	get<0>(res) = T_jacobi;
	get<1>(res) = D;

	return res;
}


VectorXd jacobiMethod(MatrixXd &A, VectorXd &b, int N, float epsilon) {
	// x_(k+1) = T * x_(k) + D.inverse() * b
	// PRIMERO, NECESITAMOS REESCRIBIR A COMO A = D - L - U

	tuple<MatrixXd, MatrixXd, MatrixXd> DLU = eigenDLU(A);
	tuple<MatrixXd, MatrixXd> TyD = matrizIteracionJacobi(DLU);

	MatrixXd T = get<0>(TyD);
	MatrixXd D_inverse = get<1>(TyD);

    double radio_espectral = metodoPotencia(T, 5000, epsilon);
    cout << "Radio espectral: " << radio_espectral << endl;
    if (abs(radio_espectral) >= 1) {
        cout << "El método diverge" << endl;
        assert(false);
    }

	int size = T.cols(); // T matriz cuadrada --> T.cols == T.rows()
	VectorXd x(size);
	for (int i = 0; i < size; i++) {
	   x[i] = rand();   // Lo llenamos de coeficientes aleatorios
	}

	VectorXd x_cmp(size); // Vector para comparar la diferencia entre x_(k+1) y x_(k)
	int k = 0;
	do {

		x_cmp = x;
		x = T * x;
        x += D_inverse * b;     // Separé las cuentas porque sino el programa hacía cualquier cosa
		k++;
			// Comparo las normas y si el error es muy chico, doy por terminada la iteración
	} while (k < N && abs((x - x_cmp).norm()) >= epsilon);

	return x;
}

VectorXd jacobiMethodXMult(MatrixXd &A, VectorXd &b, VectorXd &x_init, int d, int N, float epsilon) {
    // x_(k+1) = T * x_(k) + D.inverse() * b
    // PRIMERO, NECESITAMOS REESCRIBIR A COMO A = D - L - U

    tuple<MatrixXd, MatrixXd, MatrixXd> DLU = eigenDLU(A);
    tuple<MatrixXd, MatrixXd> TyD = matrizIteracionJacobi(DLU);

    MatrixXd T = get<0>(TyD);
    MatrixXd D_inverse = get<1>(TyD);

    double radio_espectral = metodoPotencia(T, 5000, epsilon);
    cout << "Radio espectral: " << radio_espectral << endl;
    if (abs(radio_espectral) >= 1) {
        cout << "El método diverge" << endl;
        assert(false);
    }

    int size = T.cols(); // T matriz cuadrada --> T.cols == T.rows()
    VectorXd x = x_init * d;
    x(0) += 1;

    VectorXd x_cmp(size); // Vector para comparar la diferencia entre x_(k+1) y x_(k)
    int k = 0;
    do {

        x_cmp = x;
        x = T * x;
        x += D_inverse * b;     // Separé las cuentas porque sino el programa hacía cualquier cosa
        k++;
        // Comparo las normas y si el error es muy chico, doy por terminada la iteración
    } while (k < N && abs((x - x_cmp).norm()) >= epsilon);

    return x;
}

VectorXd jacobiMethodHastaConvergencia(MatrixXd &A, VectorXd &b, float epsilon) {
    // x_(k+1) = T * x_(k) + D.inverse() * b
    // PRIMERO, NECESITAMOS REESCRIBIR A COMO A = D - L - U

    tuple<MatrixXd, MatrixXd, MatrixXd> DLU = eigenDLU(A);
    tuple<MatrixXd, MatrixXd> TyD = matrizIteracionJacobi(DLU);

    MatrixXd T = get<0>(TyD);
    MatrixXd D_inverse = get<1>(TyD);

    double radio_espectral = metodoPotencia(T, 5000, epsilon);
    cout << "Radio espectral: " << radio_espectral << endl;
    if (abs(radio_espectral) >= 1) {
        cout << "El método diverge" << endl;
        assert(false);
    }

    int size = T.cols(); // T matriz cuadrada --> T.cols == T.rows()
    VectorXd x(size);
    for (int i = 0; i < size; i++) {
        x[i] = rand();   // Lo llenamos de coeficientes aleatorios
    }

    VectorXd x_cmp(size); // Vector para comparar la diferencia entre x_(k+1) y x_(k)
    do {
        x_cmp = x;
        x = T * x;
        x += D_inverse * b;     // Separé las cuentas porque sino el programa hacía cualquier cosa
        // Comparo las normas y si el error es muy chico, doy por terminada la iteración
    } while (abs((x - x_cmp).norm()) >= epsilon);

    return x;
}

VectorXd jacobiSumatoria(MatrixXd &A, VectorXd &b, int N, float epsilon) {
    int size = A.cols();    // A matriz cuadrada --> A.cols() = A.rows()
    int k = 0;
    int err_count = 0;
    float err = 0.0;

    VectorXd x(size);
    for (int i = 0; i < size; i++) {
        x[i] = rand();   // Lo llenamos de coeficientes aleatorios
    }
    VectorXd x_cmp = x;
    do {
        x_cmp = x;
        x = b;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (j != i) {
                    x(i) -= A(i, j) * x_cmp(j);
                }
            }
            x(i) = x(i) / A(i,i);
        }
        k++;

        if (ECM(x, x_cmp) > err){
            err_count++;
        } else {
            err_count = 0;
        }

        if (err_count == 10) {
            cout << "El método diverge" << endl;
            assert(false);
        }

        err = ECM(x, x_cmp);

    } while (k < N && (x - x_cmp).norm() >= epsilon);

    return x;
}

VectorXd jacobiSumatoriaXMult(MatrixXd &A, VectorXd &b, VectorXd &x_init, int d, int N, float epsilon) {
    int size = A.cols();    // A matriz cuadrada --> A.cols() = A.rows()
    int k = 0;
    int err_count = 0;
    float err = 0.0;

    VectorXd x = x_init * d;
    x(0) += 1;
    VectorXd x_cmp = x;
    do {
        x_cmp = x;
        x = b;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (j != i) {
                    x(i) -= A(i, j) * x_cmp(j);
                }
            }
            x(i) = x(i) / A(i,i);
        }
        k++;

        if (ECM(x, x_cmp) > err){
            err_count++;
        } else {
            err_count = 0;
        }

        if (err_count == 10) {
            cout << "El método diverge" << endl;
            assert(false);
        }

        err = ECM(x, x_cmp);

    } while (k < N && (x - x_cmp).norm() >= epsilon);

    return x;
}

VectorXd jacobiSumatoriaHastaConvergencia(MatrixXd &A, VectorXd &b, float epsilon) {
    int size = A.cols();    // A matriz cuadrada --> A.cols() = A.rows()
    int err_count = 0;
    float err = 0.0;

    VectorXd x(size);
    for (int i = 0; i < size; i++) {
        x[i] = rand();   // Lo llenamos de coeficientes aleatorios
    }
    VectorXd x_cmp = x;
    do {
        x_cmp = x;
        x = b;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (j != i) {
                    x(i) -= A(i, j) * x_cmp(j);
                }
            }
            x(i) = x(i) / A(i,i);
        }

        if (ECM(x, x_cmp) > err){
            err_count++;
        } else {
            err_count = 0;
        }

        if (err_count == 10) {
            cout << "El método diverge" << endl;
            assert(false);
        }

        err = ECM(x, x_cmp);

    } while ((x - x_cmp).norm() >= epsilon);

    return x;
}