#include "GAUSS-SEIDEL.h"
#include <vector>
#include <tuple>
#include <iostream>
#include "Auxiliares.h"

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

tuple<MatrixXd, MatrixXd> matrizIteracionGS(tuple<MatrixXd, MatrixXd, MatrixXd> &DLU) {
    MatrixXd D = get<0>(DLU);
    MatrixXd L = get<1>(DLU);
    MatrixXd U = get<2>(DLU);


    int size = D.cols(); // Matrices cuadradas --> D.cols() == D.rows()
    for (int i = 0; i < size; i++) {
        if (D(i,i) == 0) {
            cout << "No se puede aplicar el método de Gauss-Seidel en esta matriz." << endl;
            assert(false);   // TODO: VER QUE PASA SI LLEGO A ESTE ASSERT
        }
    }

    MatrixXd DL = D - L;
    DL = DL.inverse();
    MatrixXd T_gs = DL * U;

    tuple<MatrixXd, MatrixXd> res;
    get<0>(res) = T_gs;
    get<1>(res) = DL;

    return res;
}


VectorXd gaussSeidelMethod(MatrixXd &A, VectorXd &b, int N, float epsilon) {
    tuple<MatrixXd, MatrixXd, MatrixXd> dlu = eigenDLU(A);
    tuple<MatrixXd, MatrixXd> TyDL = matrizIteracionGS(dlu);
    // Descompongo la tupla
    MatrixXd T = get<0>(TyDL);
    MatrixXd DL = get<1>(TyDL);

    double radio_espectral = metodoPotencia(T, 5000, epsilon);
    cout << "Radio espectral: " << radio_espectral << endl;
    if (abs(radio_espectral) >= 1) {
        cout << "El método diverge" << endl;
        assert(false);
    }

    // Ahora tengo que crear el vector resultante y asegurarme de que no sea nulo
    int size = T.cols(); // Matrices cuadradas --> T.cols() == T.rows()
    VectorXd x(size);
    for (int i = 0; i < size; i++) {
        x[i] = rand();
    }

    VectorXd x_cmp = x;
    int k = 0;
    do {

        x_cmp = x;
        x = T*x;
        x += DL*b;
        k++;

    } while (k < N && abs((x - x_cmp).norm()) >= epsilon);

    return x;

}

VectorXd gaussSeidelMethodXMult(MatrixXd &A, VectorXd &b, VectorXd &x_init, float d, int N, float epsilon) {
    tuple<MatrixXd, MatrixXd, MatrixXd> dlu = eigenDLU(A);
    tuple<MatrixXd, MatrixXd> TyDL = matrizIteracionGS(dlu);
    // Descompongo la tupla
    MatrixXd T = get<0>(TyDL);
    MatrixXd DL = get<1>(TyDL);

    double radio_espectral = metodoPotencia(T, 5000, epsilon);
    cout << "Radio espectral: " << radio_espectral << endl;
    if (abs(radio_espectral) >= 1) {
        cout << "El método diverge" << endl;
        assert(false);
    }

    // Ahora tengo que crear el vector resultante y asegurarme de que no sea nulo
    int size = T.cols(); // Matrices cuadradas --> T.cols() == T.rows()
    VectorXd x = x_init * d;
    x(0) += 1;
    VectorXd x_cmp = x;
    int k = 0;
    do {

        x_cmp = x;
        x = T*x;
        x += DL*b;
        k++;

    } while (k < N && abs((x - x_cmp).norm()) >= epsilon);

    return x;

}

VectorXd gaussSeidelMethodHastaConvergencia(MatrixXd &A, VectorXd &b, float epsilon) {
    tuple<MatrixXd, MatrixXd, MatrixXd> dlu = eigenDLU(A);
    tuple<MatrixXd, MatrixXd> TyDL = matrizIteracionGS(dlu);
    // Descompongo la tupla
    MatrixXd T = get<0>(TyDL);
    MatrixXd DL = get<1>(TyDL);

    double radio_espectral = metodoPotencia(T, 5000, epsilon);
    cout << "Radio espectral: " << radio_espectral << endl;
    if (abs(radio_espectral) >= 1) {
        cout << "El método diverge" << endl;
        assert(false);
    }

    // Ahora tengo que crear el vector resultante y asegurarme de que no sea nulo
    int size = T.cols(); // Matrices cuadradas --> T.cols() == T.rows()
    VectorXd x(size);
    for (int i = 0; i < size; i++) {
        x[i] = rand();
    }

    VectorXd x_cmp = x;
    do {
        x_cmp = x;
        x = T*x;
        x += DL*b;
    } while (abs((x - x_cmp).norm()) >= epsilon);

    return x;
}


VectorXd gaussSeidelSumatoria(MatrixXd &A, VectorXd &b, int N, float epsilon) {
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
                    if (j < i) {
                        x(i) -= A(i, j) * x(j);
                    } else {
                        x(i) -= A(i, j) * x_cmp(j);
                    }
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

VectorXd gaussSeidelSumatoriaXMult(MatrixXd &A, VectorXd &b, VectorXd &x_init, float d, int N, float epsilon) {
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
                    if (j < i) {
                        x(i) -= A(i, j) * x(j);
                    } else {
                        x(i) -= A(i, j) * x_cmp(j);
                    }
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

VectorXd gaussSeidelSumatoriaHastaConvergencia(MatrixXd &A, VectorXd &b, float epsilon) {
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
                    if (j < i) {
                        x(i) -= A(i, j) * x(j);
                    } else {
                        x(i) -= A(i, j) * x_cmp(j);
                    }
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
