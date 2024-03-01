#include <iostream>
#include <cstdio>
#include "Eigen/Dense"
#include "JACOBI.h"
#include "GAUSS-SEIDEL.h"
#include "JACOBI.cpp"
#include "GAUSS-SEIDEL.cpp"
#include "DLU.h"
#include "DLU.cpp"
#include "LU.h"
#include "LU.cpp"
#include "Auxiliares.h"
#include "Auxiliares.cpp"
#include <chrono>
#include "Experimentacion.h"
#include "Archivos.h"
#include "Archivos.cpp"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

using Eigen::MatrixXd;
using Eigen::VectorXd;

void Experimentacion_XMult_J_Mat(MatrixXd &A, VectorXd &x, float epsilon, vector<int> iters){
    int cant_it = iters.size();

    VectorXd b = A * x;
    cout << "b = " << b << endl;
    int size = x.size();
    VectorXd x_init(size);
    /* for (int i = 0; i < size; i++) {
        x_init[i] = rand();   // Lo llenamos de coeficientes aleatorios
    } */

    vector<float> errs_x_2(cant_it, 0.0);
    vector<float> errs_x_10(cant_it, 0.0);
    vector<float> errs_x_100(cant_it, 0.0);
    vector<float> errs_x_1000(cant_it, 0.0);

    for(int i = 0; i < iters.size(); i++){
        VectorXd Ax;

        VectorXd res_x_2 = jacobiMethodXMult(A, b, x, 2, iters[i], epsilon);
        Ax = A * res_x_2;
        cout << "res_x_2 para " + to_string(iters[i]) + " iteraciones: " << res_x_2 << endl;
        cout << "Ax: " << Ax << endl;
        errs_x_2[i] = ECM(Ax, b);

        VectorXd res_x_10 = jacobiMethodXMult(A, b, x, 10, iters[i], epsilon);
        Ax = A * res_x_10;
        cout << "res_x_10 para " + to_string(iters[i]) + " iteraciones: " << res_x_10 << endl;
        cout << "Ax: " << Ax << endl;
        errs_x_10[i] = ECM(Ax, b);

        VectorXd res_x_100 = jacobiMethodXMult(A, b, x, 100, iters[i], epsilon);
        Ax = A * res_x_100;
        cout << "res_x_100 para " + to_string(iters[i]) + " iteraciones: " << res_x_100 << endl;
        cout << "Ax: " << Ax << endl;
        errs_x_100[i] = ECM(Ax, b);

        VectorXd res_x_1000 = jacobiMethodXMult(A, b, x, 1000, iters[i], epsilon);
        Ax = A * res_x_1000;
        cout << "res_x_1000 para " + to_string(iters[i]) + " iteraciones: " << res_x_1000 << endl;
        cout << "Ax: " << Ax << endl;
        errs_x_1000[i] = ECM(Ax, b);
    }

    generarArchivo_Error_XMult("Error_XMult_J_Mat.txt", errs_x_2, errs_x_10, errs_x_100, errs_x_1000);
}

void Experimentacion_XMult_J_Sum(MatrixXd &A, VectorXd &x, float epsilon, vector<int> iters){
    int cant_it = iters.size();

    VectorXd b = A * x;
    int size = x.size();
    VectorXd x_init(size);
    /* for (int i = 0; i < size; i++) {
        x_init[i] = rand();   // Lo llenamos de coeficientes aleatorios
    } */

    vector<float> errs_x_2(cant_it, 0.0);
    vector<float> errs_x_10(cant_it, 0.0);
    vector<float> errs_x_100(cant_it, 0.0);
    vector<float> errs_x_1000(cant_it, 0.0);

    for(int i = 0; i < iters.size(); i++){
        VectorXd Ax;

        VectorXd res_x_2 = jacobiSumatoriaXMult(A, b, x, 2, iters[i], epsilon);
        Ax = A * res_x_2;
        errs_x_2[i] = ECM(Ax, b);

        VectorXd res_x_10 = jacobiSumatoriaXMult(A, b, x, 10, iters[i], epsilon);
        Ax = A * res_x_10;
        errs_x_10[i] = ECM(Ax, b);

        VectorXd res_x_100 = jacobiSumatoriaXMult(A, b, x, 100, iters[i], epsilon);
        Ax = A * res_x_100;
        errs_x_100[i] = ECM(Ax, b);

        VectorXd res_x_1000 = jacobiSumatoriaXMult(A, b, x, 1000, iters[i], epsilon);
        Ax = A * res_x_1000;
        errs_x_1000[i] = ECM(Ax, b);
    }

    generarArchivo_Error_XMult("Error_XMult_J_Sum.txt", errs_x_2, errs_x_10, errs_x_100, errs_x_1000);
}

void Experimentacion_XMult_GS_Mat(MatrixXd &A, VectorXd &x, float epsilon, vector<int> iters){
    int cant_it = iters.size();

    VectorXd b = A * x;
    int size = x.size();
    VectorXd x_init(size);
     /* for (int i = 0; i < size; i++) {
        x_init[i] = rand();   // Lo llenamos de coeficientes aleatorios
    } */

    vector<float> errs_x_2(cant_it, 0.0);
    vector<float> errs_x_10(cant_it, 0.0);
    vector<float> errs_x_100(cant_it, 0.0);
    vector<float> errs_x_1000(cant_it, 0.0);

    for(int i = 0; i < iters.size(); i++){
        VectorXd Ax;

        VectorXd res_x_2 = gaussSeidelMethodXMult(A, b, x, 2, iters[i], epsilon);
        Ax = A * res_x_2;
        errs_x_2[i] = ECM(Ax, b);

        VectorXd res_x_10 = gaussSeidelMethodXMult(A, b, x, 10, iters[i], epsilon);
        Ax = A * res_x_10;
        errs_x_10[i] = ECM(Ax, b);

        VectorXd res_x_100 = gaussSeidelMethodXMult(A, b, x, 100, iters[i], epsilon);
        Ax = A * res_x_100;
        errs_x_100[i] = ECM(Ax, b);

        VectorXd res_x_1000 = gaussSeidelMethodXMult(A, b, x, 1000, iters[i], epsilon);
        Ax = A * res_x_1000;
        errs_x_1000[i] = ECM(Ax, b);
    }

    generarArchivo_Error_XMult("Error_XMult_GS_Mat.txt", errs_x_2, errs_x_10, errs_x_100, errs_x_1000);
}

void Experimentacion_XMult_GS_Sum(MatrixXd &A, VectorXd &x, float epsilon, vector<int> iters){
    int cant_it = iters.size();

    VectorXd b = A * x;
    int size = x.size();
    VectorXd x_init(size);
    /* for (int i = 0; i < size; i++) {
        x_init[i] = rand();   // Lo llenamos de coeficientes aleatorios
    } */

    vector<float> errs_x_2(cant_it, 0.0);
    vector<float> errs_x_10(cant_it, 0.0);
    vector<float> errs_x_100(cant_it, 0.0);
    vector<float> errs_x_1000(cant_it, 0.0);

    for(int i = 0; i < iters.size(); i++){
        VectorXd Ax;

        VectorXd res_x_2 = gaussSeidelSumatoriaXMult(A, b, x, 2, iters[i], epsilon);
        Ax = A * res_x_2;
        errs_x_2[i] = ECM(Ax, b);

        VectorXd res_x_10 = gaussSeidelSumatoriaXMult(A, b, x, 10, iters[i], epsilon);
        Ax = A * res_x_10;
        errs_x_10[i] = ECM(Ax, b);

        VectorXd res_x_100 = gaussSeidelSumatoriaXMult(A, b, x, 100, iters[i], epsilon);
        Ax = A * res_x_100;
        errs_x_100[i] = ECM(Ax, b);

        VectorXd res_x_1000 = gaussSeidelSumatoriaXMult(A, b, x, 1000, iters[i], epsilon);
        Ax = A * res_x_1000;
        errs_x_1000[i] = ECM(Ax, b);
    }

    generarArchivo_Error_XMult("Error_XMult_GS_Sum.txt", errs_x_2, errs_x_10, errs_x_100, errs_x_1000);
}

void correrExperimentacion() {
    float e = 0.00000000000000001;

    vector<int> iteraciones{10, 50, 100, 300, 500, 1000, 2500};
    int cant_it = iteraciones.size();

    vector<int> iteraciones_xMult{10, 13, 15, 17, 20, 30, 40, 50, 100};

    vector<int> sizes{2, 4, 8, 16, 32, 64, 128, 256, 512};
    int cant_sizes = sizes.size();

    vector<float> errs_J_Mat(cant_it, 0.0);
    vector<float> errs_J_Sum(cant_it, 0.0);
    vector<float> errs_GS_Mat(cant_it, 0.0);
    vector<float> errs_GS_Sum(cant_it, 0.0);
    vector<float> errs_LU(cant_it, 0.0);

    vector<float> times_J_Mat(cant_sizes, 0.0);
    vector<float> times_J_Sum(cant_sizes, 0.0);
    vector<float> times_GS_Mat(cant_sizes, 0.0);
    vector<float> times_GS_Sum(cant_sizes, 0.0);
    vector<float> times_LU(cant_sizes, 0.0);

    vector<float> times_J_Mat_segun_cant_it(cant_it, 0.0);
    vector<float> times_J_Sum_segun_cant_it(cant_it, 0.0);
    vector<float> times_GS_Mat_segun_cant_it(cant_it, 0.0);
    vector<float> times_GS_Sum_segun_cant_it(cant_it, 0.0);
    vector<float> times_LU_segun_cant_it(cant_it, 0.0);

    vector<float> time_x_it_J_Mat(10, 0.0);
    vector<float> time_x_it_J_Sum(10, 0.0);
    vector<float> time_x_it_GS_Mat(10, 0.0);
    vector<float> time_x_it_GS_Sum(10, 0.0);
    vector<float> time_x_it_LU(10, 0.0);

    for (int i = 0; i < cant_sizes; i++) {
        int k = sizes[i];
        cout << "MATRICES TAMAÑO " << to_string(k) << endl;
        for (int s = 0; s < 10; s++) {
            cout << "ITERACIÓN N° " << to_string(s) << endl;
            MatrixXi A_temp = MatrixXi::Random(k, k); // GENERAMOS MATRIZ RANDOM DE ENTEROS
            VectorXd b = VectorXd::Random(k);

            MatrixXd A(k, k);
            for (int t = 0; t < k; t++) {
                for (int h = 0; h < k; h++) {
                    double elem = A_temp(t, h); // PASAMOS A DOUBLES
                    A(t, h) = elem;
                }
            }

            for (int t = 0; t < k; t++) {       // GENERAMOS MATRIZ EDD
                double suma = A.row(t).cwiseAbs().sum();
                suma = suma + A.col(t).cwiseAbs().sum();
                A(t, t) = suma + 1;
            }

            auto start = steady_clock::now();
            VectorXd resJ_Mat = jacobiMethodHastaConvergencia(A, b, e);
            auto end = steady_clock::now();
            times_J_Mat[i] += float(duration_cast<milliseconds>(end - start).count()) / 10;
            time_x_it_J_Mat[s] = float(duration_cast<milliseconds>(end - start).count());

            start = steady_clock::now();
            VectorXd resJ_Sum = jacobiSumatoriaHastaConvergencia(A, b, e);
            end = steady_clock::now();
            times_J_Sum[i] += float(duration_cast<milliseconds>(end - start).count()) / 10;
            time_x_it_J_Sum[s] = float(duration_cast<milliseconds>(end - start).count());

            start = steady_clock::now();
            VectorXd resGS_Mat = gaussSeidelMethodHastaConvergencia(A, b, e);
            end = steady_clock::now();
            times_GS_Mat[i] += float(duration_cast<milliseconds>(end - start).count()) / 10;
            time_x_it_GS_Mat[s] = float(duration_cast<milliseconds>(end - start).count());

            start = steady_clock::now();
            VectorXd resGS_Sum = gaussSeidelSumatoriaHastaConvergencia(A, b, e);
            end = steady_clock::now();
            times_GS_Sum[i] += float(duration_cast<milliseconds>(end - start).count()) / 10;
            time_x_it_GS_Sum[s] = float(duration_cast<milliseconds>(end - start).count());

            start = steady_clock::now();
            VectorXd res_LU = LUMethod(A, b);
            end = steady_clock::now();
            times_LU[i] += float(duration_cast<milliseconds>(end - start).count()) / 10;
            time_x_it_LU[s] = float(duration_cast<milliseconds>(end - start).count());

            for (int j = 0; j < iteraciones.size(); j++) {
                int N = iteraciones[j];
                VectorXd Ax(k);

                start = steady_clock::now();
                resJ_Mat = jacobiMethod(A, b, N, e);
                end = steady_clock::now();
                times_J_Mat_segun_cant_it[j] += float(duration_cast<milliseconds>(end - start).count()) / 10;
                Ax = A * resJ_Mat;
                errs_J_Mat[j] += ECM(Ax, b) / 10;
                cout << "Listo Jacobi Matricial para " << to_string(N) << " iteraciones" << endl;

                start = steady_clock::now();
                resJ_Sum = jacobiSumatoria(A, b, N, e);
                end = steady_clock::now();
                times_J_Sum_segun_cant_it[j] += float(duration_cast<milliseconds>(end - start).count()) / 10;
                Ax = A * resJ_Sum;
                errs_J_Sum[j] += ECM(Ax, b) / 10;
                cout << "Listo Jacobi Sumatoria para " << to_string(N) << " iteraciones" << endl;

                start = steady_clock::now();
                resGS_Mat = gaussSeidelMethod(A, b, N, e);
                end = steady_clock::now();
                times_GS_Mat_segun_cant_it[j] += float(duration_cast<milliseconds>(end - start).count()) / 10;
                Ax = A * resGS_Mat;
                errs_GS_Mat[j] += ECM(Ax, b) / 10;
                cout << "Listo Gauss-Seidel Matricial para " << to_string(N) << " iteraciones" << endl;

                start = steady_clock::now();
                resGS_Sum = gaussSeidelSumatoria(A, b, N, e);
                end = steady_clock::now();
                times_GS_Sum_segun_cant_it[j] += float(duration_cast<milliseconds>(end - start).count()) / 10;
                Ax = A * resGS_Sum;
                errs_GS_Sum[j] += ECM(Ax, b) / 10;
                cout << "Listo Gauss-Seidel Sumatoria para " << to_string(N) << " iteraciones" << endl;

                start = steady_clock::now();
                res_LU = LUMethod(A, b);
                end = steady_clock::now();
                times_LU_segun_cant_it[j] += float(duration_cast<milliseconds>(end - start).count()) / 10;
                Ax = A * res_LU;
                errs_LU[j] += ECM(Ax, b) / 10;
                cout << "Listo LU" << endl;
            }
        }

        string filename_err = "ErroresTamaño" + to_string(sizes[i]) + ".txt";
        string filename_time = "TimesTamaño" + to_string(sizes[i]) + ".txt";
        string filename_time_x_it = "TiempoPorItTamaño" + to_string(sizes[i]) + ".txt";

        generarArchivo_Error_Tiempo(filename_err, errs_J_Mat, errs_J_Sum, errs_GS_Mat, errs_GS_Sum, errs_LU);
        generarArchivo_Error_Tiempo(filename_time, times_J_Mat_segun_cant_it, times_J_Sum_segun_cant_it, times_GS_Mat_segun_cant_it, times_GS_Sum_segun_cant_it, times_LU_segun_cant_it);
        generarArchivo_Error_Tiempo(filename_time_x_it, time_x_it_J_Mat, time_x_it_J_Sum, time_x_it_GS_Mat, time_x_it_GS_Sum, time_x_it_LU);

    }

    generarArchivo_Error_Tiempo("TimesMetodosXTamaño.txt", times_J_Mat, times_J_Sum, times_GS_Mat, times_GS_Sum, times_LU);

    // Ahora vamos a generar resultados para después poder armar un gráfico que nos diga cómo cambia el error de los
    // métodos iterativos según qué tan lejos está el x inicial de la solución del sistema.
    // primero nos armamos una matriz edd de 512 x 512 y un vector random x

    MatrixXd A_5 = MatrixXd::Random(5, 5); // GENERAMOS MATRIZ RANDOM DE ENTEROS
    VectorXd x(5);
    for (int i = 0; i < 5; i++) {
        x[i] = rand();   // Lo llenamos de coeficientes aleatorios
    }

    for (int t = 0; t < 5; t++) {       // GENERAMOS MATRIZ EDD
        double suma = A_5.row(t).cwiseAbs().sum();
        suma = suma + A_5.col(t).cwiseAbs().sum();
        A_5(t, t) = suma + 1;
    }

    cout << "MATRIZ DE XMULT:" << endl;
    cout << A_5 << endl;


    Experimentacion_XMult_J_Mat(A_5, x, e, iteraciones_xMult);
    Experimentacion_XMult_J_Sum(A_5, x, e, iteraciones_xMult);
    Experimentacion_XMult_GS_Mat(A_5, x, e, iteraciones_xMult);
    Experimentacion_XMult_GS_Sum(A_5, x, e, iteraciones_xMult);

}