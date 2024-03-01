#include <iostream>
#include <fstream>
#include <cstdio>
#include <eigen3/Eigen/Dense>
#include <tuple>
#include <vector>
#include <chrono>


using namespace Eigen;
using namespace std;
using namespace std::chrono;

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;
 
MatrixXd deflacion(MatrixXd A, VectorXd v, float x) {
    // x es el lambda
    // v es el autovector asociado a lambda
    MatrixXd B = v * v.transpose();
    B = A - x * B;
    return B;
}

tuple<MatrixXd, VectorXd, double> metodoPotencia(MatrixXd A, int N, double e){
    // METODO POTENCIA
    // Suponemos que hay un autovalor maximo en modulo a todo el resto
    // Empezamos con un vector al azar NO NULO x
    // x se escribe como combinacion lineal de los autovectores de A
    // (Porque conforman una base)
    // Ax_0 = sum bj * A * vj = sum bj * lambda j * vj = x_1
    // A^k  * x_0 = b1 * lambda1^k * v1 + (sum bj * A * vj)
    // = lambda1^k * (b1 * v1 + sum (lambdaj/lambda1)^k * vj) = x_k
    // Normalizamos q_k = x_k / ||x_k|| converge al autovector v1
    // Ademas lambdak  = (q_k)t * A * q_k converge a lambda 1
    // TODO JUNTO 
    // Potencia para A, desinflamos, potencia para B ? etc etc creo que si
    int size = A.cols();
    // generamos al azar un vector no nulo de tamano size
    VectorXd x(size);
    x[0] = 1; // para asegurarnos de que no es nulo
    for (int i = 1; i < size; ++i ) {
        x[i] = rand();
    }
    // = lambda1^k * (b1 * v1 + sum (lambdaj/lambda1)^k * vj) = x_k
    // Ax_0 = sum bj * A * vj = sum bj * lambda j * vj = x_1
    // queremos multiplicar por A k veces.
    VectorXd res (size);
    // for (int i = 0; i < N; i++) {
        // res = A * x;
        // res.normalize();
        // x = res;
    // }

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
    // lambdak  = (q_k)t * A * q_k converge a lambda 1
    double lambda = x.transpose() * A * x;
    tuple <MatrixXd, VectorXd, double> result; 
    get<0>(result) = A;
    get<1>(result) = x;
    get<2>(result) = lambda;
    // res = make_tuple<A, x, res>;
    return result;
}
// Para leer matriz de un archivo
// MatrixXd A(nrows, ncols);
//   for (int i = 0; i < nrows; i++) {
//       for (int j = 0; j < ncols; j++) {
//            fin >> A(i, j);
//       }
//    }

int metodoPotenciaDeflacion(string matriz, int N, string values_filename, string vectors_filename, string times_filename){
    // double e = 0.000001;
    double e = 0.000000000001;
    ifstream fin(matriz);
    if (!fin.is_open()) {
        cerr << "Error: no se puede abrir el archivo con los datos de entrada " << matriz << endl;
        return 1;
    }
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
    vector<VectorXd> vecs;
    vector<double> vals;
    vector<float> tiempo_x_autovalor;

    unsigned t0 = clock(); // esto para calcular el tiempo total del método (no se si lo usamos igual)
    for (int i = 0; i < A.cols(); i++){
        auto start = steady_clock::now(); // para calcular el tiempo de convergencia para autovalor, yendo desde el de mayor módulo hasta el de menor
        tuple<MatrixXd, VectorXd, double> results = metodoPotencia(A, N, e);
        vecs.push_back(get<1>(results));
        vals.push_back(get<2>(results));
        A = deflacion(get<0>(results), get<1>(results), get<2>(results));
        auto end = steady_clock::now();

        float tiempo_av = float(duration_cast<nanoseconds>(end - start).count()); // Lo vamos a medir en cantidad de clocks porque en segundos da siempre 0 porque es muy rápido para los ejemplos del 1 b)
        tiempo_x_autovalor.push_back(tiempo_av);
    }
    unsigned t1 = clock();

    double tiempo_total = 0.0;
    tiempo_total = (t1 - t0) / CLOCKS_PER_SEC;

    // Imprimo autovalores en archivo
    ofstream fout(values_filename);
    if (!fout.is_open()) {
        cerr << "Error: no se puede abrir el archivo de salida " << values_filename << endl;
        return 1;
    }

    for (int i = 0; i < (vals).size(); i++) {
        fout << vals[i] << endl;
    }
    
    fout.close();

    // Imprimo autovectores en segundo archivo
    ofstream fout2(vectors_filename);
    if (!fout2.is_open()) {
        cerr << "Error: no se puede abrir el archivo de salida " << vectors_filename << endl;
        return 1;
    }

    for (int i = 0; i < vecs.size(); i++) {
        for (int j = 0; j < vecs.size(); j++) {
            fout2 << vecs[j][i] << " ";
            if (j == vecs.size() - 1) {
                fout2 << endl;
            }
        }

    }
    fout2.close();

    // Imprimo el tiempo de convergencia para cada autovalor en un tercer archivo
    ofstream fout3(times_filename);

    for (int i = 0; i < (tiempo_x_autovalor).size(); i++) {
        if(i == tiempo_x_autovalor.size() - 1){
            fout3 << tiempo_x_autovalor[i] << endl;
        } else {
            fout3 << tiempo_x_autovalor[i] << " ";
        }
    }

    fout3.close();

    return 0;
    }
