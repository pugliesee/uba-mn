#include <iostream>
#include <fstream>
#include <cstdio>
#include <eigen3/Eigen/Dense>
#include <tuple>
#include "metodoPotDef.h"
#include <vector>
#include <string>
#include <map>

int main () {
    /* RECORDAR CORRER EL MAIN.PY ANTES DE EJECUTAR ESTE ARCHIVO
    TESTEOS EJERCICIO 1
    Los archivos de prueba fueron incluídos en caso de querer "jugar" un poco con el programa
    y para generar resultados para el ejercicio 1 b)
    Los tests tienen las siguientes características:
    TEST 1: Identidad
    TEST 2: Triangular
    TEST 3: Simétrica
    TEST 4: Random
    Finalmente en el informe el gráfico que involucra a los tests 1, 2, 3, y 4 no fue incluido
    TESTS 5 A 8: Diagonales con autovalores de distinta magnitud y con diferentes distancias entre ellos
    Todos los tests son con matrices de 10 x 10 */
    /* map<int, vector<float>> tiempos_tests;
    // en tiempos_tests[9] guardo los promedio entre los tests 1, 2, 3, y 4 (promedio de los promedios)
    for(int i = 0; i < 10; i++){
        tiempos_tests[i] = vector<float>(10, 0.0); // Arreglo para tiempos promedio de cada test
    }

    for(int i = 1; i < 9; i++){
        // voy a correr 10 veces cada test
        for(int k = 0; k < 10; k++){
            string nro_test = to_string(i);
            // corremos el método para el test i
            metodoPotenciaDeflacion("Test" + nro_test + "_1B.txt", 500, "AutovaloresTest" + nro_test + ".txt", "AutovectoresTest" + nro_test + ".txt", "TiemposTest" + nro_test + ".txt");
            vector<float> tiempos_test(10, 0.0);
            // colocamos en tiempos_test los tiempos de convergencia para cada autovector para el test i
            ifstream fin("TiemposTest" + nro_test + ".txt");
            for(int j = 0; j < 10; j++){
                fin >> tiempos_test[j];
                (tiempos_tests[i])[j] += tiempos_test[j] / 10; // coloco el tiempo para el autovalor j en la posición adecuada y lo divido por 10 ya que nos vamos a quedar con el promedio entre las 10 corridas del test
            }
            fin.close();
        }
    }

    for(int s = 0; s < 10; s++){
        for(int t = 1; t < 5; t++){
            (tiempos_tests[9])[s] += (tiempos_tests[t])[s] / 4; // Me voy guardando el tiempo promedio entre los 4 primeros tests (promedio de los promedios)
        }
    }

    for(int t = 1; t < 10; t++){
        ofstream fout;
        if(t == 9){
            fout.open("TiemposPromedio.txt");
        } else {
            fout.open("TiemposTest" + to_string(t) + ".txt");
        }
        for(int i = 0; i < tiempos_tests[t].size(); i++){
            if(i == tiempos_tests[t].size() - 1){
                fout << (tiempos_tests[t])[i] << endl;
            } else {
                fout << (tiempos_tests[t])[i] << " ";
            }
        }
        fout.close();
    }

    // Para después poder armar el gráfico de error según el epsilon
    for(int t = 5; t < 9; t++){
        MatrixXd mat(10, 10);
        ifstream fin("Test" + to_string(t) + "_1B.txt");
        int nrows, ncols;
        fin >> nrows >> ncols;
        for(int i = 0; i < 10; i++){
            for(int j = 0; j < 10; j++){
                fin >> mat(i, j);
            }
        }
        fin.close();

        Eigen::EigenSolver<MatrixXd> eig(mat);
        VectorXcd vals = eig.eigenvalues();
        MatrixXcd vecs = eig.eigenvectors();

        ofstream fout("EigenvectorsTest" + to_string(t) + ".txt");
        for(int i = 0; i < 10; i++){
            for(int j = 0; j < 10; j++){
                if( j == 9){
                    fout << vecs(i, j).real() << endl;
                } else {
                    fout << vecs(i, j).real() << " ";
                }
            }
        }
        fout.close();
    }

    // Para gráfico de error según cantidad de iteraciones
    // Vamos a utilizar la matriz del test 5
    vector<int> iteraciones{10, 20, 30, 40, 50, 80, 100, 150, 200, 300, 500};

    for(int i = 0; i < iteraciones.size(); i++){
        metodoPotenciaDeflacion("Test5_1B.txt", iteraciones[i], "AutovaloresTest.txt", "AutovectoresTest5" + to_string(iteraciones[i]) + "Iters.txt", "TiemposTest5.txt");
    } */


    // EJECUTAR LA LINEA 116 PARA OBTENER LOS AUTOVALORES Y AUTOVECTORES USANDO PCA

    // En el segundo parámetro podemos elegir la cantidad de iteraciones que realiza el método para conseguir cada autovalor y autovector.
    // Nosotros utilizamos 300 a lo largo del proyecto.
    // El épsilon fue incluído directamente dentro del código. Usamos Epsilon =  1x10^(-12)
    // El tercer y cuarto parámetro fue creado por comodidad, para que los archivos de texto con los resultados lleven esos nombres.

    // metodoPotenciaDeflacion("PCATest.txt", 300, "autovalores_PCA.txt", "autovectores_PCA.txt", "timepos_PCA.txt");


    // EJECUTAR LA LINEA 124 PARA OBTENER LOS AUTOVALORES Y AUTOVECTORES USANDO 2DPCA

    //Los parámetros siguen funcionando igual, lo único que cambia es el nombre de los archivos de entrada y salida.
    // Otra vez, nosotros usamos 300 iteraciones para calcular dichos resultados.

    metodoPotenciaDeflacion("2DPCA.txt", 300, "autovalores_2DPCA.txt", "autovectores_2DPCA.txt", "tiempos_2DPCA.txt");

    // TESTEOS EJERCICIO 3 C

    // metodoPotenciaDeflacion("PCATest_reduced.txt", 300, "autovalores_PCA_reduced.txt", "autovectores_PCA_reduced.txt", "tiempos_PCA_reduced.txt");
    // metodoPotenciaDeflacion("PCATest_reduced_sinPrimero.txt", 300, "autovalores_PCA_reduced_sinPrimero.txt", "autovectores_PCA_reduced_sinPrimero.txt", "tiempos_PCA_reduced_sinPrimero.txt");
    // metodoPotenciaDeflacion("2DPCA_sinPrimero.txt", 300, "autovalores_2DPCA_sinPrimero.txt", "autovectores_2DPCA_sinPrimero.txt", "tiempos_2DPCA_sinPrimero.txt");
}