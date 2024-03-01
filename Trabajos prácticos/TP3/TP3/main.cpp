#include "Experimentacion.h"
#include "Experimentacion.cpp"

int main(){
    correrExperimentacion();

	float e = 0.00000000000000001;

	MatrixXi test1(5, 5);	// Probamos con una matriz EDD
	test1 << 10, 2, 2, 1, 0,
			3, 19, 3, 5, 6,
			1, 5, 90, 5, 48,
			0, 9, 21, 50, 4,
			1, 0, 8, 3, 71;

	MatrixXd A1 = toDouble(test1);	// VER FUNCIÓN AUXILIAR
	VectorXd x1 = VectorXd::Random(5);	// Generamos un x aleatorio

	VectorXd b1 = A1 * x1;	// Vector resultante b de tamaño 5

	// Y ahora queremos encontrar la solución con cada uno de los métodos
	cout << "Vamos con la matriz EDD" << endl;
	VectorXd resJMat1 = jacobiMethodHastaConvergencia(A1, b1, e);
	VectorXd resJSum1 = jacobiSumatoriaHastaConvergencia(A1, b1, e);
	VectorXd resGSMat1 = gaussSeidelMethodHastaConvergencia(A1, b1, e);
	VectorXd resGSSum1= gaussSeidelSumatoriaHastaConvergencia(A1, b1, e);
	VectorXd resLU1 = LUMethod(A1, b1);

    vector<float> errs_edd{ECM(A1*resJMat1, b1), ECM(A1*resJSum1, b1), ECM(A1*resGSMat1, b1), ECM(A1*resGSSum1, b1), ECM(A1*resLU1, b1)};


	MatrixXi test2(5, 5);	// MATRIZ IDENTIDAD ---> IGUAL A EDD. Va a converger para cualquier método.
	test2 << 1, 0, 0, 0, 0,
			 0, 1, 0, 0, 0,
			 0, 0, 1, 0, 0,
			 0, 0, 0, 1, 0,
             0, 0, 0, 0, 1;

    MatrixXd A2 = toDouble(test2);
    VectorXd x2 = VectorXd::Random(5);
    VectorXd b2 = A2 * x2;

    cout << "Vamos con la matriz Identidad" << endl;
    VectorXd resJMat2 = jacobiMethodHastaConvergencia(A2, b2, e);
    VectorXd resJSum2 = jacobiSumatoriaHastaConvergencia(A2, b2, e);
    VectorXd resGSMat2 = gaussSeidelMethodHastaConvergencia(A2, b2, e);
    VectorXd resGSSum2= gaussSeidelSumatoriaHastaConvergencia(A2, b2, e);
    VectorXd resLU2 = LUMethod(A2, b2);

    vector<float> errs_identidad{ECM(A2*resJMat2, b2), ECM(A2*resJSum2, b2), ECM(A2*resGSMat2, b2), ECM(A2*resGSSum2, b2), ECM(A2*resLU2, b2)};

	MatrixXi test3(3, 3); 	// MATRIZ SDP --> Sabemos que converge para Gauss-Seidel. No sabemos para Jacobi
							// En el caso de que diverga para Jacobi, se frenará el programa.

	test3 << 2, -1, 0,
			 -1, 2, -1,
			 0, -1, 2;

	MatrixXd A3 = toDouble(test3);
	VectorXd x3 = VectorXd::Random(3);
	VectorXd b3 = A3 * x3;

	cout << "Vamos con la matriz SDP" << endl;
	VectorXd resJMat3 = jacobiMethodHastaConvergencia(A3, b3, e);	// TESTEADO: CONVERGE PARA JACOBI
	VectorXd resJSum3 = jacobiSumatoriaHastaConvergencia(A3, b3, e);
	VectorXd resGSMat3 = gaussSeidelMethodHastaConvergencia(A3, b3, e);
	VectorXd resGSSum3 = gaussSeidelSumatoriaHastaConvergencia(A3, b3, e);
	VectorXd resLU3 = LUMethod(A3, b3);

    vector<float> errs_sdp{ECM(A3*resJMat3, b3), ECM(A3*resJSum3, b3), ECM(A3*resGSMat3, b3), ECM(A3*resGSSum3, b3), ECM(A3*resLU3, b3)};


	MatrixXi test4(5, 5);	// Vamos a hacer una matriz triangular
	test4 << 1, 25, 8, 4, 5,
			 0, 5, 2, 10, 8,
			 0, 0, 1, 7, 8,
			 0, 0, 0, 9, 15,
			 0, 0, 0, 0, 39;
	// Los autovalores de las matrices triangulares son equivalentes a los elementos de su diagonal.
	// Acá podemos ver que el radio espectral de la matriz es 39.
	// Aún así, lo que a nosotros nos interesa es el Radio Espectral de las matrices iterativas de Jacobi y GS
	// Es por esto que es posible que esta matriz converga con Jacobi y Gauss-Seidel
	MatrixXd A4 = toDouble(test4);
	VectorXd x4 = VectorXd::Random(5);
	VectorXd b4 = A4 * x4;

	VectorXd resJMat4 = jacobiMethodHastaConvergencia(A4, b4, e);
	VectorXd resJSum4 = jacobiSumatoriaHastaConvergencia(A4, b4, e);
	VectorXd resGSMat4 = gaussSeidelMethodHastaConvergencia(A4, b4, e);
	VectorXd resGSSum4 = gaussSeidelSumatoriaHastaConvergencia(A4, b4, e);
	VectorXd resLU4 = LUMethod(A4, b4);

    vector<float> errs_triang{ECM(A4*resJMat4, b4), ECM(A4*resJSum4, b4), ECM(A4*resGSMat4, b4), ECM(A4*resGSSum4, b4), ECM(A4*resLU4, b4)};

    generarArchivo_Error_XMult("Errores_Matrices_Particulares.txt", errs_edd, errs_identidad, errs_sdp, errs_triang);

	// AHORA QUEREMOS UNA MATRIZ QUE DIVERJA PARA VER EL BUEN FUNCIONAMIENTO DE NUESTRO PROGRAMA
	MatrixXi test5(3, 3);
	test5 << 1, 2, -2,
			 1, 1, 1,
			 2, 2, 1;	// Radio espectral 2.0004 para la matriz de iteración de Gauss-Seidel ! --> DIVERGE

	MatrixXd A5 = toDouble(test5);
	VectorXd x5 = VectorXd::Random(3);
	VectorXd b5 = A5 * x5;

	VectorXd resJMat5 = jacobiMethodHastaConvergencia(A5, b5, e);
	VectorXd resJSum5 = jacobiSumatoriaHastaConvergencia(A5, b5, e);
	VectorXd resGSMat5 = gaussSeidelMethodHastaConvergencia(A5, b5, e);	// TESTING: DIVERGE PARA GAUSS-SEIDEL
	VectorXd resGSSum5 = gaussSeidelSumatoriaHastaConvergencia(A5, b5, e);
	VectorXd resLU5 = LUMethod(A5, b5);

    // para esta no genero vector de errores ya que diverge

    return 0;
}
