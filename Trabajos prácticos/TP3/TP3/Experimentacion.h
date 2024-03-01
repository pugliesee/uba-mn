#include <iostream>
#include <cstdio>
#include "Eigen/Dense"
#include "JACOBI.h"
#include "GAUSS-SEIDEL.h"
#include "DLU.h"
#include "LU.h"
#include "Auxiliares.h"
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

using Eigen::MatrixXd;
using Eigen::VectorXd;

void Experimentacion_XMult_J_Mat(MatrixXd &A, VectorXd &x, float epsilon, vector<int> iters);

void Experimentacion_XMult_J_Sum(MatrixXd &A, VectorXd &x, float epsilon, vector<int> iters);

void Experimentacion_XMult_GS_Mat(MatrixXd &A, VectorXd &x, float epsilon, vector<int> iters);

void Experimentacion_XMult_GS_Sum(MatrixXd &A, VectorXd &x, float epsilon, vector<int> iters);

void correrExperimentacion();