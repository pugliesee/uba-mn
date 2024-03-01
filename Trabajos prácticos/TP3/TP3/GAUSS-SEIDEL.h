#include "Eigen/Dense"
#include <cstdio>
#include <tuple>
#include <vector>
#include "DLU.h"

using namespace Eigen;

using Eigen::MatrixXd;
using Eigen::VectorXd;

tuple<MatrixXd, MatrixXd> matrizIteracionGS(tuple<MatrixXd, MatrixXd, MatrixXd> &DLU);

VectorXd gaussSeidelMethod(MatrixXd &A, VectorXd &b, int N, float epsilon);

VectorXd gaussSeidelMethodXMult(MatrixXd &A, VectorXd &b, VectorXd &x_init, float d, int N, float epsilon);

VectorXd gaussSeidelMethodHastaConvergencia(MatrixXd &A, VectorXd &b, float epsilon);

VectorXd gaussSeidelSumatoria(MatrixXd &A, VectorXd &b, int N, float epsilon);

VectorXd gaussSeidelSumatoriaXMult(MatrixXd &A, VectorXd &b, VectorXd &x_init, float d, int N, float epsilon);

VectorXd gaussSeidelSumatoriaHastaConvergencia(MatrixXd &A, VectorXd &b, float epsilon);