#include "Eigen/Dense"
#include <cstdio>
#include <tuple>
#include <vector>
#include "DLU.h"

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;

tuple<MatrixXd, MatrixXd> matrizIteracionJacobi(tuple<MatrixXd, MatrixXd, MatrixXd> &DLU);

VectorXd jacobiMethod(MatrixXd &A, VectorXd &b, int N, float epsilon);

VectorXd jacobiMethodXMult(MatrixXd &A, VectorXd &b, VectorXd &x_init, int d, int N, float epsilon);

VectorXd jacobiMethodHastaConvergencia(MatrixXd &A, VectorXd &b, float epsilon);

VectorXd jacobiSumatoria(MatrixXd &A, VectorXd &b, int N, float epsilon);

VectorXd jacobiSumatoriaXMult(MatrixXd &A, VectorXd &b, VectorXd &x_init, int d, int N, float epsilon);

VectorXd jacobiSumatoriaHastaConvergencia(MatrixXd &A, VectorXd &b, float epsilon);
