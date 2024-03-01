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

MatrixXd read_txt(string matriz);

double metodoPotencia(MatrixXd &A, int N, double e);

float ECM(VectorXd x, VectorXd b);

MatrixXd toDouble(MatrixXi &entry);