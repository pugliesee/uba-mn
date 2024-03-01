#include "Eigen/Dense"
#include "Eigen/LU"
#include <cstdio>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;

VectorXd LUMethod(MatrixXd &A, VectorXd &b);