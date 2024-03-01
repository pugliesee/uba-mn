#include "Eigen/Dense"
#include <cstdio>
#include <tuple>
#include <vector>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;

bool sinCerosDiagonal(MatrixXd &A);

tuple<MatrixXd, MatrixXd, MatrixXd> eigenDLU(MatrixXd &A);


