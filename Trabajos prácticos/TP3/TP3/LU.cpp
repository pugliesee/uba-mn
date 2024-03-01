#include "Eigen/Dense"
#include "Eigen/LU"
#include <cstdio>
#include "LU.h"

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;
using Eigen::PartialPivLU;


VectorXd LUMethod(MatrixXd &A, VectorXd &b){
    PartialPivLU<MatrixXd> lu(A);
    VectorXd x = lu.solve(b);
    return x;
}