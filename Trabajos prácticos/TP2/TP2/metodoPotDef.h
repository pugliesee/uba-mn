// #include <sqlite3.h>
#include <cstdio>
#include <eigen3/Eigen/Dense>
#include <tuple>
#include <vector>
#include "metodoPotDef.cpp"


using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;

MatrixXd deflacion(MatrixXd A, VectorXd v, float x);

tuple<MatrixXd, VectorXd, double> metodoPotencia(MatrixXd A, int N, double e);

int metodoPotenciaDeflacion(string matriz, int N, string values_filename, string vectors_filename, string times_filename);
    