#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>

using namespace std;

void generarArchivo_Error_Tiempo(string filename, vector<float> res_J_Mat, vector<float> res_J_Sum, vector<float> res_GS_Mat, vector<float> res_GS_Sum, vector<float> res_LU);

void generarArchivo_Error_XMult(string filename, vector<float> res_J_Mat, vector<float> res_J_Sum, vector<float> res_GS_Mat, vector<float> res_GS_Sum);
