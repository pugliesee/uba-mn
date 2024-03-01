#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>

using namespace std;

void generarArchivo_Error_Tiempo(string filename, vector<float> res_J_Mat, vector<float> res_J_Sum, vector<float> res_GS_Mat, vector<float> res_GS_Sum, vector<float> res_LU){
    ofstream fout(filename);
    for(int n = 0; n < res_J_Mat.size(); n++){
        if(n == res_J_Mat.size() - 1){
            fout << res_J_Mat[n] << endl;
        } else {
            fout << res_J_Mat[n] << " ";
        }
    }
    for(int n = 0; n < res_J_Sum.size(); n++){
        if(n == res_J_Sum.size() - 1){
            fout << res_J_Sum[n] << endl;
        } else {
            fout << res_J_Sum[n] << " ";
        }
    }
    for(int n = 0; n < res_GS_Mat.size(); n++){
        if(n == res_GS_Mat.size() - 1){
            fout << res_GS_Mat[n] << endl;
        } else {
            fout << res_GS_Mat[n] << " ";
        }
    }
    for(int n = 0; n < res_GS_Sum.size(); n++){
        if(n == res_GS_Sum.size() - 1){
            fout << res_GS_Sum[n] << endl;
        } else {
            fout << res_GS_Sum[n] << " ";
        }
    }

    for(int n = 0; n < res_LU.size(); n++){
        if(n == res_LU.size() - 1){
            fout << res_LU[n] << endl;
        } else {
            fout << res_LU[n] << " ";
        }
    }

    fout.close();
}

void generarArchivo_Error_XMult(string filename, vector<float> res_x_2, vector<float> res_x_10, vector<float> res_x_100, vector<float> res_x_1000){
    ofstream fout(filename);
    for(int n = 0; n < res_x_2.size(); n++){
        if(n == res_x_2.size() - 1){
            fout << res_x_2[n] << endl;
        } else {
            fout << res_x_2[n] << " ";
        }
    }
    for(int n = 0; n < res_x_10.size(); n++){
        if(n == res_x_10.size() - 1){
            fout << res_x_10[n] << endl;
        } else {
            fout << res_x_10[n] << " ";
        }
    }
    for(int n = 0; n < res_x_100.size(); n++){
        if(n == res_x_100.size() - 1){
            fout << res_x_100[n] << endl;
        } else {
            fout << res_x_100[n] << " ";
        }
    }
    for(int n = 0; n < res_x_1000.size(); n++){
        if(n == res_x_1000.size() - 1){
            fout << res_x_1000[n] << endl;
        } else {
            fout << res_x_1000[n] << " ";
        }
    }

    fout.close();
}