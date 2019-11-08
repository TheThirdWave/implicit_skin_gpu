#ifndef HRBFGENERATOR_H
#define HRBFGENERATOR_H

#include <../eigen-eigen-323c052e1731/Eigen/Dense>
#include <vector>
#include <iostream>

using Eigen::MatrixXf;
using Eigen::VectorXf;

class HRBFGenerator
{
private:
    MatrixXf coefficients;
    VectorXf unknowns;
    VectorXf results;
    VectorXf mPoints;


    float smoothfunc(float x, float y, float z);
    float derivx(float x, float y, float z);
    float derivy(float x, float y, float z);
    float derivz(float x, float y, float z);
    float h00(float x, float y, float z);
    float h01(float x, float y, float z);
    float h02(float x, float y, float z);
    float h10(float x, float y, float z);
    float h11(float x, float y, float z);
    float h12(float x, float y, float z);
    float h20(float x, float y, float z);
    float h21(float x, float y, float z);
    float h22(float x, float y, float z);

public:
    HRBFGenerator(float points[], int plen, float normals[], int nlen);
    ~HRBFGenerator();


    float eval(float x, float y, float z);

    void solve();

    MatrixXf* getCoefficients();
    VectorXf* getUnknowns();
    VectorXf* getResults();
};

#endif // HRBFGENERATOR_H
