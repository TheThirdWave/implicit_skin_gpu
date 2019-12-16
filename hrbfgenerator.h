#ifndef HRBFGENERATOR_H
#define HRBFGENERATOR_H

#include <../eigen-eigen-323c052e1731/Eigen/Dense>
#include <vector>
#include <iostream>

using Eigen::MatrixXf;
using Eigen::VectorXf;

typedef float rawMat4x4[4][4];
static float identity[4][4] = {   1,0,0,0,
                                  0,1,0,0,
                                  0,0,1,0,
                                  0,0,0,1  };

class HRBFGenerator
{
private:
    MatrixXf* coefficients;
    VectorXf* unknowns;
    VectorXf* results;
    VectorXf* mPoints;
    float radius;


    bool recalc;
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
    HRBFGenerator();
    HRBFGenerator(std::vector<float> points, int plen, std::vector<float> normals, int nlen, Eigen::Vector3f startJoint, Eigen::Vector3f endJoint);
    ~HRBFGenerator();


    float eval(float x, float y, float z);
    Eigen::Vector3f grad(float x, float y, float z);

    void init(std::vector<float> points, int plen, std::vector<float> normals, int nlen, Eigen::Vector3f startJoint, Eigen::Vector3f endJoint);
    void solve();

    MatrixXf* getCoefficients();
    VectorXf* getUnknowns();
    VectorXf* getResults();
    bool getNeedRecalc();

    void setRecalc(bool r);
};

#endif // HRBFGENERATOR_H
