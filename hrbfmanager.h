#ifndef HRBFMANAGER_H
#define HRBFMANAGER_H

#include <../eigen-eigen-323c052e1731/Eigen/Dense>
#include <vector>
#include <iostream>

#include "hrbfgenerator.h"

class HRBFManager
{
private:
    bool recalc;
    int numHRBFS;
    HRBFGenerator* hrbfs;
    std::vector<float> isoVals;


public:
    HRBFManager();
    ~HRBFManager();
    void createHRBFS(int n);
    void clearHRBFS();
    bool initHRBFS(float points[], int plen, float normals[], int nlen, std::vector<int> weights[], int wlen, float jointPos[], rawMat4x4 *invMats, int jlen);
    float eval(float x, float y, float z, rawMat4x4 *invMats, int *maxIdx);
    Eigen::Vector3f grad(float x, float y, float z, rawMat4x4 *invMats, int matIdx);
    std::vector<float> adjustToHRBF(float x, float y, float z, rawMat4x4 *invMats, int idx);
    void setNeedRecalc(bool r);
    bool getNeedRecalc();
    int getNumHRBFS();
};

#endif // HRBFMANAGER_H
