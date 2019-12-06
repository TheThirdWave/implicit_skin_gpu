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


public:
    HRBFManager();
    ~HRBFManager();
    void createHRBFS(int n);
    bool initHRBFS(float points[], int plen, float normals[], int nlen, float weights, int wlen);
    bool getNeedRecalc();
    int getNumHRBFS();
};

#endif // HRBFMANAGER_H
