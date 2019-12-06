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
    bool buildHRBFS(float points[], int plen, float normals[], int nlen);
    bool getNeedRecalc();
};

#endif // HRBFMANAGER_H
