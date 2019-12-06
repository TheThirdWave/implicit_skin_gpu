#include "hrbfmanager.h"

HRBFManager::HRBFManager()
{
    numHRBFS = 0;
    recalc = true;
    hrbfs = NULL;
}

HRBFManager::~HRBFManager()
{
    delete[] hrbfs;
}

void HRBFManager::createHRBFS(int n)
{
    hrbfs = new HRBFGenerator[numHRBFS];
    numHRBFS = n;
}

bool HRBFManager::initHRBFS(float points[], int plen, float normals[], int nlen, float weights, int wlen)
{

}

bool HRBFManager::getNeedRecalc()
{
    recalc = false;
    for(int i = 0; i < numHRBFS; i++)
    {
        if(hrbfs[i].getNeedRecalc() == true) recalc = true;
    }
    return recalc;
}

int HRBFManager::getNumHRBFS()
{
    return numHRBFS;
}
