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

bool HRBFManager::initHRBFS(float points[], int plen, float normals[], int nlen, int weights[], int wlen, float jointPos[], int jlen)
{
    std::vector<float> pts[numHRBFS];
    std::vector<float> norms[numHRBFS];
    //sort out points and normals by weight
    for(int i = 0; i < wlen; i++)
    {
        pts[weights[i]].push_back(points[i*3]);
        pts[weights[i]].push_back(points[i*3+1]);
        pts[weights[i]].push_back(points[i*3+2]);
        norms[weights[i]].push_back(normals[i*3]);
        norms[weights[i]].push_back(normals[i*3+1]);
        norms[weights[i]].push_back(normals[i*3+2]);
    }
    for(int i = 0; i < numHRBFS; i++)
    {
        hrbfs[i].init(pts[i], pts[i].size(), norms[i], norms[i].size(), Eigen::Vector3f(jointPos[i*3], jointPos[i*3+1], jointPos[i*3+1]), Eigen::Vector3f(jointPos[(i+1)*3], jointPos[(i+1)*3+1], jointPos[(i+1)*3+1]));
    }
    return true;
}

float HRBFManager::eval(float x, float y, float z)
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
