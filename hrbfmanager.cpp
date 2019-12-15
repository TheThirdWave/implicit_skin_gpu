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
    isoVals.resize(plen);
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

    for(int i = 0; i < plen; i++)
    {
        isoVals[i] = eval(points[i*3], points[i*3+1], points[i*3+2], &identity);
    }
    return true;
}

float HRBFManager::eval(float x, float y, float z, rawMat4x4* invMats)
{
    float fs[numHRBFS];
    float final = 0;
    //get input position values from each hrbf.
    for(int i = 0; i < numHRBFS; i++)
    {
        Eigen::Vector4f hrbfSpace;
        hrbfSpace << x, y, z, 1.0f;
        Eigen::Matrix4f invMat;
        invMat <<   invMats[i][0][0], invMats[i][0][1], invMats[i][0][2], invMats[i][0][3],
                    invMats[i][1][0], invMats[i][1][1], invMats[i][1][2], invMats[i][1][3],
                    invMats[i][2][0], invMats[i][2][1], invMats[i][2][2], invMats[i][2][3],
                    invMats[i][3][0], invMats[i][3][1], invMats[i][3][2], invMats[i][3][3];
        hrbfSpace = invMat * hrbfSpace;
        fs[i] = hrbfs[i].eval(hrbfSpace(0), hrbfSpace(1), hrbfSpace(2));
    }
    //blend the values to get the final global function value.
    //(right now I'm just taking the max).
    for(int i = 0; i < numHRBFS-1; i++)
    {
        if(fs[i] > final) final = fs[i];
    }
    return final;
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
