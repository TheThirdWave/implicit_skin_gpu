#include "hrbfmanager.h"

#define PROJECTIONNUM 0.35
#define MINISODIFFERENCE 0.001
#define PI 3.1415926
#define DISCONTINUITYANGLE 55

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
    numHRBFS = n;
    hrbfs = new HRBFGenerator[numHRBFS];
    for(int i = 0; i < numHRBFS; i++)
    {
        hrbfs[i] = HRBFGenerator();
    }
}

bool HRBFManager::initHRBFS(float points[], int plen, float normals[], int nlen, int weights[], int wlen, float jointPos[], int jlen)
{
    std::vector<float> pts[numHRBFS];
    std::vector<float> norms[numHRBFS];
    std::cout << "plen: " << plen << std::endl;
    std::cout << "nlen: " << nlen << std::endl;
    std::cout << "wlen: " << wlen << std::endl;
    std::cout << "jlen: " << jlen << std::endl;
    std::cout << "numhrbfs: " << numHRBFS << std::endl;
    isoVals.resize(plen);
    //sort out points and normals by weight
    for(int i = 0; i < wlen; i++)
    {
        //std::cout << "vdata: " << i << std::endl;
        //std::cout << weights[i] << std::endl;
        pts[weights[i]].push_back(points[i*3]);
        pts[weights[i]].push_back(points[i*3+1]);
        pts[weights[i]].push_back(points[i*3+2]);
        norms[weights[i]].push_back(normals[i*3]);
        norms[weights[i]].push_back(normals[i*3+1]);
        norms[weights[i]].push_back(normals[i*3+2]);
    }
    for(int i = 0; i < numHRBFS; i++)
    {
        std::cout << "inithrbf: " << i << std::endl;
        hrbfs[i].init(pts[i], pts[i].size(), norms[i], norms[i].size(), Eigen::Vector3f(jointPos[i*3], jointPos[i*3+1], jointPos[i*3+1]), Eigen::Vector3f(jointPos[(i+1)*3], jointPos[(i+1)*3+1], jointPos[(i+1)*3+1]));
        std::cout << "solvehrbf:" << i << std::endl;
        hrbfs[i].solve();
    }

    int* emptyIdx = new int;
    for(int i = 0; i < plen; i++)
    {
        std::cout << "get isoval: " << i << std::endl;
        isoVals[i] = eval(points[i*3], points[i*3+1], points[i*3+2], &identity, emptyIdx);
    }
    return true;
}

float HRBFManager::eval(float x, float y, float z, rawMat4x4* invMats, int* maxIdx)
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
        if(fs[i] > final)
        {
            final = fs[i];
            *maxIdx = i;
        }
    }
    return final;
}

Eigen::Vector3f HRBFManager::grad(float x, float y, float z, rawMat4x4 *invMats, int matIdx)
{
    int i = matIdx;
    Eigen::Vector4f hrbfSpace;
    hrbfSpace << x, y, z, 1.0f;
    Eigen::Matrix4f invMat;
    invMat <<   invMats[i][0][0], invMats[i][0][1], invMats[i][0][2], invMats[i][0][3],
                invMats[i][1][0], invMats[i][1][1], invMats[i][1][2], invMats[i][1][3],
                invMats[i][2][0], invMats[i][2][1], invMats[i][2][2], invMats[i][2][3],
                invMats[i][3][0], invMats[i][3][1], invMats[i][3][2], invMats[i][3][3];
    hrbfSpace = invMat * hrbfSpace;
    Eigen::Vector3f grad = hrbfs[i].grad(hrbfSpace(0), hrbfSpace(1), hrbfSpace(2));
    return grad;
}

std::vector<float> HRBFManager::adjustToHRBF(float x, float y, float z, rawMat4x4 *invMats, int idx)
{
    std::cout << x << " " << y << " " << z << std::endl;
    fflush(stdout);
    int* maxIdx = new int;
    float iso = eval(x, y, z, invMats, maxIdx);
    float isoDiff = iso - isoVals[idx];
    Eigen::Vector3f pt(x,y,z);
    Eigen::Vector3f ptGrad;
    Eigen::Vector3f oldGrad;

    //Get the gradient at the point we're evaluating.
    ptGrad = grad(x, y, z, invMats, *maxIdx);
    oldGrad = ptGrad;

    //We'll need to get the difference between the gradient angles of the current iteration and the previous iteration, it should be 0 for the first iteration though.
    float gAngle = ptGrad.dot(oldGrad);
    gAngle = gAngle/(ptGrad.norm() * oldGrad.norm());
    gAngle = std::acos(gAngle);
    gAngle = gAngle/(2*PI) * 360;
    //we move the point closer to it's iso value using Newton iterations. Break when we're close enough or hit an angle discontinuity
    //(The angle discontinuity suggests we've hit another HRBF and so we should stop to prevent self-intersection)
    while(isoDiff > MINISODIFFERENCE && gAngle < DISCONTINUITYANGLE)
    {
        pt = pt + PROJECTIONNUM * (isoDiff) * (ptGrad/(ptGrad.norm()*ptGrad.norm()));

        iso = eval(x, y, z, invMats, maxIdx);
        isoDiff = iso - isoVals[idx];
        oldGrad = ptGrad;
        ptGrad = grad(x, y, z, invMats, *maxIdx);
        float gAngle = ptGrad.dot(oldGrad);
        gAngle = gAngle/(ptGrad.norm() * oldGrad.norm());
        gAngle = std::acos(gAngle);
        gAngle = gAngle/(2*PI) * 360;
    }
    return std::vector<float> { pt(0), pt(1), pt(2) };
}

bool HRBFManager::getNeedRecalc()
{
    std::cout << "GET NEED RECALC: ";
    recalc = false;
    if(numHRBFS <= 0) recalc = true;
    for(int i = 0; i < numHRBFS; i++)
    {
        if(hrbfs[i].getNeedRecalc() == true) recalc = true;
    }
    std::cout << recalc << std::endl;
    return recalc;
}

int HRBFManager::getNumHRBFS()
{
    return numHRBFS;
}
