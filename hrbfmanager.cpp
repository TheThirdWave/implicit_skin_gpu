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

void HRBFManager::clearHRBFS()
{
    delete[] hrbfs;
    isoVals.clear();
    numHRBFS = 0;
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
        hrbfs[i].init(pts[i], pts[i].size(), norms[i], norms[i].size(), Eigen::Vector3f(jointPos[i*3], jointPos[i*3+1], jointPos[i*3+2]), Eigen::Vector3f(jointPos[(i+1)*3], jointPos[(i+1)*3+1], jointPos[(i+1)*3+2]));
        std::cout << "solvehrbf:" << i << std::endl;
        hrbfs[i].solve();
    }

    int* emptyIdx = new int;
    for(int i = 0; i < plen; i++)
    {
        //std::cout << "get isoval: " << i << std::endl;
        isoVals[i] = eval(points[i*3], points[i*3+1], points[i*3+2], &identity, emptyIdx);
    }
    return true;
}

float HRBFManager::eval(float x, float y, float z, rawMat4x4* invMats, int* maxIdx)
{
    float fs[numHRBFS];
    float final = -1;
    //get input position values from each hrbf.
    for(int i = 0; i < numHRBFS; i++)
    {
        if(hrbfs[i].getNumMPoints() <= 0) continue;

        Eigen::Vector4f hrbfSpace;
        hrbfSpace << x, y, z, 1.0f;
        Eigen::Matrix4f invMat;
        invMat <<   invMats[i][0][0], invMats[i][0][1], invMats[i][0][2], invMats[i][0][3],
                    invMats[i][1][0], invMats[i][1][1], invMats[i][1][2], invMats[i][1][3],
                    invMats[i][2][0], invMats[i][2][1], invMats[i][2][2], invMats[i][2][3],
                    invMats[i][3][0], invMats[i][3][1], invMats[i][3][2], invMats[i][3][3];
        hrbfSpace = hrbfSpace.transpose() * invMat;
        std::cout << "EVAL " << i << " HRBFSPACE: " << hrbfSpace << std::endl;
        fs[i] = hrbfs[i].eval(hrbfSpace(0), hrbfSpace(1), hrbfSpace(2));
        //std::cout << fs[i] << " ";
    }
    //std::cout << std::endl;
    //blend the values to get the final global function value.
    //(right now I'm just taking the max).
    for(int i = 0; i < numHRBFS; i++)
    {
        if(hrbfs[i].getNumMPoints() <= 0) continue;
        if(fs[i] > final)
        {
            final = fs[i];
            *maxIdx = i;
            //std::cout << maxIdx << std::endl;
            //std::cout << *maxIdx << std::endl;
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
    //std::cout << "matIdx: " << matIdx << std::endl;
    invMat <<   invMats[i][0][0], invMats[i][0][1], invMats[i][0][2], invMats[i][0][3],
                invMats[i][1][0], invMats[i][1][1], invMats[i][1][2], invMats[i][1][3],
                invMats[i][2][0], invMats[i][2][1], invMats[i][2][2], invMats[i][2][3],
                invMats[i][3][0], invMats[i][3][1], invMats[i][3][2], invMats[i][3][3];
    hrbfSpace = hrbfSpace.transpose() * invMat;
    Eigen::Vector3f grad = hrbfs[i].grad(hrbfSpace(0), hrbfSpace(1), hrbfSpace(2));
    return grad;
}

std::vector<float> HRBFManager::adjustToHRBF(float x, float y, float z, rawMat4x4 *invMats, int idx)
{
    std::cout << "IN POINT: " << x << " " << y << " " << z << std::endl;
    fflush(stdout);
    int* maxIdx = new int;
    //std::cout << maxIdx << std::endl;
    //get iso value for current point.
    float iso = eval(x, y, z, invMats, maxIdx);
    std::cout << "EVALOUT: " << iso << std::endl;
    fflush(stdout);
    float isoDiff = iso - isoVals[idx];
    std::cout << "ISODIFF: " << isoDiff << std::endl;
    //translate input point to hrbf space.
    Eigen::Vector4f hrbfSpace;
    Eigen::Vector4f worldSpace;
    hrbfSpace << x, y, z, 1.0f;
    Eigen::Matrix4f invMat;
    invMat <<   invMats[*maxIdx][0][0], invMats[*maxIdx][0][1], invMats[*maxIdx][0][2], invMats[*maxIdx][0][3],
                invMats[*maxIdx][1][0], invMats[*maxIdx][1][1], invMats[*maxIdx][1][2], invMats[*maxIdx][1][3],
                invMats[*maxIdx][2][0], invMats[*maxIdx][2][1], invMats[*maxIdx][2][2], invMats[*maxIdx][2][3],
                invMats[*maxIdx][3][0], invMats[*maxIdx][3][1], invMats[*maxIdx][3][2], invMats[*maxIdx][3][3];
    std::cout << "INVMAT: " << invMat << std::endl;
    hrbfSpace = hrbfSpace.transpose() * invMat;
    std::cout << "HRBFSPACE: " << hrbfSpace << std::endl;
    Eigen::Matrix4f origMat = invMat.inverse();
    std::cout << "ORIGMAT: " << origMat << std::endl;
    Eigen::Vector3f pt(x,y,z);
    Eigen::Vector3f ptGrad;
    Eigen::Vector3f oldGrad;

    //Get the gradient at the point we're evaluating.
    //std::cout << maxIdx << std::endl;
    //std::cout << *maxIdx << std::endl;
    fflush(stdout);
    ptGrad = grad(x, y, z, invMats, *maxIdx);
    oldGrad = ptGrad;

    //We'll need to get the difference between the gradient angles of the current iteration and the previous iteration, it should be 0 for the first iteration though.
    float gAngle = ptGrad.dot(oldGrad);
    gAngle = gAngle/(ptGrad.norm() * oldGrad.norm());
    gAngle = std::acos(gAngle);
    gAngle = gAngle/(2*PI) * 360;
    if(std::isnan(gAngle)) gAngle = 0;
    std::cout << "GANGLE: " << gAngle << std::endl;
    //we move the point closer to it's iso value using Newton iterations. Break when we're close enough or hit an angle discontinuity
    //(The angle discontinuity suggests we've hit another HRBF and so we should stop to prevent self-intersection)
    while(isoDiff > MINISODIFFERENCE && gAngle < DISCONTINUITYANGLE)
    {
        std::cout << "ISODIFF: " << isoDiff << std::endl;
        pt = pt + PROJECTIONNUM * (isoDiff) * (ptGrad/(ptGrad.norm()*ptGrad.norm()));

        //after adjusting the hrbfspace point, translate pt back to world space since that's what eval expects.
        hrbfSpace << pt(0), pt(1), pt(2), 1.0f;
        std::cout << "HRBFSPACE: " << hrbfSpace << std::endl;
        worldSpace = hrbfSpace.transpose() * origMat;
        std::cout << "WORLDSPACE: " << worldSpace << std::endl;

        iso = eval(worldSpace(0),worldSpace(1),worldSpace(2), invMats, maxIdx);
        std::cout << "EVALOUT: " << iso << std::endl;
        fflush(stdout);
        isoDiff = iso - isoVals[idx];
        std::cout << "ISODIFF: " << isoDiff << std::endl;
        oldGrad = ptGrad;
        //std::cout << *maxIdx << std::endl;
        ptGrad = grad(worldSpace(0),worldSpace(1),worldSpace(2), invMats, *maxIdx);
        float gAngle = ptGrad.dot(oldGrad);
        gAngle = gAngle/(ptGrad.norm() * oldGrad.norm());
        gAngle = std::acos(gAngle);
        gAngle = gAngle/(2*PI) * 360;
    }
    return std::vector<float> { pt(0), pt(1), pt(2) };
}

void HRBFManager::setNeedRecalc(bool r)
{
    for(int i = 0; i < numHRBFS; i++)
    {
        hrbfs[i].setRecalc(true);
    }
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
