#include "hrbfgenerator.h"

#define VECTORLEN 3

using namespace Eigen;

HRBFGenerator::HRBFGenerator()
{
    recalc = true;
}

HRBFGenerator::HRBFGenerator(float points[], int plen, float normals[], int nlen)
{
    if(plen/3 != nlen/3)
    {
        std::cout << "HRBFGEN ERROR: input arrays of different length!" << std::endl;
        return;
    }
    int sidelen = plen/3;
    coefficients.resize(sidelen*4, sidelen*4);
    unknowns.resize(sidelen*4);
    results.resize(sidelen*4);
    mPoints.resize(sidelen*3);

    for(int i = 0; i < 4*sidelen; i += 4)
    {
        int idxi = i / 4 * VECTORLEN;
        for(int j = 0; j < 4*sidelen; j += 4)
        {
            int idxj = j / 4 * VECTORLEN;
            coefficients(i,j) = smoothfunc(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i,j+1) = -derivx(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i,j+2) = -derivy(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i,j+3) = -derivz(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+1,j) = derivx(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+1,j+1) = -h00(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+1,j+2) = -h01(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+1,j+3) = -h02(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+2,j) = derivy(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+2,j+1) = -h10(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+2,j+2) = -h11(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+2,j+3) = -h12(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+3,j) = derivz(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+3,j+1) = -h20(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+3,j+2) = -h21(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+3,j+3) = -h22(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
        }
        results(i) = 0;
        results(i+1) = normals[idxi];
        results(i+2) = normals[idxi+1];
        results(i+3) = normals[idxi+2];
        mPoints(idxi) = points[idxi];
        mPoints(idxi+1) = points[idxi+1];
        mPoints(idxi+2) = points[idxi+2];
    }
    recalc = false;
    return;
}

HRBFGenerator::~HRBFGenerator()
{
    coefficients.resize(0,0);
    unknowns.resize(0);
    results.resize(0);
}

void HRBFGenerator::init(float points[], int plen, float normals[], int nlen)
{
    if(plen/3 != nlen/3)
    {
        std::cout << "HRBFGEN ERROR: input arrays of different length!" << std::endl;
        return;
    }
    int sidelen = plen/3;
    coefficients.resize(sidelen*4, sidelen*4);
    unknowns.resize(sidelen*4);
    results.resize(sidelen*4);
    mPoints.resize(sidelen*3);

    for(int i = 0; i < 4*sidelen; i += 4)
    {
        int idxi = i / 4 * VECTORLEN;
        for(int j = 0; j < 4*sidelen; j += 4)
        {
            int idxj = j / 4 * VECTORLEN;
            coefficients(i,j) = smoothfunc(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i,j+1) = -derivx(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i,j+2) = -derivy(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i,j+3) = -derivz(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+1,j) = derivx(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+1,j+1) = -h00(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+1,j+2) = -h01(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+1,j+3) = -h02(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+2,j) = derivy(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+2,j+1) = -h10(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+2,j+2) = -h11(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+2,j+3) = -h12(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+3,j) = derivz(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+3,j+1) = -h20(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+3,j+2) = -h21(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            coefficients(i+3,j+3) = -h22(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
        }
        results(i) = 0;
        results(i+1) = normals[idxi];
        results(i+2) = normals[idxi+1];
        results(i+3) = normals[idxi+2];
        mPoints(idxi) = points[idxi];
        mPoints(idxi+1) = points[idxi+1];
        mPoints(idxi+2) = points[idxi+2];
    }
    recalc = false;
    return;
}

float HRBFGenerator::eval(float x, float y, float z)
{
    Vector3f p(x, y, z);
    float out = 0;
    for(int i = 0; i < mPoints.size()/3; i++)
    {
        int mpidx = i * 3;
        int cidx = i * 4;
        Vector3f vk(mPoints(mpidx), mPoints(mpidx+1), mPoints(mpidx+2));
        float alpha = unknowns(cidx);
        Vector3f beta(unknowns(cidx+1), unknowns(cidx+2), unknowns(cidx+3));
        Vector3f diff = p - vk;
        Vector3f grad(derivx(diff(0), diff(1), diff(2)), derivy(diff(0), diff(1), diff(2)), derivz(diff(0), diff(1), diff(2)));
        out += alpha*smoothfunc(diff(0), diff(1), diff(2)) - beta.dot(grad);

    }
    return out;
}

void HRBFGenerator::solve()
{
    unknowns = coefficients.colPivHouseholderQr().solve(results);
}

float HRBFGenerator::smoothfunc(float x, float y, float z)
{
    return pow(sqrt(x*x+y*y+z*z),3);
}

float HRBFGenerator::derivx(float x, float y, float z)
{
    return 3*x*sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::derivy(float x, float y, float z)
{
    return 3*y*sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::derivz(float x, float y, float z)
{
    return 3*z*sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::h00(float x, float y, float z)
{
    if(x == 0 && y == 0 && z == 0) return 0;
    return 3*(2*x*x+y*y+z*z)/sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::h01(float x, float y, float z)
{
    if(x == 0 && y == 0 && z == 0) return 0;
    return 3*x*y/sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::h02(float x, float y, float z)
{
    if(x == 0 && y == 0 && z == 0) return 0;
    return 3*x*z/sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::h10(float x, float y, float z)
{
    if(x == 0 && y == 0 && z == 0) return 0;
    return 3*x*y/sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::h11(float x, float y, float z)
{
    if(x == 0 && y == 0 && z == 0) return 0;
    return 3*(x*x+2*y*y+z*z)/sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::h12(float x, float y, float z)
{
    if(x == 0 && y == 0 && z == 0) return 0;
    return 3*z*y/sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::h20(float x, float y, float z)
{
    if(x == 0 && y == 0 && z == 0) return 0;
    return 3*x*z/sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::h21(float x, float y, float z)
{
    if(x == 0 && y == 0 && z == 0) return 0;
    return 3*y*z/sqrt(x*x+y*y+z*z);
}

float HRBFGenerator::h22(float x, float y, float z)
{
    if(x == 0 && y == 0 && z == 0) return 0;
    return 3*(x*x+y*y+2*z*z)/sqrt(x*x+y*y+z*z);
}

MatrixXf* HRBFGenerator::getCoefficients()
{
    return &coefficients;
}

VectorXf* HRBFGenerator::getUnknowns()
{
    return &unknowns;
}

VectorXf* HRBFGenerator::getResults()
{
    return &results;
}

bool HRBFGenerator::getNeedRecalc()
{
    return recalc;
}

void HRBFGenerator::setRecalc(bool r)
{
    recalc = r;
}
