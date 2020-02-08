#include "hrbfgenerator.h"

#define VECTORLEN 3
#define CULLDISTANCE 0.05

using namespace Eigen;

HRBFGenerator::HRBFGenerator()
{
    coefficients = new MatrixXf(0, 0);
    std::cout << coefficients<< std::endl;
    fflush(stdout);
    unknowns = new VectorXf(0);
    results = new VectorXf(0);
    mPoints = new VectorXf(0);
    recalc = true;
}

HRBFGenerator::HRBFGenerator(std::vector<float> points, int plen, std::vector<float> normals, int nlen, Vector3f startJoint, Vector3f endJoint)
{
    std::cout << "START HRBF CONSTRUCTOR"<< std::endl;
    fflush(stdout);
    if(plen/3 != nlen/3)
    {
        std::cout << "HRBFGEN ERROR: input arrays of different length!" << std::endl;
        return;
    }
    int sidelen = plen/3;
    coefficients->resize(sidelen*4, sidelen*4);
    unknowns->resize(sidelen*4);
    results->resize(sidelen*4);
    mPoints->resize(sidelen*3);
    float largestR = 0;
    float smallestR = -1;


    //get the direction vector of the joint.
    Eigen::Vector3f direction = (endJoint-startJoint);
    Eigen::Vector3f directionNorm = direction / direction.norm();

    //find largest distance and cull points that are too close for reasons.
    for(int i = 0; i < plen; i = i + 3)
    {
        //find the distance of the point from the rigging bone.
        Eigen::Vector3f pt = Eigen::Vector3f(points[i], points[i+1], points[i+2]);
        //map the point onto the plane defined by the start of the joint and it's direction.
        Eigen::Vector3f toPt = pt - startJoint;
        Eigen::Vector3f ptPlane = pt - (toPt.dot(directionNorm) * directionNorm);
        //the length of the vector between the start joint and the mapped point is the radial distance from the bone.
        float dist = (ptPlane-endJoint).norm();
        if(dist > largestR) largestR = dist;
        if(dist < smallestR || smallestR == -1) smallestR = dist;
        float a = (endJoint-startJoint).norm();
        float cullPlane = (toPt.transpose() * direction);
        cullPlane = cullPlane/(a * a);
        if(cullPlane <= CULLDISTANCE || cullPlane >= 1.0 - CULLDISTANCE)
        {
            points.erase(points.begin()+i, points.begin()+i+3);
            plen -= 3;
            normals.erase(normals.begin()+i, normals.begin()+i+3);
            nlen -= 3;
        }
    }
    //update radius
    radius = largestR;

    //Add cap points to each end of the HRBF to allow for smooth deformation at the joints.
    Eigen::Vector3f capPoint = endJoint + (directionNorm * smallestR);
    points.push_back(capPoint(0));
    points.push_back(capPoint(1));
    points.push_back(capPoint(2));
    plen += 3;
    normals.push_back(directionNorm(0));
    normals.push_back(directionNorm(1));
    normals.push_back(directionNorm(2));
    nlen += 3;
    capPoint = startJoint - (directionNorm * smallestR);
    points.push_back(capPoint(0));
    points.push_back(capPoint(1));
    points.push_back(capPoint(2));
    plen += 3;
    normals.push_back(directionNorm(0));
    normals.push_back(directionNorm(1));
    normals.push_back(directionNorm(2));
    nlen += 3;


    //BUILD ALL THE BIG MATRICES YOU NEED TO SOLVE FOR THE STUFF.
    for(int i = 0; i < 4*sidelen; i += 4)
    {
        int idxi = i / 4 * VECTORLEN;
        for(int j = 0; j < 4*sidelen; j += 4)
        {
            int idxj = j / 4 * VECTORLEN;
            ((*coefficients))(i,j) = smoothfunc(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i,j+1) = -derivx(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i,j+2) = -derivy(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i,j+3) = -derivz(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+1,j) = derivx(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+1,j+1) = -h00(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+1,j+2) = -h01(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+1,j+3) = -h02(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+2,j) = derivy(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+2,j+1) = -h10(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+2,j+2) = -h11(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+2,j+3) = -h12(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+3,j) = derivz(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+3,j+1) = -h20(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+3,j+2) = -h21(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+3,j+3) = -h22(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
        }
        ((*results))(i) = 0;
        (*results)(i+1) = normals[idxi];
        (*results)(i+2) = normals[idxi+1];
        (*results)(i+3) = normals[idxi+2];
        (*mPoints)(idxi) = points[idxi];
        (*mPoints)(idxi+1) = points[idxi+1];
        (*mPoints)(idxi+2) = points[idxi+2];
    }
    recalc = false;
    return;
}

HRBFGenerator::~HRBFGenerator()
{
    coefficients->resize(0,0);
    unknowns->resize(0);
    results->resize(0);
    mPoints->resize(0);
}

void HRBFGenerator::init(std::vector<float> points, int plen, std::vector<float> normals, int nlen, Eigen::Vector3f startJoint, Eigen::Vector3f endJoint)
{
    if(plen/3 != nlen/3)
    {
        std::cout << "HRBFGEN ERROR: input arrays of different length!" << std::endl;
        return;
    }
    std::cout << "START HRBF INIT"<< std::endl;
    float largestR = 0;
    float smallestR = -1;


    std::cout << "get direction vector of joint."  << std::endl;
    fflush(stdout);
    //get the direction vector of the joint.
    Eigen::Vector3f direction = (endJoint-startJoint);
    Eigen::Vector3f directionNorm = direction / direction.norm();

    std::cout << "find largest distance and cull points that are too close for reasons." << std::endl;
    fflush(stdout);
    //find largest distance and cull points that are too close for reasons.
    for(int i = 0; i < plen; i = i + 3)
    {
        //find the distance of the point from the rigging bone.
        Eigen::Vector3f pt = Eigen::Vector3f(points[i], points[i+1], points[i+2]);
        std::cout << "PT: " << pt << std::endl;
        //map the point onto the plane defined by the start of the joint and it's direction.
        Eigen::Vector3f toPt = pt - startJoint;
        Eigen::Vector3f ptPlane = pt - (toPt.dot(directionNorm) * directionNorm);
        //the length of the vector between the start joint and the mapped point is the radial distance from the bone.
        float dist = (ptPlane-startJoint).norm();
        if(dist > largestR) largestR = dist;
        if(dist < smallestR || smallestR == -1) smallestR = dist;
        float a = (endJoint-startJoint).norm();
        float cullPlane = (toPt.transpose() * direction);
        cullPlane = cullPlane/(a * a);
        std::cout << "cullPlane: " << cullPlane << std::endl;
        if(cullPlane <= CULLDISTANCE || cullPlane >= 1.0 - CULLDISTANCE)
        {
            std::cout << "directionNorm: " << directionNorm << std::endl;
            std::cout << "topt: " << toPt << std::endl;
            std::cout << "ptPlane: " << ptPlane << std::endl;
            std::cout << "dist: " << dist << std::endl;
            std::cout << "largestR: " << largestR << std::endl;
            std::cout << "smallestR: " << smallestR << std::endl;
            std::cout << "a: " << a << std::endl;
            points.erase(points.begin()+i, points.begin()+i+3);
            plen -= 3;
            normals.erase(normals.begin()+i, normals.begin()+i+3);
            nlen -= 3;
            i -= 3;
        }
    }

    std::cout << "update radius" << std::endl;
    fflush(stdout);
    //update radius
    radius = largestR;
    std::cout << "RADIUS: " << radius << std::endl;

    std::cout << "add cap points" << std::endl;
    fflush(stdout);
    //Add cap points to each end of the HRBF to allow for smooth deformation at the joints.
    Eigen::Vector3f capPoint = endJoint + (directionNorm * smallestR);
    std::cout << "CPBack: " << capPoint << std::endl;
    points.push_back(capPoint(0));
    points.push_back(capPoint(1));
    points.push_back(capPoint(2));
    plen += 3;
    normals.push_back(directionNorm(0));
    normals.push_back(directionNorm(1));
    normals.push_back(directionNorm(2));
    nlen += 3;
    capPoint = startJoint - (directionNorm * smallestR);
    std::cout << "CPFront: " << capPoint << std::endl;
    points.push_back(capPoint(0));
    points.push_back(capPoint(1));
    points.push_back(capPoint(2));
    plen += 3;
    normals.push_back(-directionNorm(0));
    normals.push_back(-directionNorm(1));
    normals.push_back(-directionNorm(2));
    nlen += 3;

    fflush(stdout);
    int sidelen = plen/3;
    std::cout << "sidelen: "<< sidelen << std::endl;
    std::cout << "StartJoint: " << startJoint << std::endl;
    std::cout << "EndJoint: " << endJoint << std::endl;
    fflush(stdout);
    std::cout << coefficients << std::endl;
    if(coefficients != nullptr)
    {
        std::cout << "DELETING COEFFICIENTS"<< std::endl;
        fflush(stdout);
        delete coefficients;
        coefficients = nullptr;
    }
    coefficients = new MatrixXf(sidelen*4, sidelen*4);
    std::cout << "resize coeff"<< std::endl;
    fflush(stdout);
    if(unknowns != nullptr)
    {
        delete unknowns;
        unknowns = nullptr;
    }
    unknowns = new VectorXf(sidelen*4);
    std::cout << "resize unkown"<< std::endl;
    fflush(stdout);
    if(results != nullptr)
    {
        delete results;
        results = nullptr;
    }
    results = new VectorXf(sidelen*4);
    std::cout << "resize results"<< std::endl;
    fflush(stdout);
    if(mPoints != nullptr)
    {
        delete mPoints;
        mPoints = nullptr;
    }
    mPoints = new VectorXf(sidelen*3);
    std::cout << "resize mPoints"<< std::endl;
    fflush(stdout);


    std::cout << "build matrices" << std::endl;
    fflush(stdout);
    //BUILD ALL THE BIG MATRICES YOU NEED TO SOLVE FOR THE STUFF.
    for(int i = 0; i < 4*sidelen; i += 4)
    {
        int idxi = i / 4 * VECTORLEN;
        for(int j = 0; j < 4*sidelen; j += 4)
        {
            //std::cout << i << " " << j << std::endl;
            int idxj = j / 4 * VECTORLEN;
            (*coefficients)(i,j) = smoothfunc(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i,j+1) = -derivx(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i,j+2) = -derivy(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i,j+3) = -derivz(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+1,j) = derivx(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+1,j+1) = -h00(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+1,j+2) = -h01(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+1,j+3) = -h02(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+2,j) = derivy(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+2,j+1) = -h10(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+2,j+2) = -h11(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+2,j+3) = -h12(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+3,j) = derivz(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+3,j+1) = -h20(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+3,j+2) = -h21(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
            (*coefficients)(i+3,j+3) = -h22(points[idxi] - points[idxj], points[idxi+1] - points[idxj+1], points[idxi+2] - points[idxj+2]);
        }
        (*results)(i) = 0;
        (*results)(i+1) = normals[idxi];
        (*results)(i+2) = normals[idxi+1];
        (*results)(i+3) = normals[idxi+2];
        (*mPoints)(idxi) = points[idxi];
        (*mPoints)(idxi+1) = points[idxi+1];
        (*mPoints)(idxi+2) = points[idxi+2];
    }
    recalc = false;
    return;
}

float HRBFGenerator::eval(float x, float y, float z)
{
    //std::cout << "STARTING EVAL!!" << std::endl;
    //std::cout << "size = " << mPoints->size() << std::endl;
    fflush(stdout);
    Vector3f p(x, y, z);
    float out = 0;
    for(int i = 0; i < mPoints->size()/3; i++)
    {
        int mpidx = i * 3;
        int cidx = i * 4;
        Vector3f vk((*mPoints)(mpidx), (*mPoints)(mpidx+1), (*mPoints)(mpidx+2));
        float alpha = (*unknowns)(cidx);
        Vector3f beta((*unknowns)(cidx+1), (*unknowns)(cidx+2), (*unknowns)(cidx+3));
        Vector3f diff = p - vk;
        Vector3f grad(derivx(diff(0), diff(1), diff(2)), derivy(diff(0), diff(1), diff(2)), derivz(diff(0), diff(1), diff(2)));
        //std::cout << alpha << std::endl;
        //std::cout << beta(0) << "," << beta(1) << "," << beta(2) << " " << diff(0) << "," << diff(1) << "," << diff(2) << " " << grad(0) << "," << grad(1) << "," << grad(2) << std::endl;
        out += alpha*smoothfunc(diff(0), diff(1), diff(2)) - beta.dot(grad);
        //std::cout << out << std::endl;

    }
    std::cout << out << std::endl;
    if(out < -radius) return 1;
    if(out > radius) return 0;
    return (-3.0/16.0)*pow(out/radius, 5) + (5.0/8.0)*pow(out/radius, 3) -(15.0/16.0)*(out/radius) + 0.5;
}

Vector3f HRBFGenerator::grad(float x, float y, float z)
{
    //std::cout << "STARTING GRAD EVAL!!" << std::endl;
    //std::cout << "size = " << mPoints->size() << std::endl;
    fflush(stdout);
    if(mPoints->size() <= 0) return Vector3f(0,0,0);
    Vector3f p(x, y, z);
    std::cout << p << std::endl;
    fflush(stdout);
    Vector3f out(0, 0, 0);
    //std::cout << out << std::endl;
    fflush(stdout);
    for(int i = 0; i < mPoints->size()/3; i++)
    {
        int mpidx = i * 3;
        int cidx = i * 4;
        //std::cout << "idxs: " << mpidx << " " << cidx << std::endl;
        fflush(stdout);
        Vector3f vk((*mPoints)(mpidx), (*mPoints)(mpidx+1), (*mPoints)(mpidx+2));
        //std::cout << vk << std::endl;
        fflush(stdout);
        float alpha = (*unknowns)(cidx);
        //std::cout << alpha << std::endl;
        Vector3f beta((*unknowns)(cidx+1), (*unknowns)(cidx+2), (*unknowns)(cidx+3));
        //std::cout << beta << std::endl;
        Vector3f diff = p - vk;
        //std::cout << diff << std::endl;
        Vector3f grad(derivx(diff(0), diff(1), diff(2)), derivy(diff(0), diff(1), diff(2)), derivz(diff(0), diff(1), diff(2)));
        //std::cout << grad << std::endl;
        Matrix3f hess;
        hess << h00(diff(0), diff(1), diff(2)), h01(diff(0), diff(1), diff(2)), h02(diff(0), diff(1), diff(2)),
                h10(diff(0), diff(1), diff(2)), h11(diff(0), diff(1), diff(2)), h12(diff(0), diff(1), diff(2)),
                h20(diff(0), diff(1), diff(2)), h21(diff(0), diff(1), diff(2)), h22(diff(0), diff(1), diff(2));
        //std::cout << alpha << std::endl;
        //std::cout << beta(0) << "," << beta(1) << "," << beta(2) << " " << diff(0) << "," << diff(1) << "," << diff(2) << " " << grad(0) << "," << grad(1) << "," << grad(2) << std::endl;
        out += alpha*grad - hess*beta;
        //std::cout << out << std::endl;

    }
    std::cout << "GRAD OUT: " << out << std::endl;

    return out;
}

void HRBFGenerator::solve()
{
    if(coefficients->cols() > 0)
    {
        std::cout << "coefficient cols: " << coefficients->cols() << std::endl;
        *unknowns = coefficients->colPivHouseholderQr().solve(*results);
    }
    else
    {
        std::cout << "ERROR: tried to solve empty matrix!" << std::endl;
    }
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
    return coefficients;
}

VectorXf* HRBFGenerator::getUnknowns()
{
    return unknowns;
}

VectorXf* HRBFGenerator::getResults()
{
    return results;
}

int HRBFGenerator::getNumMPoints()
{
    return mPoints->size();
}

bool HRBFGenerator::getNeedRecalc()
{
    return recalc;
}

void HRBFGenerator::setRecalc(bool r)
{
    recalc = r;
}
