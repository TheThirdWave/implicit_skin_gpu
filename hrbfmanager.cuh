#ifndef HRBFMANAGER_H
#define HRBFMANAGER_H

#include <../eigen-3.3.7/Eigen/Dense>
#include <vector>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime_api.h>
#define GLM_FORCE_CUDA
#include <../glm/glm/glm.hpp>
#include <../glm/glm/gtc/type_ptr.hpp>

#include "hrbfgenerator.cuh"
#include "hrbfgeneratordev.cuh"

class HRBFManager
{
private:
    bool recalc;
    int numHRBFS;
    HRBFGenerator* hrbfs;
    HRBFGeneratorDev* hrbf_devs;
    std::vector<float> isoVals;
    float* isoVals_dev;
    int numIsoVals = 0;


public:
    HRBFManager();
    ~HRBFManager();

    void createHRBFS(int n);
    void clearHRBFS();
    bool initHRBFS(float points[], int plen, float normals[], int nlen, std::vector<int> weights[], int wlen, float jointPos[], rawMat4x4 *invMats, int jlen);
    float eval(float x, float y, float z, rawMat4x4 *invMats, int *maxIdx);
    __device__ static float eval_dev(float x, float y, float z, rawMat4x4* invMats, int* maxIdx, HRBFGeneratorDev* hrbfs, int numHRBFS, int idx);
    Eigen::Vector3f grad(float x, float y, float z, rawMat4x4 *invMats, int matIdx);
    __device__ static glm::vec3 grad_dev(float x, float y, float z, rawMat4x4* invMats, int matIdx, HRBFGeneratorDev* hrbfs, int numHRBFS, int idx);
    std::vector<float> adjustToHRBF(float x, float y, float z, float nx, float ny, float nz, rawMat4x4 *invMats, int idx);
    void adjustAllPoints(std::vector<float>* pts, std::vector<float>* norms, rawMat4x4* invMats, std::vector<int>* indicies, float* out, int outSize);
    void setNeedRecalc(bool r);
    bool getNeedRecalc();
    int getNumHRBFS();
};

#endif // HRBFMANAGER_H
