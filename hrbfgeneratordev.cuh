#ifndef HRBFGENERATORDEV
#define HRBFGENERATORDEV

#include <../eigen-3.3.7/Eigen/dense>
#include <cuda.h>
#define GLM_FORCE_CUDA
#include <../glm/glm/glm.hpp>
#include <../glm/glm/gtc/type_ptr.hpp>>
#include <vector>
//#include <cuda.h>
//#include <cuda_runtime_api.h>
//#include <device_launch_parameters.h>
#include <iostream>
#include "SerialFunctions.h"
#include "Vector3fDev.cuh"
#include "HRBFGenerator.cuh"

using Eigen::Vector3f;
using Eigen::Matrix3f;

class HRBFGeneratorDev 
{
private:

public:
    float* unknowns_dev;
    float* mPoints_dev;
    int unknSize = 0;
    int mPSize = 0;
    float radius = -1;

    ~HRBFGeneratorDev();

    __host__ void init(HRBFGenerator* hrbf);
    __host__ void freeBufs();

    __device__ float eval(float x, float y, float z, int idx);
    __device__ glm::vec3 grad(float x, float y, float z);

    static __device__ float smoothfunc_dev(float x, float y, float z);
    static __device__ float derivx_dev(float x, float y, float z);
    static __device__ float derivy_dev(float x, float y, float z);
    static __device__ float derivz_dev(float x, float y, float z);
    static __device__ float h00_dev(float x, float y, float z);
    static __device__ float h01_dev(float x, float y, float z);
    static __device__ float h02_dev(float x, float y, float z);
    static __device__ float h10_dev(float x, float y, float z);
    static __device__ float h11_dev(float x, float y, float z);
    static __device__ float h12_dev(float x, float y, float z);
    static __device__ float h20_dev(float x, float y, float z);
    static __device__ float h21_dev(float x, float y, float z);
    static __device__ float h22_dev(float x, float y, float z);

};

#endif //HRBFGENERATORDEV
