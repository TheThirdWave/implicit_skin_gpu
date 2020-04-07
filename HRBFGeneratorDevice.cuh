/*

__device__ float derivx_dev(float x, float y, float z);
__device__ float derivy_dev(float x, float y, float z);
__device__ float derivz_dev(float x, float y, float z);
__device__ float smoothfunc_dev(float x, float y, float z);
__global__ void getPointEvals(float* unknowns, float* mPoints, float* outs);*/