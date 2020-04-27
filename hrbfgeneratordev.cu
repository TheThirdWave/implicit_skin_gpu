#include "hrbfgeneratordev.cuh"

HRBFGeneratorDev::~HRBFGeneratorDev()
{
    if (mPSize > 0) cudaFree(mPoints_dev);
    if (unknSize > 0) cudaFree(unknowns_dev);
}

__host__ void HRBFGeneratorDev::init(HRBFGenerator* hrbf)
{
    VectorXf* origMPoints = hrbf->getMPoints();
    VectorXf* origUnknowns = hrbf->getUnknowns();
    std::cout << "MPSize: " << hrbf->getNumMPoints() << std::endl;
    for (int i = 0; i < hrbf->getNumMPoints() / 3; i++)
    {
        std::cout << "HRBF MPOINTS[" << i << "]: " << origMPoints->coeff(i * 3) << " " << origMPoints->coeff(i * 3 + 1) << " " << origMPoints->coeff(i * 3 + 2) << std::endl;
    }
    std::cout << "UnknownSize: " << hrbf->getNumUnknowns() << std::endl;
    for (int i = 0; i < hrbf->getNumUnknowns() / 4; i++)
    {
        std::cout << "HRBF UNKNOWNS[" << i << "]: " << origUnknowns->coeff(i * 4) << " " << origUnknowns->coeff(i * 4 + 1) << " " << origUnknowns->coeff(i * 4 + 2) << " " << origUnknowns->coeff(i * 4 + 3) << std::endl;
    }
    float* mPointsHold = hrbf->getMPoints()->data();
    float* unknownsHold = hrbf->getUnknowns()->data();
    int numMPoints = hrbf->getNumMPoints();
    int numUnknowns = hrbf->getNumUnknowns();

    if (numUnknowns != unknSize)
    {
        if(unknSize != 0) cudaFree(unknowns_dev);
        unknSize = numUnknowns;
        cudaMallocManaged(&unknowns_dev, unknSize * sizeof(float));
    }
    if (numMPoints != mPSize)
    {
        if (mPSize != 0) cudaFree(mPoints_dev);
        mPSize = numMPoints;
        cudaMallocManaged(&mPoints_dev, mPSize * sizeof(float));
    }
    memcpy(unknowns_dev, unknownsHold, unknSize * sizeof(float));
    memcpy(mPoints_dev, mPointsHold, mPSize * sizeof(float));
    for (int i = 0; i < numMPoints/3; i++)
    {
        std::cout << "HRBF_DEV MPOINTS[" << i << "]: " << mPoints_dev[i*3] << " " << mPoints_dev[i*3+1] << " " << mPoints_dev[i*3+2] << std::endl;
    }
    for (int i = 0; i < unknSize/4; i++)
    {
        std::cout << "HRBF_DEV UNKNOWNS[" << i << "]: " << unknowns_dev[i*4] << " " << unknowns_dev[i*4+1] << " " << unknowns_dev[i*4+2] << " " << unknowns_dev[i*4+3] << std::endl;
    }
    radius = hrbf->getRadius();
}

__host__ void HRBFGeneratorDev::freeBufs()
{
    unknSize = 0;
    mPSize = 0;
    cudaFree(&unknowns_dev);
    cudaFree(&mPoints_dev);
}

__device__ float HRBFGeneratorDev::eval(float x, float y, float z, int idx)
{
    //std::cout << "STARTING EVAL!!" << std::endl;
    //std::cout << "size = " << mPoints->size() << std::endl;
    //fflush(stdout);
    glm::vec3 p(x, y, z);
    float out = 0;

    //std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
    //std::cout << "p: " << p(0) << " " << p(1) << " " << p(2) << std::endl;

    for(int i = 0; i < mPSize / 3; i++)
    {
        int mpidx = i * 3;
        int cidx = i * 4;
        glm::vec3 vk(mPoints_dev[mpidx], mPoints_dev[mpidx+1], mPoints_dev[mpidx+2]);
        //if (idx == 44) printf("THREAD %d:   vk: %f %f %f\n", idx, vk.x, vk.y, vk.z);
        //std::cout << "mpidx: " << mpidx << std::endl;
        //std::cout << "mPoints: " << mPoints_dev[mpidx] << " " << mPoints_dev[mpidx + 1] << " " << mPoints_dev[mpidx + 2] << std::endl;
        float alpha = unknowns_dev[cidx];
        //if (idx == 44) printf("THREAD %d:   alpha: %f\n", idx, alpha);
        glm::vec3 beta(unknowns_dev[cidx+1], unknowns_dev[cidx+2], unknowns_dev[cidx+3]);
        //if (idx == 44) printf("THREAD %d:   beta: %f %f %f\n", idx, beta.x, beta.y, beta.z);
        glm::vec3 diff(p.x - vk.x, p.y - vk.y, p.z - vk.z);// = p - vk;
        //if (idx == 44) printf("THREAD %d:   diff: %f %f %f\n", idx, diff.x, diff.y, diff.z);
        glm::vec3 grad(derivx_dev(diff.x, diff.y, diff.z), derivy_dev(diff.x, diff.y, diff.z), derivz_dev(diff.x, diff.y, diff.z));
        //if (idx == 44) printf("THREAD %d:   grad: %f %f %f\n", idx, grad.x, grad.y, grad.z);
        //std::cout << alpha << std::endl;
        //std::cout << "beta: " << beta(0) << "," << beta(1) << "," << beta(2) << " " << "diff: " << diff(0) << "," << diff(1) << "," << diff(2) << " " << "grad: " << grad(0) << "," << grad(1) << "," << grad(2) << std::endl;
        out += alpha*smoothfunc_dev(diff.x, diff.y, diff.z) - glm::dot(beta, grad);
        //if (idx == 44) printf("THREAD %d:   out: %f\n", idx, out);
        //std::cout << "out: " << out << std::endl;
    }

    //if (idx == 44) printf("THREAD %d:   OUT FINAL: %f\n", idx, out);
    if (out < -radius) return 1;
    if (out > radius) return 0;
    return (-3.0 / 16.0) * pow(out / radius, 5) + (5.0 / 8.0) * pow(out / radius, 3) - (15.0 / 16.0) * (out / radius) + 0.5;
}

__device__ glm::vec3 HRBFGeneratorDev::grad(float x, float y, float z)
{
    //std::cout << "STARTING GRAD EVAL!!" << std::endl;
    //std::cout << "size = " << mPoints->size() << std::endl;
    //fflush(stdout);
    if (mPSize <= 0) return glm::vec3(0, 0, 0);
    glm::vec3 p(x, y, z);
    //std::cout << p << std::endl;
    //fflush(stdout);
    glm::vec3 out(0, 0, 0);
    //std::cout << out << std::endl;
    //fflush(stdout);
    for (int i = 0; i < mPSize / 3; i++)
    {
        int mpidx = i * 3;
        int cidx = i * 4;
        //std::cout << "idxs: " << mpidx << " " << cidx << std::endl;
        //fflush(stdout);
        glm::vec3 vk(mPoints_dev[mpidx], mPoints_dev[mpidx + 1], mPoints_dev[mpidx + 2]);
        //std::cout << vk << std::endl;
        //fflush(stdout);
        float alpha = unknowns_dev[cidx];
        //std::cout << alpha << std::endl;
        glm::vec3 beta(unknowns_dev[cidx + 1], unknowns_dev[cidx + 2], unknowns_dev[cidx + 3]);
        //std::cout << beta << std::endl;
        glm::vec3 diff = p - vk;
        //std::cout << diff << std::endl;
        glm::vec3 grad(derivx_dev(diff.x, diff.y, diff.z), derivy_dev(diff.x, diff.y, diff.z), derivz_dev(diff.x, diff.y, diff.z));
        //std::cout << grad << std::endl;
        glm::mat3x3 hess;
        hess[0][0] = h00_dev(diff.x, diff.y, diff.z);
        hess[0][1] = h01_dev(diff.x, diff.y, diff.z);
        hess[0][2] = h02_dev(diff.x, diff.y, diff.z);
        hess[1][0] = h10_dev(diff.x, diff.y, diff.z);
        hess[1][1] = h11_dev(diff.x, diff.y, diff.z);
        hess[1][2] = h12_dev(diff.x, diff.y, diff.z);
        hess[2][0] = h20_dev(diff.x, diff.y, diff.z);
        hess[2][1] = h21_dev(diff.x, diff.y, diff.z);
        hess[2][2] = h22_dev(diff.x, diff.y, diff.z);
        //std::cout << alpha << std::endl;
        //std::cout << beta(0) << "," << beta(1) << "," << beta(2) << " " << diff(0) << "," << diff(1) << "," << diff(2) << " " << grad(0) << "," << grad(1) << "," << grad(2) << std::endl;
        out += alpha * grad - hess * beta;
        //std::cout << out << std::endl;

    }
    //std::cout << "GRAD OUT: " << out << std::endl;

    return out;
}

__device__ float HRBFGeneratorDev::smoothfunc_dev(float x, float y, float z)
{
    return powf(sqrtf(x * x + y * y + z * z), 3);
}

__device__ float HRBFGeneratorDev::derivx_dev(float x, float y, float z)
{
    return 3 * x * sqrtf(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::derivy_dev(float x, float y, float z)
{
    return 3 * y * sqrtf(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::derivz_dev(float x, float y, float z)
{
    return 3 * z * sqrtf(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::h00_dev(float x, float y, float z)
{
    if (x == 0 && y == 0 && z == 0) return 0;
    return 3 * (2 * x * x + y * y + z * z) / sqrt(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::h01_dev(float x, float y, float z)
{
    if (x == 0 && y == 0 && z == 0) return 0;
    return 3 * x * y / sqrt(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::h02_dev(float x, float y, float z)
{
    if (x == 0 && y == 0 && z == 0) return 0;
    return 3 * x * z / sqrt(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::h10_dev(float x, float y, float z)
{
    if (x == 0 && y == 0 && z == 0) return 0;
    return 3 * x * y / sqrt(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::h11_dev(float x, float y, float z)
{
    if (x == 0 && y == 0 && z == 0) return 0;
    return 3 * (x * x + 2 * y * y + z * z) / sqrt(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::h12_dev(float x, float y, float z)
{
    if (x == 0 && y == 0 && z == 0) return 0;
    return 3 * z * y / sqrt(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::h20_dev(float x, float y, float z)
{
    if (x == 0 && y == 0 && z == 0) return 0;
    return 3 * x * z / sqrt(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::h21_dev(float x, float y, float z)
{
    if (x == 0 && y == 0 && z == 0) return 0;
    return 3 * y * z / sqrt(x * x + y * y + z * z);
}

__device__ float HRBFGeneratorDev::h22_dev(float x, float y, float z)
{
    if (x == 0 && y == 0 && z == 0) return 0;
    return 3 * (x * x + y * y + 2 * z * z) / sqrt(x * x + y * y + z * z);
}