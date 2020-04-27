#include "hrbfmanager.cuh"

#define PROJECTIONNUM 0.35
#define MINISODIFFERENCE 0.001
#define PI 3.1415926
#define DISCONTINUITYANGLE 55
#define MAXHRBFPOINTS 50
#define MAXTHREADS 256

HRBFManager::HRBFManager()
{
    numHRBFS = 0;
    recalc = true;
    hrbfs = NULL;
}

HRBFManager::~HRBFManager()
{
    delete[] hrbfs;
    hrbf_devs->freeBufs();
    cudaFree(&hrbf_devs);
    isoVals.clear();
    numHRBFS = 0;
    cudaFree(isoVals_dev);
}

void HRBFManager::createHRBFS(int n)
{
    numHRBFS = n;
    hrbfs = new HRBFGenerator[numHRBFS];
    for(int i = 0; i < numHRBFS; i++)
    {
        hrbfs[i] = HRBFGenerator();
    }
    cudaMallocManaged((void**)&hrbf_devs, numHRBFS * sizeof(HRBFGeneratorDev));
}

void HRBFManager::clearHRBFS()
{
    delete[] hrbfs;
    hrbf_devs->freeBufs();
    cudaFree(&hrbf_devs);
    isoVals.clear();
    numHRBFS = 0;
    cudaFree(isoVals_dev);
}

bool HRBFManager::initHRBFS(float points[], int plen, float normals[], int nlen, std::vector<int> weights[], int wlen, float jointPos[], rawMat4x4 *invMats, int jlen)
{
    std::vector<float>* pts = new std::vector<float>[numHRBFS];
    std::vector<float>* norms = new std::vector<float>[numHRBFS];
    //std::cout << "plen: " << plen << std::endl;
    //std::cout << "nlen: " << nlen << std::endl;
    //std::cout << "wlen: " << wlen << std::endl;
    //std::cout << "jlen: " << jlen << std::endl;
    //std::cout << "numhrbfs: " << numHRBFS << std::endl;
    isoVals.resize(wlen);
    //sort out points and normals by weight
    for(int i = 0; i < wlen; i++)
    {
        //std::cout << "vdata: " << i << std::endl;
        //std::cout << weights[i] << std::endl;
        for(int j = 0; j < weights[i].size(); j++)
        {
            pts[weights[i][j]].push_back(points[i*3]);
            pts[weights[i][j]].push_back(points[i*3+1]);
            pts[weights[i][j]].push_back(points[i*3+2]);
            norms[weights[i][j]].push_back(normals[i*3]);
            norms[weights[i][j]].push_back(normals[i*3+1]);
            norms[weights[i][j]].push_back(normals[i*3+2]);
        }
    }
    //If there are more than 50 vertices in an HRBf, we randomly cull them until we get to 50.
    for(int i =0; i < wlen; i++)
    {
        if(weights[i].size() > MAXHRBFPOINTS)
        {
            int numToCull = weights[i].size() - MAXHRBFPOINTS;
            for(int j = 0; j < numToCull; j++)
            {
                int randPt = std::rand() % weights[i].size();
                weights[i].erase(weights[i].begin() + randPt);
            }
        }
    }
    for(int i = 0; i < numHRBFS; i++)
    {
        //std::cout << "inithrbf: " << i << std::endl;
        hrbfs[i].init(pts[i], pts[i].size(), norms[i], norms[i].size(), Eigen::Vector3f(jointPos[i*3], jointPos[i*3+1], jointPos[i*3+2]), Eigen::Vector3f(jointPos[(i+1)*3], jointPos[(i+1)*3+1], jointPos[(i+1)*3+2]));
        //std::cout << "solvehrbf:" << i << std::endl;
        hrbfs[i].solve();
        hrbf_devs[i].init(&hrbfs[i]);
    }

    int* emptyIdx = new int;
    for(int i = 0; i < wlen; i++)
    {
        //std::cout << "get isoval: " << i << std::endl;
        isoVals[i] = eval(points[i*3], points[i*3+1], points[i*3+2], invMats, emptyIdx);
    }

    if(numIsoVals == 0) cudaMallocManaged((void**)&isoVals_dev, isoVals.size() * sizeof(float));
    memcpy(isoVals_dev, isoVals.data(), isoVals.size() * sizeof(float));
    numIsoVals = isoVals.size();
    delete [] pts;
    delete [] norms;
    return true;
}

float HRBFManager::eval(float x, float y, float z, rawMat4x4* invMats, int* maxIdx)
{
    float* fs = new float[numHRBFS];
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
        //std::cout << "INVMAT: " << invMat << std::endl;
        hrbfSpace = hrbfSpace.transpose() * invMat;
        //std::cout << "EVAL " << i << " HRBFSPACE: " << hrbfSpace << std::endl;
        fs[i] = hrbfs[i].eval(hrbfSpace(0), hrbfSpace(1), hrbfSpace(2));
        //std::cout << "HOST FS[" << i << "] " << fs[i] << std::endl;
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
    delete [] fs;
    return final;
}

__device__ float HRBFManager::eval_dev(float x, float y, float z, rawMat4x4* invMats, int* maxIdx, HRBFGeneratorDev* hrbfs, int numHRBFS, int idx)
{
    float* fs = new float[numHRBFS];
    float final = -1;
    //get input position values from each hrbf.
    //if(idx == 0) printf("THREAD: %d   start eval_dev   numHRBFS: %d\n", idx, numHRBFS);
    for (int i = 0; i < numHRBFS; i++)
    {
        //if (idx == 0) printf("THREAD: %d    hrbfs[%d].mPsize: %d\n", idx, i, hrbfs[i].mPSize);
        if (hrbfs[i].mPSize <= 0) continue;

        //printf("THREAD: %d    Make hrbfSpace\n", idx);
        glm::vec4 hrbfSpace(x, y, z, 1.0f);
        //hrbfSpace << x, y, z, 1.0f;
        //if (idx == 0) printf("THREAD: %d    xyz: %f %f %f\n", idx, x, y, z);
        //if (idx == 0) printf("THREAD: %d    hrbfSpace: %f %f %f\n", idx, hrbfSpace.x, hrbfSpace.y, hrbfSpace.z);
        //printf("THREAD: %d   Make invMat\n", idx);
        glm::mat4 invMat;
        //if(idx == 44) printf("THREAD: %d   INVMAT: \n", idx);
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                invMat[j][k] = invMats[i][j][k];
                //if(idx == 44) printf("%f ", invMat[j][k]);
            }
            //if(idx == 44) printf("\n");
        }
       // printf("THREAD: %d   hrbfSpace transpose\n", idx);
        hrbfSpace = invMat * hrbfSpace;
        //if (idx == 44) printf("THREAD: %d   hrbfSpace transpose: %f %f %f\n", idx, hrbfSpace.x, hrbfSpace.y, hrbfSpace.z);
        //std::cout << "EVAL " << i << " HRBFSPACE: " << hrbfSpace << std::endl;
        fs[i] = hrbfs[i].eval(hrbfSpace.x, hrbfSpace.y, hrbfSpace.z, idx);
        //if (idx == 0) printf("THREAD: %d   fs[%d]: %f\n", idx, i, fs[i]);
    }
    //std::cout << std::endl;
    //blend the values to get the final global function value.
    //(right now I'm just taking the max).
    for (int i = 0; i < numHRBFS; i++)
    {
        if (hrbfs[i].mPSize <= 0) continue;
        if (fs[i] > final)
        {
            final = fs[i];
            *maxIdx = i;
            //std::cout << maxIdx << std::endl;
            //std::cout << *maxIdx << std::endl;
        }
    }
    //if(idx == 44) printf("THREAD: %d   EVALOUT: %f\n", idx, final);
    delete[] fs;
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

__device__ glm::vec3 HRBFManager::grad_dev(float x, float y, float z, rawMat4x4* invMats, int matIdx, HRBFGeneratorDev* hrbfs, int numHRBFS, int idx)
{
    int i = matIdx;
    glm::vec4 hrbfSpace(x, y, z, 1.0);
    glm::mat4x4 invMatGLM;
    for (int j = 0; j < 4; j++)
    {
        for (int k = 0; k < 4; k++)
        {
            invMatGLM[j][k] = invMats[i][j][k];
        }
    }
    hrbfSpace = invMatGLM * hrbfSpace;
    //printf("THREAD: %d,   START GENDEV_GRAD\n", idx);
    glm::vec3 gradGLM = hrbfs[i].grad(hrbfSpace.x, hrbfSpace.y, hrbfSpace.z);
    return gradGLM;
}

std::vector<float> HRBFManager::adjustToHRBF(float x, float y, float z, float nx, float ny, float nz, rawMat4x4 *invMats, int idx)
{
    std::cout << "IN POINT: " << x << " " << y << " " << z << std::endl;
    //fflush(stdout);
    int* maxIdx = new int;
    //std::cout << maxIdx << std::endl;
    //get iso value for current point.
    float iso = eval(x, y, z, invMats, maxIdx);
    std::cout << "EVALOUT: " << iso << std::endl;
    //fflush(stdout);
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
    //std::cout << "ORIGMAT: " << origMat << std::endl;
    worldSpace = hrbfSpace.transpose() * origMat;
    std::cout << "WORLDSPACE: " << worldSpace << std::endl;
    Eigen::Vector3f pt;
    Eigen::Vector3f oldPt;//this is in worldspace. because I suck.
    Eigen::Vector4f pt4(hrbfSpace(0), hrbfSpace(1), hrbfSpace(2), 1.0f);
    Eigen::Vector3f ptGrad;
    Eigen::Vector3f oldGrad;

    //Get the gradient at the point we're evaluating.
    //std::cout << maxIdx << std::endl;
    //std::cout << *maxIdx << std::endl;
    //fflush(stdout);
    ptGrad = grad(x, y, z, invMats, *maxIdx);
    std::cout << "PTGRAD: " << ptGrad << std::endl;
    oldGrad = ptGrad;

    //We'll need to get the difference between the gradient angles of the current iteration and the previous iteration, it should be 0 for the first iteration though.
    float gAngle = ptGrad.dot(oldGrad);
    gAngle = gAngle/(ptGrad.norm() * oldGrad.norm());
    gAngle = std::acos(gAngle);
    gAngle = gAngle/(2.0f*PI) * 360.0f;
    if(std::isnan(gAngle)) gAngle = 0;
    std::cout << "GANGLE: " << gAngle << std::endl;
    int loopCount = 0;
    //we move the point closer to it's iso value using Newton iterations. Break when we're close enough or hit an angle discontinuity
    //(The angle discontinuity suggests we've hit another HRBF and so we should stop to prevent self-intersection)
    while(std::abs(isoDiff) > MINISODIFFERENCE && gAngle < DISCONTINUITYANGLE)
    {
        std::cout << "LOOPCOUNT: " << loopCount << std::endl;
        loopCount++;
        pt << hrbfSpace(0), hrbfSpace(1), hrbfSpace(2);
        //VERTEX PROJECTION STEP.
        pt = pt + PROJECTIONNUM * (isoDiff) * (ptGrad/(ptGrad.norm()*ptGrad.norm()));

        //after adjusting the hrbfspace point, translate pt back to world space since that's what eval expects.
        pt4 << pt(0), pt(1), pt(2), 1.0f;
        std::cout << "HRBFSPACE: " << hrbfSpace << std::endl;
        oldPt << worldSpace(0), worldSpace(1), worldSpace(2);
        worldSpace = pt4.transpose() * origMat;
        std::cout << "WORLDSPACE: " << worldSpace << std::endl;

        //RELAXATION STEP (TO FIX THE BULLSHIT THAT THE VERTEX PROJECTION STEP DOES)
        float interp = std::max(0.0, (double)(1 - std::pow(std::abs(isoDiff) - 1, 4)));
        pt << worldSpace(0), worldSpace(1), worldSpace(2);
        pt = (1-interp)*pt + interp*oldPt;

        worldSpace << pt(0), pt(1), pt(2), 1.0;

        std::cout << "RELAXED WORLDSPACE: " << worldSpace << std::endl;
        hrbfSpace << worldSpace(0), worldSpace(1), worldSpace(2), 1.0f;

        iso = eval(worldSpace(0),worldSpace(1),worldSpace(2), invMats, maxIdx);
        std::cout << "EVALOUT: " << iso << std::endl;

        //use the inverse matrix identified by the eval to get the hrbfspace coordinates of the point.
        invMat <<   invMats[*maxIdx][0][0], invMats[*maxIdx][0][1], invMats[*maxIdx][0][2], invMats[*maxIdx][0][3],
                    invMats[*maxIdx][1][0], invMats[*maxIdx][1][1], invMats[*maxIdx][1][2], invMats[*maxIdx][1][3],
                    invMats[*maxIdx][2][0], invMats[*maxIdx][2][1], invMats[*maxIdx][2][2], invMats[*maxIdx][2][3],
                    invMats[*maxIdx][3][0], invMats[*maxIdx][3][1], invMats[*maxIdx][3][2], invMats[*maxIdx][3][3];
        //std::cout << "INVMAT: " << invMat << std::endl;
        hrbfSpace = hrbfSpace.transpose() * invMat;
        std::cout << "HRBFSPACE: " << hrbfSpace << std::endl;
        origMat = invMat.inverse();
        //std::cout << "ORIGMAT: " << origMat << std::endl;

        fflush(stdout);
        isoDiff = iso - isoVals[idx];
        std::cout << "ISODIFF: " << isoDiff << std::endl;
        oldGrad = ptGrad;
        std::cout << *maxIdx << std::endl;
        ptGrad = grad(worldSpace(0),worldSpace(1),worldSpace(2), invMats, *maxIdx);
        std::cout << "PTGRAD: " << ptGrad << std::endl;
        float gAngle = ptGrad.dot(oldGrad);
        std::cout << "GANGLE: " << gAngle << std::endl;
        gAngle = gAngle/(ptGrad.norm() * oldGrad.norm());
        std::cout << "PRGRADLEN: " << ptGrad.norm() << " OLDGRADLEN: " << oldGrad.norm() << " PT*OLD: " << oldGrad.norm() * ptGrad.norm() << std::endl;
        std::cout << "GANGLE: " << gAngle << std::endl;
        gAngle = acos(gAngle);
        std::cout << "GANGLE: " << gAngle << std::endl;
        gAngle = gAngle/(2.0f*PI) * 360.0f;
        std::cout << "GANGLE: " << gAngle << std::endl;
        if (isnan(gAngle)) gAngle = 0;
        std::cout << "GANGLE: " << gAngle << std::endl;
    }
    return std::vector<float> { worldSpace(0), worldSpace(1), worldSpace(2) };
}

__global__ void adjustToHRBF_dev(float* pts, int ptsSize, float* norms, int normsSize, int* indicies, int indiciesSize, float* isoVals, int numIsoVals, rawMat4x4* invMats, HRBFGeneratorDev* hrbf_devs, int numHRBFS, float* out, int outSize)
{
    //First we get the index of the point this thread is going to operate on.  We return if the thread index is larger than the number of points.
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);
    int stride = blockDim.x * gridDim.x;
    /*if (idx == 44)
    {
        for (int i = 0; i < ptsSize / 3; i++) printf("THREAD %d:   POINT %d: %f  %f  %f\n", idx, i, pts[i * 3], pts[i * 3 + 1], pts[i * 3 + 2]);
        for (int i = 0; i < normsSize / 3; i++) printf("THREAD %d:   NORMAL %d: %f  %f  %f\n", idx, i, norms[i * 3], norms[i * 3 + 1], norms[i * 3 + 2]);
        for (int i = 0; i < indiciesSize; i++) printf("THREAD %d:   INDEX %d: %d\n", idx, i, indicies[i]);
        for (int i = 0; i < numIsoVals; i++) printf("THREAD %d:   ISOVAL %d: %f\n", idx, i, isoVals[i]);
        for (int i = 0; i < numHRBFS; i++)
        {
            printf("THREAD %d:  INVMAT %d:\n", idx, i);
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    printf("%f", invMats[i][j][k]);
                }
                printf("\n");
            }
        }
    }*/
    for (int i = idx; i < indiciesSize; i += stride)
    {
        //printf("THREAD: %d START LOOP.\n", idx);
        float x = pts[idx * 3];
        float y = pts[idx * 3 + 1];
        float z = pts[idx * 3 + 2];
        //printf("THREAD: %d  PT: %f %f %f\n", idx, x, y, z);


        //std::cout << "IN POINT: " << x << " " << y << " " << z << std::endl;
        //fflush(stdout);
        int* maxIdx = new int;
        //std::cout << maxIdx << std::endl;
        //get iso value for current point.
        if(idx == 44) printf("THREAD: %d  Start eval_dev,  xyz: %f %f %f\n", idx, x, y, z);
        float iso = HRBFManager::eval_dev(x, y, z, invMats, maxIdx, hrbf_devs, numHRBFS, idx);
        if(idx == 44) printf("THREAD: %d  iso: %f\n", idx, iso);

        float isoDiff = iso - isoVals[idx];
        if(idx == 44) printf("THREAD: %d  isoDiff: %f\n", idx, isoDiff);

        //translate input point to hrbf space.
        glm::vec4 hrbfSpace(x, y, z, 1.0f);
        if(idx == 44) printf("THREAD: %d  hrbfSpace: %f %f %f\n", idx, hrbfSpace.x, hrbfSpace.y, hrbfSpace.z);
        glm::vec4 worldSpace;

        glm::mat4x4 invMatGLM;//We need to switch to glm to get the invese because Eigen still doesn't quite work with CUDA yet (and I ain't switching everything over)
        if (idx == 44) printf("THREAD: %d   INVMAT: \n", idx);
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                invMatGLM[i][j] = invMats[*maxIdx][i][j];
                if (idx == 44) printf("%f ", invMatGLM[i][j]);
            }
            if (idx == 44) printf("\n");
        }

        //glm::mat4x4 invMatGLM = glm::mat4x4(0.0);
        glm::mat4x4 doubleInvMatGLM = glm::inverse(invMatGLM);

        //std::cout << "INVMAT: " << invMat << std::endl;
        hrbfSpace = invMatGLM * hrbfSpace;
        //std::cout << "HRBFSPACE: " << hrbfSpace << std::endl;
        glm::mat4x4 origMat = doubleInvMatGLM;
        //std::cout << "ORIGMAT: " << origMat << std::endl;
        worldSpace = origMat * hrbfSpace;
        if (idx == 44) printf("THREAD: %d   worldspace: %f %f %f\n", idx, worldSpace.x, worldSpace.y, worldSpace.z);
        //std::cout << "WORLDSPACE: " << worldSpace << std::endl;
        glm::vec3 pt;
        glm::vec3 oldPt;//this is in worldspace. because I suck.
        glm::vec4 pt4(hrbfSpace.x, hrbfSpace.y, hrbfSpace.z, 1.0f);
        glm::vec3 ptGrad;
        glm::vec3 oldGrad;

        //Get the gradient at the point we're evaluating.
        //std::cout << maxIdx << std::endl;
        //std::cout << *maxIdx << std::endl;
        //fflush(stdout);
        //printf("THREAD: %d   START GRAD_DEV\n", idx);
        ptGrad = HRBFManager::grad_dev(x, y, z, invMats, *maxIdx, hrbf_devs, numHRBFS, idx);
        if (idx == 44) printf("THREAD: %d ptGrad: %f %f %f\n", idx, ptGrad.x, ptGrad.y, ptGrad.z);
        //printf("THREAD: %d   END GRAD_DEV\n", idx);
        oldGrad = ptGrad;

        //We'll need to get the difference between the gradient angles of the current iteration and the previous iteration, it should be 0 for the first iteration though.
        float gAngle = glm::dot(ptGrad, oldGrad);
        gAngle = gAngle / (glm::length(ptGrad) * glm::length(oldGrad));
        gAngle = acos(gAngle);
        gAngle = gAngle / (2.0f * PI) * 360.0f;
        if (isnan(gAngle)) gAngle = 0;
        //std::cout << "GANGLE: " << gAngle << std::endl;
        //we move the point closer to it's iso value using Newton iterations. Break when we're close enough or hit an angle discontinuity
        //(The angle discontinuity suggests we've hit another HRBF and so we should stop to prevent self-intersection)
        if (idx == 44) printf("THREAD %d:   GANGLE: %f\n", idx, gAngle);
        int loopCount = 0;
        while (abs(isoDiff) > MINISODIFFERENCE && gAngle < DISCONTINUITYANGLE)
        {
            if (idx == 44) printf("THREAD %d:   LOOPCOUNT: %d\n", idx, loopCount);
            loopCount++;
            //printf("THREAD: %d  isoDiff: %f\n", idx, isoDiff);
            pt.x = hrbfSpace.x;
            pt.y = hrbfSpace.y;
            pt.z = hrbfSpace.z;
            //VERTEX PROJECTION STEP.
            glm::vec3 mult = (ptGrad / (glm::length(ptGrad) * glm::length(ptGrad)));
            pt = pt + (float)PROJECTIONNUM * isoDiff * mult;

            //after adjusting the hrbfspace point, translate pt back to world space since that's what eval expects.
            pt4.x = pt.x;
            pt4.y = pt.y;
            pt4.z = pt.z;
            pt4.w = 1.0f;
            if (idx == 44) printf("THREAD: %d  hrbfSpace: %f %f %f\n", idx, hrbfSpace.x, hrbfSpace.y, hrbfSpace.z);
            //std::cout << "HRBFSPACE: " << hrbfSpace << std::endl;
            oldPt.x = worldSpace.x;
            oldPt.y = worldSpace.y;
            oldPt.z = worldSpace.z;
            worldSpace = origMat * pt4;
            if (idx == 44) printf("THREAD: %d   worldspace: %f %f %f\n", idx, worldSpace.x, worldSpace.y, worldSpace.z);
            //std::cout << "WORLDSPACE: " << worldSpace << std::endl;

            //RELAXATION STEP (TO FIX THE BULLSHIT THAT THE VERTEX PROJECTION STEP DOES)
            float interp = fmax(0.0, (double)(1 - pow(std::abs(isoDiff) - 1, 4)));
            pt.x = worldSpace.x;
            pt.y = worldSpace.y;
            pt.z = worldSpace.z;
            pt = (1 - interp) * pt + interp * oldPt;

            worldSpace.x = pt.x;
            worldSpace.y = pt.y;
            worldSpace.z = pt.z;
            worldSpace.w = 1.0;

            //std::cout << "RELAXED WORLDSPACE: " << worldSpace << std::endl;
            if(idx == 44) printf("THREAD: %d   worldspace relaxed: %f %f %f\n", idx, worldSpace.x, worldSpace.y, worldSpace.z);
            hrbfSpace.x = worldSpace.x;
            hrbfSpace.y = worldSpace.y;
            hrbfSpace.z = worldSpace.z;
            hrbfSpace.w = 1.0f;

            iso = HRBFManager::eval_dev(worldSpace.x, worldSpace.y, worldSpace.z, invMats, maxIdx, hrbf_devs, numHRBFS, idx);
            if (idx == 44) printf("THREAD: %d  iso: %f\n", idx, iso);
            //std::cout << "EVALOUT: " << iso << std::endl;

            //use the inverse matrix identified by the eval to get the hrbfspace coordinates of the point.
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    invMatGLM[i][j] = invMats[*maxIdx][i][j];
                    if (idx == 44) printf("%f ", invMatGLM[i][j]);
                }
                if (idx == 44) printf("\n");
            }

            doubleInvMatGLM = glm::inverse(invMatGLM);
            //std::cout << "INVMAT: " << invMat << std::endl;
            hrbfSpace = invMatGLM * hrbfSpace;
            //std::cout << "HRBFSPACE: " << hrbfSpace << std::endl;
            origMat = doubleInvMatGLM;
            //std::cout << "ORIGMAT: " << origMat << std::endl;

            //fflush(stdout);
            isoDiff = iso - isoVals[idx];
            if (idx == 44) printf("THREAD: %d isodiff: %f\n", idx, isoDiff);
            //std::cout << "ISODIFF: " << isoDiff << std::endl;
            oldGrad = ptGrad;
            //std::cout << *maxIdx << std::endl;
            ptGrad = HRBFManager::grad_dev(worldSpace.x, worldSpace.y, worldSpace.z, invMats, *maxIdx, hrbf_devs, numHRBFS, idx);
            if (idx == 44) printf("THREAD: %d ptGrad: %f %f %f\n", idx, ptGrad.x, ptGrad.y, ptGrad.z);
            float gAngle = glm::dot(ptGrad, oldGrad);
            if (idx == 44) printf("THREAD: %d   GANGLE: %f\n", idx, gAngle);
            gAngle = gAngle / (glm::length(ptGrad) * glm::length(oldGrad));
            if (idx == 44) printf("THREAD: %d   PTGRADLEN: %f OLDGRADLEN: %f PT*OLD: %f\n", idx, glm::length(ptGrad), glm::length(oldGrad), glm::length(ptGrad) * glm::length(oldGrad));
            if (idx == 44) printf("THREAD: %d   GANGLE: %f\n", idx, gAngle);
            gAngle = std::acos(gAngle);
            if (idx == 44) printf("THREAD: %d   GANGLE: %f\n", idx, gAngle);
            gAngle = gAngle / (2.0f * PI) * 360.0f;
            if (idx == 44) printf("THREAD: %d   GANGLE: %f\n", idx, gAngle);
            if (isnan(gAngle)) gAngle = 0;
            if (idx == 44) printf("THREAD: %d   GANGLE: %f\n", idx, gAngle);
        }
        if(idx == 44) printf("THREAD: %d   OUT: %f, %f, %f\n", idx, worldSpace.x, worldSpace.y, worldSpace.z);
        out[idx * 3] = worldSpace.x;
        out[idx * 3 + 1] = worldSpace.y;
        out[idx * 3 + 2] = worldSpace.z;
        
    }
    return;
}

void HRBFManager::adjustAllPoints(std::vector<float>* pts, std::vector<float>* norms, rawMat4x4* invMats, std::vector<int>* indicies, float* out, int outSize)
{
    float* pts_dev;
    float* norms_dev;
    rawMat4x4* invMats_dev;
    int* indicies_dev;
    int numPts = pts->size();
    int numNorms = norms->size();
    int numIndicies = indicies->size();
    cudaMallocManaged(&pts_dev, pts->size() * sizeof(float));
    cudaMallocManaged(&norms_dev, norms->size() * sizeof(float));
    cudaMallocManaged(&invMats_dev, numHRBFS * sizeof(rawMat4x4));
    cudaMallocManaged(&indicies_dev, indicies->size() * sizeof(float));
    memcpy(pts_dev, pts->data(), pts->size() * sizeof(float));
    memcpy(norms_dev, norms->data(), norms->size() * sizeof(float));
    for (int i = 0; i < numHRBFS; i++)
    {
        memcpy(invMats_dev[i], invMats[i], sizeof(rawMat4x4));
    }
    memcpy(indicies_dev, indicies->data(), indicies->size() * sizeof(float));

    for (int i = 0; i < numPts/3; i++) std::cout << "Point" << i << ": " << pts_dev[i * 3] << " " << pts_dev[i * 3 + 1] << pts_dev[i * 3 + 2] << std::endl;
    for (int i = 0; i < numNorms/3; i++) std::cout << "Norm" << i << ": " << norms_dev[i*3] << norms_dev[i*3+1] << norms_dev[i*3+2] << std::endl;
    for (int i = 0; i < numIndicies; i++) std::cout << "index" << i << ": " << indicies_dev[i] << std::endl;
    for (int i = 0; i < numHRBFS; i++)
    {
        std::cout << "invMat" << i << ":" << std::endl;;
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                std::cout << invMats_dev[i][j][k];
            }
            std::cout << std::endl;
        }
    }

    std::cout << "outSize:" << outSize << std::endl;
    int numBlocks = (outSize + MAXTHREADS - 1) / MAXTHREADS;
    std::cout << "numBlocks: " << numBlocks << std::endl;
    std::vector<float> hostOut = adjustToHRBF((*pts)[44 * 3], (*pts)[44 * 3 + 1], (*pts)[44 * 3 + 2], (*norms)[44 * 3], (*norms)[44 * 3 + 1], (*norms)[44 * 3 + 2], invMats, 44);
    std::cout << "FINAL HOSTADJUST: " << hostOut[0] << " " << hostOut[1] << " " << hostOut[2] << std::endl;
    adjustToHRBF_dev <<< numBlocks, MAXTHREADS >>> (pts_dev, numPts, norms_dev, numNorms, indicies_dev, numIndicies, isoVals_dev, numIsoVals, invMats_dev, hrbf_devs, numHRBFS, out, outSize);

    cudaDeviceSynchronize();

    //for (int i = 0; i < numPts / 3; i++) std::cout << "OUTPoint" << i << ": " << out[i * 3] << " " << out[i * 3 + 1] << " " << out[i * 3 + 2] << std::endl;


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
    //std::cout << "GET NEED RECALC: ";
    recalc = false;
    if(numHRBFS <= 0) recalc = true;
    for(int i = 0; i < numHRBFS; i++)
    {
        if(hrbfs[i].getNeedRecalc() == true) recalc = true;
    }
    //std::cout << recalc << std::endl;
    return recalc;
}

int HRBFManager::getNumHRBFS()
{
    return numHRBFS;
}
