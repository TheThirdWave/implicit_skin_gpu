#ifndef IMPLICITSKIN_H
#define IMPLICITSKIN_H

//Have to include CUDA before Maya, and redefine some types because they both use the same names for certain things.
#include <cuda_runtime_api.h>
#define short2 CUDA_short2
#define short3 CUDA_short3
#define long2 CUDA_long2
#define long3 CUDA_long3
#define int2 CUDA_int2
#define int3 CUDA_int3
#define float2 CUDA_float2
#define float3 CUDA_float3
#define double2 CUDA_double2
#define double3 CUDA_double3
#define double4 CUDA_double4

#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MStatus.h>
#include <maya/MItMeshVertex.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MArgList.h>
#include <maya/MObject.h>
#include <maya/MGlobal.h>
#include <maya/MPxCommand.h>
#include <maya/MPxNode.h>
#include <maya/MItGeometry.h>
#include <maya/MMatrix.h>
#include <maya/MPointArray.h>
#include <maya/MTransformationMatrix.h>

#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnMesh.h>

#include <maya/MPxDeformerNode.h>

#include "hrbfgenerator.cuh"
#include "hrbfmanager.h"

#include <maya/MFnPlugin.h>
#include <maya/MTypeId.h>
#include <maya/MMatrixArray.h>
#include <maya/MStringArray.h>
#include <maya/MPxSkinCluster.h>
#include <maya/MItGeometry.h>
#include <maya/MPoint.h>
#include <maya/MFnMatrixData.h>

class ImplicitSkin : public MPxSkinCluster
{
public:

    static  void*   creator();
    static  MStatus initialize();
    // Deformation function
    //
    virtual MStatus deform(MDataBlock&    block,
                           MItGeometry&   iter,
                           const MMatrix& mat,
                           unsigned int multiIndex);
    static const MTypeId id;
    void initHRBFS();

    static MObject aHRBFRecalc;

private:
  HRBFManager* hrbfs = NULL;

  std::vector<float> getConnectedVerts(MObject meshObj, int vertIdx);
};

#endif // IMPLICITSKIN_H

