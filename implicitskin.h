#ifndef IMPLICITSKIN_H
#define IMPLICITSKIN_H

#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MStatus.h>
#include <maya/MArgList.h>
#include <maya/MObject.h>
#include <maya/MGlobal.h>
#include <maya/MPxCommand.h>
#include <maya/MPxNode.h>
#include <maya/MItGeometry.h>
#include <maya/MMatrix.h>
#include <maya/MPointArray.h>

#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnMesh.h>

#include <maya/MPxDeformerNode.h>

#include "hrbfgenerator.h"
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

private:
  HRBFManager hrbfs;
};

#endif // IMPLICITSKIN_H

