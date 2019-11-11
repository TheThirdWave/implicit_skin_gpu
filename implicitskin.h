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

class ImplicitSkin : public MPxDeformerNode {
 public:
  ImplicitSkin() {}
  virtual MStatus deform(MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex);
  static void* creator();
  static MStatus initialize();

  static MTypeId id;
  static MObject aBlendMesh;
  static MObject aBlendWeight;
  static MObject iPoint;
  static MObject oPoint;

private:
  HRBFGenerator hrbfgen;
};

#endif // IMPLICITSKIN_H

