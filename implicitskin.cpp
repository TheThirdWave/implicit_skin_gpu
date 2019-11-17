#include "implicitskin.h"
#include <maya/MFnPlugin.h>

MTypeId ImplicitSkin::id(0x00000451);
MObject ImplicitSkin::aBlendMesh;
MObject ImplicitSkin::aBlendWeight;
MObject ImplicitSkin::oPoint;
MObject ImplicitSkin::iPoint;

void* ImplicitSkin::creator() { return new ImplicitSkin; }


MStatus ImplicitSkin::deform(MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex) {

  MStatus status;

  cout << "TEST OUT!" << endl;
  fflush(stdout);

  //get the envelope and blend weight
  float env = data.inputValue(envelope).asFloat();
  float blendWeight = data.inputValue(aBlendWeight).asFloat();
  blendWeight *= env;

  //get the blend mesh
  MObject oBlendMesh = data.inputValue(aBlendMesh).asMesh();
  if(oBlendMesh.isNull()){
      //No blend mesh attached so exit node
      return MS::kSuccess;
  }



  //Get the blend points
  MFnMesh fnBlendMesh(oBlendMesh, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MPointArray blendPoints;
  fnBlendMesh.getPoints(blendPoints);

  MPoint iPt(data.inputValue(iPoint).asFloat3());

  MDataHandle oPtHandle = data.outputValue(oPoint);


  if(hrbfgen.getNeedRecalc())
  {
      cout << "START RECALC!" << endl;
      fflush(stdout);
      int nPoints = itGeo.exactCount();
      float* pts = new float[nPoints * 3];
      float* norms = new float[nPoints * 3];

      MPoint pt;
      int idx;
      for (; !itGeo.isDone(); itGeo.next()){
          //Get the input point
          pt = itGeo.position();
          idx = itGeo.index();
          pts[idx*3] = pt.x;
          norms[idx*3] = pt.x;
          pts[idx*3+1] = pt.y;
          norms[idx*3+1] = pt.y;
          pts[idx*3+2] = pt.z;
          norms[idx*3+2] = pt.z;
      }
      hrbfgen.init(pts, nPoints * 3, norms, nPoints * 3);
      hrbfgen.solve();


  }

  float outNum = hrbfgen.eval(iPt.x, iPt.y, iPt.z);


  cout << iPt.x << ", " << iPt.y << ", " << iPt.z << endl;
  cout << outNum << endl;
  fflush(stdout);
  oPtHandle.set(outNum);

  MPoint pt;
  float w = 0.0f;
  for (; !itGeo.isDone(); itGeo.next()){
      //Get the input point
      pt = itGeo.position();
      //Get the painted weight value
      w = weightValue(data, mIndex, itGeo.index());
      //Perform the deformation
      pt = pt + (blendPoints[itGeo.index()] - pt) * blendWeight * w;
      //Set the new output point
      itGeo.setPosition(pt);
  }

  return MS::kSuccess;
}

MStatus ImplicitSkin::initialize(){
    MFnTypedAttribute tAttr;
    MFnNumericAttribute nAttr;

    aBlendMesh = tAttr.create("blendMesh", "blendMesh", MFnNumericData::kMesh);
    addAttribute(aBlendMesh);
    attributeAffects(aBlendMesh, outputGeom);

    aBlendWeight = nAttr.create("blendWeight", "bw", MFnNumericData::kFloat);
    nAttr.setKeyable(true);
    addAttribute(aBlendWeight);
    attributeAffects(aBlendWeight, outputGeom);

    oPoint = nAttr.create("oPoint", "opt", MFnNumericData::kFloat);
    nAttr.setWritable(false);
    nAttr.setStorable(false);
    addAttribute(oPoint);

    iPoint = nAttr.create("iPoint", "ipt", MFnNumericData::k3Float);
    nAttr.setKeyable(true);
    addAttribute(iPoint);
    attributeAffects(iPoint, oPoint);

    //Make the deformer weights paintable
    MGlobal::executeCommand("makePaintable -attrType multiFloat -sm deformer blendNode weights;");

    return MS::kSuccess;
}

MStatus initializePlugin(MObject obj) {
  MFnPlugin plugin(obj, "C.J. Dopheide", "1.0", "Any");
  MStatus status = plugin.registerNode("implicitSkinning", ImplicitSkin::id, ImplicitSkin::creator, ImplicitSkin::initialize, MPxNode::kDeformerNode);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  return status;
}

MStatus uninitializePlugin(MObject obj) {
  MFnPlugin plugin(obj);
  MStatus status = plugin.deregisterNode(ImplicitSkin::id);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  return status;
}
