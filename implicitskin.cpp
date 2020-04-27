#include "implicitskin.h"
#include <maya/MFnPlugin.h>


const MTypeId ImplicitSkin::id( 0x00080030 );
MObject ImplicitSkin::aHRBFRecalc;

void* ImplicitSkin::creator()
{
    return new ImplicitSkin();
}

MStatus ImplicitSkin::initialize()
{
    MFnNumericAttribute nHRBFRecalc;

    aHRBFRecalc = nHRBFRecalc.create("Recalc", "recalc", MFnNumericData::kBoolean);
    nHRBFRecalc.setWritable(true);
    nHRBFRecalc.setKeyable(true);
    addAttribute(aHRBFRecalc);

    return MStatus::kSuccess;
}

MStatus ImplicitSkin::deform( MDataBlock& block,
                      MItGeometry& iter,
                      const MMatrix& /*m*/,
                      unsigned int multiIndex)
//
// Method: deform
//
// Description:   Deforms the point with a simple smooth skinning algorithm
//
// Arguments:
//   block      : the datablock of the node
//   iter       : an iterator for the geometry to be deformed
//   m          : matrix to transform the point into world space
//   multiIndex : the index of the geometry that we are deforming
//
//
{
    MStatus returnStatus;


    // get input mesh.
    MArrayDataHandle meshHandle = block.outputArrayValue( input );
    MObject meshObject = meshHandle.outputValue().child( inputGeom ).asMesh();
    if(meshObject.isNull())
    {
        std::cout << "meshObject is Null! (no mesh attached to deformer)" << endl;
    }
    std::cout << meshObject.apiTypeStr() << endl;
    MFnMesh meshData;
    meshData.setObject(meshObject);

    // get the influence transforms
    //

    MArrayDataHandle transformsHandle = block.inputArrayValue( matrix );
    int numTransforms = transformsHandle.elementCount();
    if ( numTransforms == 0 ) {
        return MS::kSuccess;
    }
    MMatrixArray transforms;
    for ( int i=0; i<numTransforms; ++i ) {
        transforms.append( MFnMatrixData( transformsHandle.inputValue().data() ).matrix() );
        transformsHandle.next();
    }
    MArrayDataHandle bindHandle = block.inputArrayValue( bindPreMatrix );
    rawMat4x4* inverseMatrices = new rawMat4x4[numTransforms];
    MMatrixArray invTransforms;
    if ( bindHandle.elementCount() > 0 ) {
        for ( int i=0; i<numTransforms; ++i ) {
            transforms[i] = MFnMatrixData( bindHandle.inputValue().data() ).matrix() * transforms[i];
            invTransforms.append(MFnMatrixData( bindHandle.inputValue().data() ).matrix());
            transforms[i].inverse().get(inverseMatrices[i]);
            //std::cout << "INVMAT" << i << ": " << std::endl;
            for(int j = 0; j < 4; j++)
            {
                for(int k = 0; k < 4; k++)
                {
                    std::cout << inverseMatrices[i][j][k] << " ";
                }
                std::cout << std::endl;
            }
            bindHandle.next();
        }
    }

    MArrayDataHandle weightListHandle = block.inputArrayValue( weightList );
    if ( weightListHandle.elementCount() == 0 ) {
        // no weights - nothing to do
        return MS::kSuccess;
    }

    if(block.inputValue(aHRBFRecalc).asBool() == true)
    {
        MDataHandle reset = block.inputValue(aHRBFRecalc);
        reset.setBool(false);
        reset.setClean();
        hrbfs->clearHRBFS();
        hrbfs->setNeedRecalc(true);
    }

    int nPoints = iter.exactCount();
    initHRBFS();
    if(hrbfs->getNeedRecalc())
    {
        cout << "START RECALC!" << endl;
        //fflush(stdout);

        //create new HRBFgenerator objects if none exist
        if(hrbfs->getNumHRBFS() == 0)
        {
            cout << "Creating HRBFS!" << endl;
            hrbfs->createHRBFS(numTransforms - 1);
            //fflush(stdout);
        }

        float* pts = new float[nPoints * 3];
        float* norms = new float[nPoints * 3];
        std::vector<int>* indicies = new std::vector<int>[nPoints];
        float* transformPos = new float[numTransforms * 3];
        MPoint pt;
        MVector norm;
        int idx;
        int count = 0;

        //get vertex position, normals, and joint transform data into the hrbf manager.
        for ( ; !iter.isDone(); iter.next()) {
            pt = iter.position();
            idx = iter.index();
            meshData.getVertexNormal(idx, false, norm, MSpace::kObject);
            pts[idx*3] = pt.x;
            pts[idx*3+1] = pt.y;
            pts[idx*3+2] = pt.z;
            norms[idx*3] = norm.x;
            norms[idx*3+1] = norm.y;
            norms[idx*3+2] = norm.z;
            //sort out which hrbf the vertex shall be associated with using it's weight data (the joint it's most weighted to shall be the hrbf it's associated with)
            // get the weights for this point
            MArrayDataHandle weightsHandle = weightListHandle.inputValue().child( weights );
            // compute the skinning
            double weight = -1;
            //std::cout << idx << " " << weightsHandle.elementCount() << ": ";
            for ( int i=0; i<numTransforms; ++i ) {
                if ( MS::kSuccess == weightsHandle.jumpToElement( i ) ) {

                    double holdWeight = weightsHandle.inputValue().asDouble();

                    if(i < numTransforms-1)//(weight < holdWeight && i < (numTransforms-1))
                    {
                        weight = holdWeight;
                        indicies[idx].push_back(i);
                    }
                    //std::cout << i << ": " << weightsHandle.inputValue().asDouble() << " ";
                }
            }
            weightListHandle.next();
            count++;
        }
        weightListHandle.jumpToElement(0);
        iter.reset();

        for(int i = 0; i < numTransforms; i++)
        {
            MVector p = MTransformationMatrix(invTransforms[i].inverse()).getTranslation(MSpace::kWorld);
            //std::cout << "pT: " << p << std::endl;
            transformPos[i * 3] = p.x;
            transformPos[i * 3 + 1] = p.y;
            transformPos[i * 3 + 2] = p.z;
        }
        std::cout << "INIT HRBFS!" << std::endl;
        hrbfs->initHRBFS(pts, nPoints * 3, norms, nPoints * 3, indicies, nPoints, transformPos, inverseMatrices, numTransforms);


    }

    std::vector<float> points;
    std::vector<float> normals;
    std::vector<int> indicies;
    points.resize(nPoints * 3);
    normals.resize(nPoints * 3);
    indicies.resize(nPoints);

    // Iterate through each point in the geometry.
    //
    MStatus status;
    int count = 0;
    cout << "START ITERATION!" << endl;
    //fflush(stdout);
    for ( ; !iter.isDone(); iter.next()) {
        MPoint pt = iter.position();
        //std::cout << "PTIDX: " << iter.index() << std::endl;
        MObject item = iter.currentItem();
        //std::cout << "curItem: " << item.apiTypeStr() << std::endl;
        MFnSingleIndexedComponent fnCom(item, &status);
        int numComs;
        fnCom.getCompleteData(numComs);
        //std::cout << "num components: " << numComs << std::endl;
        MIntArray components;
        fnCom.getElements(components);
        //std::cout << "Components: ";
        //std::cout << std::endl;
        MVector norm;
        meshData.getVertexNormal(iter.index(), false, norm, MSpace::kObject);

        MPoint skinned;
        // get the weights for this point
        MArrayDataHandle weightsHandle = weightListHandle.inputValue().child( weights );
        // compute the skinning
        for ( int i=0; i<numTransforms; ++i ) {
            if ( MS::kSuccess == weightsHandle.jumpToElement( i ) ) {
                skinned += ( pt * transforms[i] ) * weightsHandle.inputValue().asDouble();
                //std::cout << "tidx: " << i << std::endl;
            }
        }
        int pIdx = count * 3;
        //std::cout << "saving point" << std::endl;
        points[pIdx] = skinned.x;
        points[pIdx + 1] = skinned.y;
        points[pIdx + 2] = skinned.z;
        //std::cout << "saving normal" << std::endl;
        normals[pIdx] = norm.x;
        normals[pIdx + 1] = norm.y;
        normals[pIdx + 2] = norm.z;
        //std::cout << "saving index" << std::endl;
        indicies[count] = iter.index();

        weightListHandle.next();
        count++;
    }
    //std::cout << "iter reset" << std::endl;
    iter.reset();
    count = 0;
    float* adjPts;
    cudaMallocManaged((void**)&adjPts, (nPoints * 3) * sizeof(float));
    std::cout << "ADJUST ALL POINTS" << std::endl;
    hrbfs->adjustAllPoints(&points, &normals, inverseMatrices, &indicies, adjPts, nPoints * 3);
    for (; !iter.isDone(); iter.next())
    {
        int pIdx = count * 3;
        MPoint adjPt(adjPts[pIdx], adjPts[pIdx + 1], adjPts[pIdx + 2]);
        //std::cout << "OUT POINT " << pIdx << ": " << adjPts[pIdx] << " " << adjPts[pIdx + 1] << " " << adjPts[pIdx + 2] << std::endl;
        //std::cout << "ADJPOINT: " << adjPt << std::endl;
        iter.setPosition(adjPt);
        count++;
    }
    cudaFree(adjPts);
    delete [] inverseMatrices;
    return returnStatus;
}
// standard initialization procedures
//
MStatus initializePlugin( MObject obj )
{

    MStatus result;
    MFnPlugin plugin( obj, PLUGIN_COMPANY, "3.0", "Any");
    result = plugin.registerNode(
        "implicitSkinning" ,
        ImplicitSkin::id ,
        &ImplicitSkin::creator ,
        &ImplicitSkin::initialize ,
        MPxNode::kSkinCluster
        );
    return result;
}
MStatus uninitializePlugin( MObject obj )
{
    MStatus result;
    MFnPlugin plugin( obj );
    result = plugin.deregisterNode( ImplicitSkin::id );
    return result;
}

void ImplicitSkin::initHRBFS()
{
    if(hrbfs == NULL){
        std::cout << "making HRBFManager" << std::endl;
        std::cout << hrbfs << std::endl;
        hrbfs = new HRBFManager();
    }
}

std::vector<float> ImplicitSkin::getConnectedVerts(MObject meshObj, int vertIdx)
{
    MStatus status;
    MItMeshVertex vertIt(meshObj, &status);
    if(status == MS::kFailure)
    {
        std::cout << "vertIt ERROR: " << status.errorString() << std::endl;
    }
    int oldIdx;
    vertIt.setIndex(vertIdx, oldIdx);
    MIntArray connectedIds;
    vertIt.getConnectedVertices(connectedIds);
    //std::cout << "VERT NUM CONNECTED: " << connectedIds.length() << std::endl;
    std::vector<float> points;
    points.resize(connectedIds.length()*3);
    for(int i = 0; i < connectedIds.length(); i++)
    {
        vertIt.setIndex(connectedIds[i], oldIdx);
        MPoint curPt = vertIt.position();
        points[i*3] = curPt.x;
        points[i*3+1] = curPt.y;
        points[i*3+2] = curPt.z;
    }
    return points;
}
