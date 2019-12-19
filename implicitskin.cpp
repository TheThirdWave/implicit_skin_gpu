#include "implicitskin.h"
#include <maya/MFnPlugin.h>


const MTypeId ImplicitSkin::id( 0x00080030 );

void* ImplicitSkin::creator()
{
    return new ImplicitSkin();
}

MStatus ImplicitSkin::initialize()
{
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
    rawMat4x4 inverseMatrices[numTransforms];
    if ( bindHandle.elementCount() > 0 ) {
        for ( int i=0; i<numTransforms; ++i ) {
            transforms[i] = MFnMatrixData( bindHandle.inputValue().data() ).matrix() * transforms[i];
            MFnMatrixData( bindHandle.inputValue().data() ).matrix().get(inverseMatrices[i]);
            bindHandle.next();
        }
    }
    MArrayDataHandle weightListHandle = block.inputArrayValue( weightList );
    if ( weightListHandle.elementCount() == 0 ) {
        // no weights - nothing to do
        return MS::kSuccess;
    }

    initHRBFS();
    if(hrbfs->getNeedRecalc())
    {
        cout << "START RECALC!" << endl;
        fflush(stdout);

        //create new HRBFgenerator objects if none exist
        if(hrbfs->getNumHRBFS() == 0)
        {
            cout << "Creating HRBFS!" << endl;
            hrbfs->createHRBFS(numTransforms - 1);
            fflush(stdout);
        }

        int nPoints = iter.exactCount();
        float* pts = new float[nPoints * 3];
        float* norms = new float[nPoints * 3];
        int* indicies = new int[nPoints];
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
            for ( int i=0; i<numTransforms - 1; ++i ) {
                if ( MS::kSuccess == weightsHandle.jumpToElement( i ) && weight < weightsHandle.inputValue().asDouble()) {
                     weight = weightsHandle.inputValue().asDouble();
                     indicies[idx] = i;
                }
                std::cout << weightsHandle.inputValue().asDouble();
            }
            std::cout << std::endl;
            count++;
        }
        iter.reset();

        for(int i = 0; i < numTransforms; i++)
        {
            MVector p = MTransformationMatrix(transforms[i]).getTranslation(MSpace::kWorld);
            std::cout << "pT: " << p << std::endl;
            transformPos[i * 3] = p.x;
            transformPos[i * 3 + 1] = p.y;
            transformPos[i * 3 + 2] = p.z;
        }
        std::cout << "INIT HRBFS!" << std::endl;
        hrbfs->initHRBFS(pts, nPoints * 3, norms, nPoints * 3, indicies, nPoints, transformPos, numTransforms);


    }

    // Iterate through each point in the geometry.
    //
    cout << "START ITERATION!" << endl;
    fflush(stdout);
    for ( ; !iter.isDone(); iter.next()) {
        MPoint pt = iter.position();
        MPoint skinned;
        // get the weights for this point
        MArrayDataHandle weightsHandle = weightListHandle.inputValue().child( weights );
        // compute the skinning
        for ( int i=0; i<numTransforms; ++i ) {
            if ( MS::kSuccess == weightsHandle.jumpToElement( i ) ) {
                skinned += ( pt * transforms[i] ) * weightsHandle.inputValue().asDouble();
            }
        }
        //adjust position using hrbfs to caculate self-intersections.
        std::vector<float> adjustedPt = hrbfs->adjustToHRBF(skinned.x, skinned.y, skinned.z, inverseMatrices, iter.index());
        MPoint adjPt(adjustedPt[0], adjustedPt[1], adjustedPt[2]);
        // Set the final position.
        iter.setPosition( adjPt );
        // advance the weight list handle
        weightListHandle.next();
    }
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
    std::cout << "making HRBFManager" << std::endl;
    std::cout << hrbfs << std::endl;
    if(hrbfs == NULL){

        hrbfs = new HRBFManager();
    }
}
