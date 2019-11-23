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
    if ( bindHandle.elementCount() > 0 ) {
        for ( int i=0; i<numTransforms; ++i ) {
            transforms[i] = MFnMatrixData( bindHandle.inputValue().data() ).matrix() * transforms[i];
            bindHandle.next();
        }
    }
    MArrayDataHandle weightListHandle = block.inputArrayValue( weightList );
    if ( weightListHandle.elementCount() == 0 ) {
        // no weights - nothing to do
        return MS::kSuccess;
    }
    // Iterate through each point in the geometry.
    //
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

        // Set the final position.
        iter.setPosition( skinned );
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
