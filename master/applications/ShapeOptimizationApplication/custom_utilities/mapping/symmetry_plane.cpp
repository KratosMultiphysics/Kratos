// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "shape_optimization_application.h"
#include "symmetry_base.h"
#include "symmetry_plane.h"

// ==============================================================================

namespace Kratos
{

SymmetryPlane::SymmetryPlane(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters Settings)
: SymmetryBase(rOriginModelPart, rDestinationModelPart, Settings)
{
    mPlanePoint = mSettings["point"].GetVector();
    mPlaneNormal = mSettings["normal"].GetVector();
	KRATOS_ERROR_IF(norm_2(mPlaneNormal) < std::numeric_limits<double>::epsilon()) << "SymmetryPlane: 'normal' vector norm is 0!" << std::endl;
    mPlaneNormal /= norm_2(mPlaneNormal);

    mReflectionMatrix = IdentityMatrix(3) - (2*outer_prod(mPlaneNormal, mPlaneNormal));

    mOriginNodes.resize(mrOriginModelPart.Nodes().size());
    block_for_each(mrOriginModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
        const int mapping_id = rNode.GetValue(MAPPING_ID);
        mOriginNodes[mapping_id] = &rNode;
    });

    mDestinationNodes.resize(mrDestinationModelPart.Nodes().size());
    block_for_each(mrDestinationModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
        const int mapping_id = rNode.GetValue(MAPPING_ID);
        mDestinationNodes[mapping_id] = &rNode;
    });
}

SymmetryPlane::NodeVectorType& SymmetryPlane::GetOriginSearchNodes() {
    return mOriginNodes;
}

std::vector<std::pair<SymmetryPlane::array_3d, bool>> SymmetryPlane::GetDestinationSearchNodes(const size_t MappingId) {
    return {
        std::make_pair(mDestinationNodes[MappingId]->Coordinates(), false),
        std::make_pair(ReflectPoint(mDestinationNodes[MappingId]->Coordinates()), true),
    };
}

void SymmetryPlane::TransformationMatrix(const size_t DestinationMappingId, const size_t OriginMappingId, BoundedMatrix<double, 3, 3>& Matrix) const
{
    noalias(Matrix) = mReflectionMatrix;
    return;
}

SymmetryPlane::array_3d SymmetryPlane::ReflectPoint(const array_3d& Coords) const {
    array_3d tmp = Coords - mPlanePoint;
    tmp = prod(mReflectionMatrix, tmp);
    return tmp + mPlanePoint;
}

}  // namespace Kratos.
