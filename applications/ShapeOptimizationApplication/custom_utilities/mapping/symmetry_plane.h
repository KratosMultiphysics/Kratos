// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

#ifndef SYMMETRY_PLANE_H
#define SYMMETRY_PLANE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "symmetry_base.h"

// ==============================================================================

namespace Kratos
{

class SymmetryPlane : public SymmetryBase
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(SymmetryPlane);

    SymmetryPlane(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters Settings)
    : SymmetryBase(rOriginModelPart, rDestinationModelPart, Settings)
    {
        mPlanePoint = mSettings["point"].GetVector();
        mPlaneNormal = mSettings["normal"].GetVector();
        mPlaneNormal /= norm_2(mPlaneNormal);

        mReflectionMatrix = IdentityMatrix(3) - (2*outer_prod(mPlaneNormal, mPlaneNormal));

        mOriginNodes.resize(mrOriginModelPart.Nodes().size());
        for (auto& r_node_i : mrOriginModelPart.Nodes()){
            const int mapping_id = r_node_i.GetValue(MAPPING_ID);
            mOriginNodes[mapping_id] = &r_node_i;
        }

        mDestinationNodes.resize(mrDestinationModelPart.Nodes().size());
        for (auto& r_node_i : mrDestinationModelPart.Nodes()){
            const int mapping_id = r_node_i.GetValue(MAPPING_ID);
            mDestinationNodes[mapping_id] = &r_node_i;
        }
    }

    NodeVector& GetOriginSearchNodes() override {
        return mOriginNodes;
    }

    std::vector<std::pair<array_3d, bool>> GetDestinationSearchNodes(const size_t MappingId) override {
        return {
            std::make_pair(mDestinationNodes[MappingId]->Coordinates(), false),
            std::make_pair(ReflectPoint(mDestinationNodes[MappingId]->Coordinates()), true),
        };
    }

    void TransformationMatrix(const size_t DestinationMappingId, const size_t OriginMappingId, BoundedMatrix<double, 3, 3>& Matrix) const override
    {
        noalias(Matrix) = mReflectionMatrix;
        return;
    }

    array_3d ReflectPoint(const array_3d& Coords){
        array_3d tmp = Coords - mPlanePoint;
        tmp = prod(mReflectionMatrix, tmp);
        return tmp + mPlanePoint;
    }

    array_3d mPlanePoint;
    array_3d mPlaneNormal;

    NodeVector mOriginNodes;
    NodeVector mDestinationNodes;
    Matrix mReflectionMatrix;
}; // Class SymmetryPlane

}  // namespace Kratos.

#endif // SYMMETRY_PLANE_H
