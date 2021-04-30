// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

#ifndef SYMMETRY_REVOLUTION_H
#define SYMMETRY_REVOLUTION_H

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
#include "symmetry_base.h"

// ==============================================================================

namespace Kratos
{

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) SymmetryRevolution : public SymmetryBase
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(SymmetryRevolution);

    SymmetryRevolution(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters Settings);

    NodeVectorType& GetOriginSearchNodes() override;

    std::vector<std::pair<array_3d, bool>> GetDestinationSearchNodes(const size_t MappingId) override;

    void TransformationMatrix(const size_t DestinationMappingId, const size_t OriginMappingId, BoundedMatrix<double, 3, 3>& Matrix) const override;

    NodeTypePointer GetTransformedNode(const NodeType& rNode);

    array_3d mPoint;
    array_3d mAxis;
    array_3d mPlaneVector1; // vector in the plane orthogonal to the mAxis

    NodeVectorType mOriginNodes;
    NodeVectorType mDestinationNodes;
    NodeVectorType mTransformedOriginNodes;
    NodeVectorType mTransformedDestinationNodes;


}; // Class SymmetryRevolution

}  // namespace Kratos.

#endif // SYMMETRY_REVOLUTION_H
