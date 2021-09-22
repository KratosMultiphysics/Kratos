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

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "symmetry_base.h"

// ==============================================================================

namespace Kratos
{

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) SymmetryPlane : public SymmetryBase
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(SymmetryPlane);

    SymmetryPlane(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters Settings);

    NodeVectorType& GetOriginSearchNodes() override;

    std::vector<std::pair<array_3d, bool>> GetDestinationSearchNodes(const size_t MappingId) override;

    void TransformationMatrix(const size_t DestinationMappingId, const size_t OriginMappingId, BoundedMatrix<double, 3, 3>& Matrix) const override;

    array_3d ReflectPoint(const array_3d& Coords) const;

    array_3d mPlanePoint;
    array_3d mPlaneNormal;

    NodeVectorType mOriginNodes;
    NodeVectorType mDestinationNodes;
    Matrix mReflectionMatrix;
}; // Class SymmetryPlane

}  // namespace Kratos.

#endif // SYMMETRY_PLANE_H
