// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

#ifndef SYMMETRY_BASE_H
#define SYMMETRY_BASE_H

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

// ==============================================================================

namespace Kratos
{

class SymmetryBase
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(SymmetryBase);

    typedef Node NodeType;
    typedef NodeType::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVectorType;
    typedef array_1d<double,3> array_3d;

    SymmetryBase(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters Settings)
    : mrOriginModelPart(rOriginModelPart), mrDestinationModelPart(rDestinationModelPart), mSettings(Settings)
    {
    }

    virtual ~SymmetryBase() {};

    virtual NodeVectorType& GetOriginSearchNodes() = 0;

    virtual std::vector<std::pair<array_3d, bool>> GetDestinationSearchNodes(const size_t MappingId) = 0;

    virtual void TransformationMatrix(const size_t DestinationMappingId, const size_t OriginMappingId, BoundedMatrix<double, 3, 3>& Matrix) const = 0;

    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    Parameters mSettings;

}; // Class SymmetryBase

}  // namespace Kratos.

#endif // SYMMETRY_BASE_H
