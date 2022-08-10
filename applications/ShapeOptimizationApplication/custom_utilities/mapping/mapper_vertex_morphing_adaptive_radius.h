// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_ADAPTIVE_RADIUS_H
#define MAPPER_VERTEX_MORPHING_ADAPTIVE_RADIUS_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

*/

template<class TBaseVertexMorphingMapper>
class MapperVertexMorphingAdaptiveRadius : public TBaseVertexMorphingMapper
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = TBaseVertexMorphingMapper;

    // Type definitions for better reading later
    using NodeType = Node <3>;

    /// Pointer definition of MapperVertexMorphingAdaptiveRadius
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingAdaptiveRadius);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingAdaptiveRadius(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        Parameters MapperSettings);

    /// Destructor.
    virtual ~MapperVertexMorphingAdaptiveRadius() = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    void Update() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    // Initialized by class constructor
    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    double mFilterRadiusFactor;

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateNeighbourBasedFilterRadius();

    void SmoothenNeighbourBasedFilterRadius();

    void CalculateAdaptiveVertexMorphingRadius();

    double GetVertexMorphingRadius(const NodeType& rNode) const override;

    ///@}

}; // Class MapperVertexMorphingAdaptiveRadius

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_ADAPTIVE_RADIUS_H
