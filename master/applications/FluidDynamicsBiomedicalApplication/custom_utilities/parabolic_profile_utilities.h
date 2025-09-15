//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah
//                   Ruben Zorilla
//

#if !defined(KRATOS_PARABOLIC_PROFILE_UTILITIES_H )
#define  KRATOS_PARABOLIC_PROFILE_UTILITIES_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/function_parser_utility.h"

// Application includes

namespace Kratos
{

///@addtogroup FluidDynamicsBiomedicalApplication
///@{

///@name Kratos Classes
///@{

/// A set of functions to compute the Wall Shear Stress (WSS)
class KRATOS_API(FLUID_DYNAMICS_BIOMEDICAL_APPLICATION) ParabolicProfileUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParabolicProfileUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ParabolicProfileUtilities);

    using NodeType = typename ModelPart::NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Deleted default constructor
    ParabolicProfileUtilities() = delete;

    /// Deleted copy constructor.
    ParabolicProfileUtilities(ParabolicProfileUtilities const& rOther) = delete;

    /// Destructor.
    ~ParabolicProfileUtilities() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    ParabolicProfileUtilities& operator=(ParabolicProfileUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    static double CalculateInletArea(const ModelPart& rModelPart);

    static ModelPart& CreateAndFillInletAuxiliaryVolumeModelPart(ModelPart& rInletModelPart);

    static void CalculateWallParallelDistance(
        ModelPart& rWallModelPart,
        ModelPart& rFluidModelPart,
        const std::size_t WallDistanceLevels);

    static void ImposeParabolicInlet(
        ModelPart &rModelPart,
        const double MaxParabolaValue,
        const double MaxValueFactor = 1.0);

    static void ImposeParabolicInlet(
        ModelPart &rModelPart,
        const GenericFunctionUtility::Pointer rMaxParabolaValue,
        const double MaxValueFactor = 1.0);

    static void FreeParabolicInlet(ModelPart& rModelPart);

    ///@}
private:
    ///@name Private Operations
    ///@{

    template<class TInputType>
    static void ImposeParabolicProfile(
        ModelPart &rModelPart,
        const TInputType& rMaxParabolaValue,
        const double MaxValueFactor);

    template<class TInputType>
    static double GetMaxParabolaValue(
        const double Time,
        const NodeType& rNode,
        TInputType& rMaxParabolaValue);

    static double CalculateBoundingBoxCharacteristicLength(const ModelPart& rModelPart);

    ///@}
};  // Class ParabolicProfileUtilities

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PARABOLIC_PROFILE_UTILITIES_H  defined
