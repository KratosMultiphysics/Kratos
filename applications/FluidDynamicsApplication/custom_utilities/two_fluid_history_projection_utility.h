//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main author:    Uxue Chasco
//                  Ruben Zorrilla
//

#if !defined(KRATOS_TWO_FLUID_HISTORY_PROJECTION_UTILITY_H_INCLUDED )
#define  KRATOS_TWO_FLUID_HISTORY_PROJECTION_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) TwoFluidHistoryProjectionUtility
{
public:

    ///@name Type Definitions
    ///@{

    using ParticleDataContainerType = std::vector<std::pair<array_1d<double,3>, array_1d<double,3>>>;

    /// Pointer definition of TwoFluidHistoryProjectionUtility
    KRATOS_CLASS_POINTER_DEFINITION(TwoFluidHistoryProjectionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    TwoFluidHistoryProjectionUtility() {};

    /// Destructor.
    ~TwoFluidHistoryProjectionUtility() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static void CalculateHistoryProjection(
        ModelPart& rModelPart,
        const bool ComputeNodalH = false);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static void FlagElementsAndNodes(
        ModelPart& rModelPart,
        const std::size_t NumberOfExtraLayers);

    static ParticleDataContainerType SeedAndConvectParticles(ModelPart& rModelPart);

    static void CalculateLagrangianVelocityInterpolation(
        ModelPart& rModelPart,
        const ParticleDataContainerType& rParticleData);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TwoFluidHistoryProjectionUtility& operator=(TwoFluidHistoryProjectionUtility const& rOther);

    /// Copy constructor.
    TwoFluidHistoryProjectionUtility(TwoFluidHistoryProjectionUtility const& rOther);

    ///@}
}; // Class TwoFluidHistoryProjectionUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const TwoFluidHistoryProjectionUtility& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TWO_FLUID_HISTORY_PROJECTION_UTILITY_H_INCLUDED  defined
