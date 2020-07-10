//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_APPLICATION_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/kratos_application.h"

// Application includes

// incompressible potential flow elements
#include "custom_elements/incompressible_potential_flow/incompressible_potential_flow_velocity_element.h"
#include "custom_elements/incompressible_potential_flow/incompressible_potential_flow_pressure_element.h"

// stabilized generic convection diffusion reaction elements
#include "custom_elements/convection_diffusion_reaction_cross_wind_stabilized_element.h"
#include "custom_elements/convection_diffusion_reaction_element.h"
#include "custom_elements/convection_diffusion_reaction_residual_based_flux_corrected_element.h"

// k-epsilon turbulence model element data
#include "custom_elements/convection_diffusion_reaction_element_data/k_epsilon/epsilon_element_data.h"
#include "custom_elements/convection_diffusion_reaction_element_data/k_epsilon/k_element_data.h"

namespace Kratos
{
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

/// Short class definition.
/** Detail class definition.
 */
class KRATOS_API(RANS_APPLICATION) KratosRANSApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosRANSApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosRANSApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosRANSApplication();

    /// Destructor.
    ~KratosRANSApplication() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

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
    std::string Info() const override
    {
        return "KratosRANSApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in my application");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());

        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    /// incompressible potential flow elements
    const IncompressiblePotentialFlowVelocityElement<2, 3> mIncompressiblePotentialFlowVelocity2D;
    const IncompressiblePotentialFlowVelocityElement<3, 4> mIncompressiblePotentialFlowVelocity3D;
    const IncompressiblePotentialFlowPressureElement<2, 3> mIncompressiblePotentialFlowPressure2D;
    const IncompressiblePotentialFlowPressureElement<3, 4> mIncompressiblePotentialFlowPressure3D;

    /// k-epsilon turbulence model elements
    /// Algebraic flux correction based elements
    const ConvectionDiffusionReactionElement<2, 3, KEpsilonElementData::KElementData<2>> mRansKEpsilonKAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, KEpsilonElementData::KElementData<3>> mRansKEpsilonKAFC3D;

    const ConvectionDiffusionReactionElement<2, 3, KEpsilonElementData::EpsilonElementData<2>> mRansKEpsilonEpsilonAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, KEpsilonElementData::EpsilonElementData<3>> mRansKEpsilonEpsilonAFC3D;

    /// Residual based flux corrected elements
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KEpsilonElementData::KElementData<2>> mRansKEpsilonKRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KEpsilonElementData::KElementData<3>> mRansKEpsilonKRFC3D;

    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KEpsilonElementData::EpsilonElementData<2>> mRansKEpsilonEpsilonRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KEpsilonElementData::EpsilonElementData<3>> mRansKEpsilonEpsilonRFC3D;

    /// Cross wind stabilization based elements
    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, KEpsilonElementData::KElementData<2>> mRansKEpsilonKCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, KEpsilonElementData::KElementData<3>> mRansKEpsilonKCWD3D;

    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, KEpsilonElementData::EpsilonElementData<2>> mRansKEpsilonEpsilonCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, KEpsilonElementData::EpsilonElementData<3>> mRansKEpsilonEpsilonCWD3D;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
    KratosRANSApplication& operator=(KratosRANSApplication const& rOther);

    /// Copy constructor.
    KratosRANSApplication(KratosRANSApplication const& rOther);

    ///@}

}; // Class KratosRANSApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_APPLICATION_H_INCLUDED  defined
