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

// stabilized generic convection diffusion reaction elements
#include "custom_elements/convection_diffusion_reaction_cross_wind_stabilized_element.h"
#include "custom_elements/convection_diffusion_reaction_element.h"
#include "custom_elements/convection_diffusion_reaction_residual_based_flux_corrected_element.h"

// k-epsilon turbulence model element data
#include "custom_elements/element_data/evm_k_epsilon_high_re/epsilon_element_data.h"
#include "custom_elements/element_data/evm_k_epsilon_high_re/k_element_data.h"

// k-omega turbulence model element data
#include "custom_elements/element_data/evm_k_omega/k_element_data.h"
#include "custom_elements/element_data/evm_k_omega/omega_element_data.h"

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

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    /// k-epsilon turbulence model elements
    /// Algebraic flux correction based elements
    const ConvectionDiffusionReactionElement<2, 3, EvmKEpsilonHighReElementData::KElementData<2>> mRansEvmKEpsilonHighReKAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, EvmKEpsilonHighReElementData::KElementData<3>> mRansEvmKEpsilonHighReKAFC3D;

    const ConvectionDiffusionReactionElement<2, 3, EvmKEpsilonHighReElementData::EpsilonElementData<2>> mRansEvmKEpsilonHighReEpsilonAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, EvmKEpsilonHighReElementData::EpsilonElementData<3>> mRansEvmKEpsilonHighReEpsilonAFC3D;

    /// Residual based flux corrected elements
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, EvmKEpsilonHighReElementData::KElementData<2>> mRansEvmKEpsilonHighReKRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, EvmKEpsilonHighReElementData::KElementData<3>> mRansEvmKEpsilonHighReKRFC3D;

    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, EvmKEpsilonHighReElementData::EpsilonElementData<2>> mRansEvmKEpsilonHighReEpsilonRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, EvmKEpsilonHighReElementData::EpsilonElementData<3>> mRansEvmKEpsilonHighReEpsilonRFC3D;

    /// Cross wind stabilization based elements
    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, EvmKEpsilonHighReElementData::KElementData<2>> mRansEvmKEpsilonHighReKCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, EvmKEpsilonHighReElementData::KElementData<3>> mRansEvmKEpsilonHighReKCWD3D;

    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, EvmKEpsilonHighReElementData::EpsilonElementData<2>> mRansEvmKEpsilonHighReEpsilonCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, EvmKEpsilonHighReElementData::EpsilonElementData<3>> mRansEvmKEpsilonHighReEpsilonCWD3D;

    /// k-omega turbulence model elements
    /// Algebraic flux correction based elements
    const ConvectionDiffusionReactionElement<2, 3, EvmKOmegaElementData::KElementData<2>> mRansEvmKOmegaKAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, EvmKOmegaElementData::KElementData<3>> mRansEvmKOmegaKAFC3D;

    const ConvectionDiffusionReactionElement<2, 3, EvmKOmegaElementData::OmegaElementData<2>> mRansEvmKOmegaOmegaAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, EvmKOmegaElementData::OmegaElementData<3>> mRansEvmKOmegaOmegaAFC3D;

    /// Residual based flux corrected elements
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, EvmKOmegaElementData::KElementData<2>> mRansEvmKOmegaKRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, EvmKOmegaElementData::KElementData<3>> mRansEvmKOmegaKRFC3D;

    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, EvmKOmegaElementData::OmegaElementData<2>> mRansEvmKOmegaOmegaRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, EvmKOmegaElementData::OmegaElementData<3>> mRansEvmKOmegaOmegaRFC3D;

    /// Cross wind stabilization based elements
    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, EvmKOmegaElementData::KElementData<2>> mRansEvmKOmegaKCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, EvmKOmegaElementData::KElementData<3>> mRansEvmKOmegaKCWD3D;

    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, EvmKOmegaElementData::OmegaElementData<2>> mRansEvmKOmegaOmegaCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, EvmKOmegaElementData::OmegaElementData<3>> mRansEvmKOmegaOmegaCWD3D;


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
