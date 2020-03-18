//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//
//  Main authors:    Riccardo Rossi, Inigo Lopez and Marc Nunez
//

#if !defined(KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_H_INCLUDED )
#define  KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "custom_elements/compressible_potential_flow_element.h"
#include "custom_elements/incompressible_potential_flow_element.h"
#include "custom_elements/compressible_perturbation_potential_flow_element.h"
#include "custom_elements/incompressible_perturbation_potential_flow_element.h"
#include "custom_elements/embedded_incompressible_potential_flow_element.h"
#include "custom_elements/embedded_compressible_potential_flow_element.h"
#include "custom_conditions/potential_wall_condition.h"
#include "custom_elements/adjoint_analytical_incompressible_potential_flow_element.h"
#include "custom_elements/adjoint_finite_difference_potential_flow_element.h"
#include "custom_conditions/adjoint_potential_wall_condition.h"
namespace Kratos {

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) KratosCompressiblePotentialFlowApplication : public KratosApplication {
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of KratosCompressiblePotentialFlowApplication
	KRATOS_CLASS_POINTER_DEFINITION(KratosCompressiblePotentialFlowApplication);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosCompressiblePotentialFlowApplication();

	/// Destructor.
	~KratosCompressiblePotentialFlowApplication() override{}

	///@}
	///@name Operations
	///@{

	void Register() override;

	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	std::string Info() const override
    {
		return "KratosCompressiblePotentialFlowApplication";
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
  		KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );

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

private:
    ///@name Member Variables
    ///@{

    const IncompressiblePotentialFlowElement<2, 3> mIncompressiblePotentialFlowElement2D3N;
    const IncompressiblePotentialFlowElement<3, 4> mIncompressiblePotentialFlowElement3D4N;
    const CompressiblePotentialFlowElement<2, 3> mCompressiblePotentialFlowElement2D3N;
    const CompressiblePotentialFlowElement<3, 4> mCompressiblePotentialFlowElement3D4N;
    const IncompressiblePerturbationPotentialFlowElement<2, 3> mIncompressiblePerturbationPotentialFlowElement2D3N;
    const IncompressiblePerturbationPotentialFlowElement<3, 4> mIncompressiblePerturbationPotentialFlowElement3D4N;
    const CompressiblePerturbationPotentialFlowElement<2, 3> mCompressiblePerturbationPotentialFlowElement2D3N;
    const CompressiblePerturbationPotentialFlowElement<3, 4> mCompressiblePerturbationPotentialFlowElement3D4N;
    const AdjointAnalyticalIncompressiblePotentialFlowElement<IncompressiblePotentialFlowElement<2, 3>> mAdjointAnalyticalIncompressiblePotentialFlowElement2D3N;
    const AdjointFiniteDifferencePotentialFlowElement<IncompressiblePotentialFlowElement<2,3>> mAdjointIncompressiblePotentialFlowElement2D3N;
    const AdjointFiniteDifferencePotentialFlowElement<CompressiblePotentialFlowElement<2,3>> mAdjointCompressiblePotentialFlowElement2D3N;
    const EmbeddedIncompressiblePotentialFlowElement<2,3> mEmbeddedIncompressiblePotentialFlowElement2D3N;
    const EmbeddedIncompressiblePotentialFlowElement<3,4> mEmbeddedIncompressiblePotentialFlowElement3D4N;
    const EmbeddedCompressiblePotentialFlowElement<2,3> mEmbeddedCompressiblePotentialFlowElement2D3N;
    const EmbeddedCompressiblePotentialFlowElement<3,4> mEmbeddedCompressiblePotentialFlowElement3D4N;
    const AdjointFiniteDifferencePotentialFlowElement<EmbeddedIncompressiblePotentialFlowElement<2,3>> mAdjointEmbeddedIncompressiblePotentialFlowElement2D3N;
    const AdjointFiniteDifferencePotentialFlowElement<EmbeddedCompressiblePotentialFlowElement<2,3>> mAdjointEmbeddedCompressiblePotentialFlowElement2D3N;


    const PotentialWallCondition<2,2> mPotentialWallCondition2D2N;
    const PotentialWallCondition<3,3> mPotentialWallCondition3D3N;
    const AdjointPotentialWallCondition<PotentialWallCondition<2,2>> mAdjointPotentialWallCondition2D2N;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosCompressiblePotentialFlowApplication& operator=(KratosCompressiblePotentialFlowApplication const& rOther);

    /// Copy constructor.
    KratosCompressiblePotentialFlowApplication(KratosCompressiblePotentialFlowApplication const& rOther);

    ///@}

}; // Class KratosCompressiblePotentialFlowApplication

///@}

}  // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_H_INCLUDED  defined
