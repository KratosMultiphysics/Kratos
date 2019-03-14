//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
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
#include "custom_conditions/potential_wall_condition.h"

namespace Kratos {

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KratosCompressiblePotentialFlowApplication : public KratosApplication {
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

    const CompressiblePotentialFlowElement<2,3> mCompressiblePotentialFlowElement2D3N;
    const CompressiblePotentialFlowElement<3,4> mCompressiblePotentialFlowElement3D4N;

    const IncompressiblePotentialFlowElement<2,3> mIncompressiblePotentialFlowElement2D3N;

    const PotentialWallCondition<2,2> mPotentialWallCondition2D2N;
    const PotentialWallCondition<3,3> mPotentialWallCondition3D3N;

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
