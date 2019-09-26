//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


#if !defined(KRATOS_FLUID_TRANSPORT_APPLICATION_H_INCLUDED )
#define  KRATOS_FLUID_TRANSPORT_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

// Application includes
#include "fluid_transport_application_variables.h"

#include "custom_elements/steady_convection_diffusion_FIC_element.hpp"
#include "custom_elements/transient_convection_diffusion_FIC_element.hpp"
#include "custom_elements/transient_convection_diffusion_FIC_explicit_element.hpp"
#include "custom_elements/transient_convection_diffusion_PFEM2_FIC_element.hpp"


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
class KratosFluidTransportApplication : public KratosApplication
{

public:
	///@name Type Definitions
	///@{


	/// Pointer definition of KratosFluidTransportApplication
	KRATOS_CLASS_POINTER_DEFINITION(KratosFluidTransportApplication);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosFluidTransportApplication();

	/// Destructor.
	~KratosFluidTransportApplication() override {}


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
		return "KratosFluidTransportApplication";
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

	const SteadyConvectionDiffusionFICElement<2,3> mSteadyConvectionDiffusionFICElement2D3N;
	const SteadyConvectionDiffusionFICElement<2,4> mSteadyConvectionDiffusionFICElement2D4N;
	const SteadyConvectionDiffusionFICElement<3,4> mSteadyConvectionDiffusionFICElement3D4N;
	const SteadyConvectionDiffusionFICElement<3,8> mSteadyConvectionDiffusionFICElement3D8N;

	const TransientConvectionDiffusionFICElement<2,3> mTransientConvectionDiffusionFICElement2D3N;
	const TransientConvectionDiffusionFICElement<2,4> mTransientConvectionDiffusionFICElement2D4N;
	const TransientConvectionDiffusionFICElement<3,4> mTransientConvectionDiffusionFICElement3D4N;
	const TransientConvectionDiffusionFICElement<3,8> mTransientConvectionDiffusionFICElement3D8N;

	const TransientConvectionDiffusionFICExplicitElement<2,3> mTransientConvectionDiffusionFICExplicitElement2D3N;
	const TransientConvectionDiffusionFICExplicitElement<2,4> mTransientConvectionDiffusionFICExplicitElement2D4N;
	const TransientConvectionDiffusionFICExplicitElement<3,4> mTransientConvectionDiffusionFICExplicitElement3D4N;
	const TransientConvectionDiffusionFICExplicitElement<3,8> mTransientConvectionDiffusionFICExplicitElement3D8N;

	const TransientConvectionDiffusionPFEM2FICElement<2,3> mTransientConvectionDiffusionPFEM2FICElement2D3N;
	const TransientConvectionDiffusionPFEM2FICElement<2,4> mTransientConvectionDiffusionPFEM2FICElement2D4N;
	const TransientConvectionDiffusionPFEM2FICElement<3,4> mTransientConvectionDiffusionPFEM2FICElement3D4N;
	const TransientConvectionDiffusionPFEM2FICElement<3,8> mTransientConvectionDiffusionPFEM2FICElement3D8N;

	// const Elem2D   mElem2D;
	// const Elem3D   mElem3D;

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
	KratosFluidTransportApplication& operator=(KratosFluidTransportApplication const& rOther);

	/// Copy constructor.
	KratosFluidTransportApplication(KratosFluidTransportApplication const& rOther);


	///@}

}; // Class KratosFluidTransportApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_MY_LAPLACIAN_APPLICATION_H_INCLUDED  defined
