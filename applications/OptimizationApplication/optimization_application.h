//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

#if !defined(KRATOS_OPTIMIZATION_APPLICATION_H_INCLUDED )
#define  KRATOS_OPTIMIZATION_APPLICATION_H_INCLUDED

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

#include <string>
#include <iostream>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/kratos_application.h"

/* ELEMENTS */
#include "custom_elements/helmholtz_surf_shape_element.h"
#include "custom_elements/helmholtz_surf_thickness_element.h"
#include "custom_elements/helmholtz_bulk_shape_element.h"
#include "custom_elements/helmholtz_bulk_element.h"

/* ADJOINT ELEMENTS */
#include "custom_elements/adjoint_small_displacement_element.h"

/* CONDITIONS */
#include "custom_conditions/helmholtz_surf_shape_condition.h"

// ==============================================================================

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
	class KRATOS_API(OPTIMIZATION_APPLICATION) KratosOptimizationApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{


		/// Pointer definition of KratosOptimizationApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosOptimizationApplication);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		KratosOptimizationApplication();

		/// Destructor.
		~KratosOptimizationApplication() override {}


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
			return "KratosOptimizationApplication";
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



		//       static const ApplicationCondition  msApplicationCondition;

		///@}
		///@name Member Variables
		///@{
		/* ELEMENTS */

		const HelmholtzSurfShapeElement mHelmholtzSurfShape3D3N;
		const HelmholtzSurfThicknessElement mHelmholtzSurfThickness3D3N;
		const HelmholtzBulkShapeElement mHelmholtzBulkShape3D4N;
		const HelmholtzBulkElement mHelmholtzBulkTopology3D4N;

		/* ADJ ELEMENTS */
		const AdjointSmallDisplacementElement mAdjointSmallDisplacementElement3D4N;

		/* CONDITIONS*/
		// Surface conditions
		const HelmholtzSurfShapeCondition mHelmholtzSurfShapeCondition3D3N;

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
		KratosOptimizationApplication& operator=(KratosOptimizationApplication const& rOther);

		/// Copy constructor.
		KratosOptimizationApplication(KratosOptimizationApplication const& rOther);


		///@}

	}; // Class KratosOptimizationApplication

	///@}


	///@name Type Definitions
	///@{


	///@}
	///@name Input and output
	///@{

	///@}


}  // namespace Kratos.

#endif // KRATOS_OPTIMIZATION_APPLICATION_H_INCLUDED  defined


