// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_TOPOLOGYOPTIMIZATION_APPLICATION_H_INCLUDED )
#define  KRATOS_TOPOLOGYOPTIMIZATION_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>

// External includes

// Core applications
#include "topology_optimization_application.h"

#include "custom_elements/small_displacement_simp_element.hpp"
#include "includes/variables.h"

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

//elements


namespace Kratos
{

	///@name Kratos Globals
	///@{

    // Variables definition with Python connection
    KRATOS_DEFINE_VARIABLE( double, E_MIN )
    KRATOS_DEFINE_VARIABLE( double, E_0 )
    KRATOS_DEFINE_VARIABLE( double, PENAL )
    KRATOS_DEFINE_VARIABLE( double, X_PHYS )
    KRATOS_DEFINE_VARIABLE( double, X_PHYS_OLD )
    KRATOS_DEFINE_VARIABLE( double, DCDX )
    KRATOS_DEFINE_VARIABLE( double, DVDX )
    KRATOS_DEFINE_VARIABLE( double, SOLID_VOID )
    KRATOS_DEFINE_VARIABLE( double, LOCAL_STRAIN_ENERGY )


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
	class KratosTopologyOptimizationApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{


		/// Pointer definition of KratosTopologyOptimizationApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosTopologyOptimizationApplication);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		KratosTopologyOptimizationApplication();

		/// Destructor.
		~KratosTopologyOptimizationApplication() override {}


		///@}
		///@name Operators
		///@{


		///@}
		///@name Operations
		///@{

		virtual void Register();



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
		virtual std::string Info() const
		{
			return "KratosTopologyOptimizationApplication";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}

		///// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
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

		///@}
		///@name Member Variables
		///@{

        //small_displacement
      	const SmallDisplacementSIMPElement mSmallDisplacementSIMPElement3D3N; // dummy element for surface representation
        const SmallDisplacementSIMPElement mSmallDisplacementSIMPElement3D4N;
        const SmallDisplacementSIMPElement mSmallDisplacementSIMPElement3D8N;

//        Extra elements to be added in the future
//        const SmallDisplacementSIMPElement mSmallDisplacementSIMPElement3D6N;
//        const SmallDisplacementSIMPElement mSmallDisplacementSIMPElement3D10N;
//        const SmallDisplacementSIMPElement mSmallDisplacementSIMPElement3D15N;
//        const SmallDisplacementSIMPElement mSmallDisplacementSIMPElement3D20N;
//        const SmallDisplacementSIMPElement mSmallDisplacementSIMPElement3D27N;


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
		KratosTopologyOptimizationApplication& operator=(KratosTopologyOptimizationApplication const& rOther);

		/// Copy constructor.
		KratosTopologyOptimizationApplication(KratosTopologyOptimizationApplication const& rOther);


		///@}

	}; // Class KratosTopologyOptimizationApplication

	///@}


	///@name Type Definitions
	///@{


	///@}
	///@name Input and output
	///@{

	///@}


}  // namespace Kratos.

#endif // KRATOS_TOPOLOGYOPTIMIZATION_APPLICATION_H_INCLUDED  defined


