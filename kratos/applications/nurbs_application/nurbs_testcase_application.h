//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_NURBSTESTCASE_APPLICATION_H_INCLUDED )
#define  KRATOS_NURBSTESTCASE_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "geometries/nurbs_3d.h"
#include "geometries/nurbs_2d.h"
#include "includes/variables.h"
#include "custom_elements/nurbs_poisson_2d.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 


    KRATOS_DEFINE_VARIABLE(int, NURBS_ID);
//    KRATOS_DEFINE_VARIABLE(double, DISTANCE);
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(NURBS_COORDINATES);




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
	class KratosNurbsTestcaseApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosNurbsTestcaseApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosNurbsTestcaseApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosNurbsTestcaseApplication();

		/// Destructor.
		virtual ~KratosNurbsTestcaseApplication(){}


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
			return "KratosNurbsTestcaseApplication";
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



		//       static const ApplicationCondition  msApplicationCondition; 

		///@} 
        ///@name Member Variables
		///@{ 

                NurbsPatchGeometry3D< Node<3> > mNurbsPatchGeometry3D;
                const NurbsPoisson2D   mNurbsPoisson2D;


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
		KratosNurbsTestcaseApplication& operator=(KratosNurbsTestcaseApplication const& rOther);

		/// Copy constructor.
		KratosNurbsTestcaseApplication(KratosNurbsTestcaseApplication const& rOther);


		///@}    

	}; // Class KratosNurbsTestcaseApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_NURBSTESTCASE_APPLICATION_H_INCLUDED  defined 


