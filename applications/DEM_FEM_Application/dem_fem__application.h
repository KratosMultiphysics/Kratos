//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_DEM_FEM__APPLICATION_H_INCLUDED )
#define  KRATOS_DEM_FEM__APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


#include "includes/variables.h"



//////////////Cfeng:120416 from structual application below


#include "includes/serializer.h"








namespace Kratos
{

      KRATOS_DEFINE_VARIABLE(Vector, PRESTRESS)
      KRATOS_DEFINE_VARIABLE( Matrix, GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
      KRATOS_DEFINE_VARIABLE(double,  PARTICLE_NUMBER_OF_NEIGHBOURS)
      
  
      KRATOS_DEFINE_VARIABLE(double,  DEM_FEM_CONVERGENCE_RATIO)

      KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FAILURE_ID)
      KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FORCE)
      KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_IF_INITIAL_CONTACT)
      KRATOS_DEFINE_VARIABLE(WeakPointerVector<Condition >,     NEIGHBOUR_RIGID_FACES)

      KRATOS_DEFINE_VARIABLE(int, IF_BOUNDARY_ELEMENT)
	  KRATOS_DEFINE_VARIABLE(int, ROTATION_OPTION)
      KRATOS_DEFINE_VARIABLE(Vector, IF_BOUNDARY_FACE)


      KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_MOMENT );
      KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_ROTATION_ANGLE );
      KRATOS_DEFINE_VARIABLE(double,  PARTICLE_MOMENT_OF_INERTIA);

	
	class KratosDEM_FEM_Application : public KratosApplication
	{
	public:


		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosDEM_FEM_Application
		KRATOS_CLASS_POINTER_DEFINITION(KratosDEM_FEM_Application);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosDEM_FEM_Application();

		/// Destructor.
		virtual ~KratosDEM_FEM_Application(){}


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
			return "KratosDEM_FEM_Application";
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
// 		const Elem2D   mElem2D; 
// 		const Elem3D   mElem3D; 


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
		KratosDEM_FEM_Application& operator=(KratosDEM_FEM_Application const& rOther);

		/// Copy constructor.
		KratosDEM_FEM_Application(KratosDEM_FEM_Application const& rOther);


		///@}    

	}; // Class KratosDEM_FEM_Application

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_DEM_FEM__APPLICATION_H_INCLUDED  defined 


