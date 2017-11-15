// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#if !defined(KRATOS_SHAPEOPTIMIZATION_APPLICATION_H_INCLUDED )
#define  KRATOS_SHAPEOPTIMIZATION_APPLICATION_H_INCLUDED

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
#include "includes/define.h"
#include "includes/kratos_application.h"

// elements
#include "custom_elements/small_displacement_analytic_sensitivity_element.hpp"

//conditions
#include "custom_conditions/shape_optimization_condition.h"

// Variables
#include "includes/variables.h"

// ==============================================================================

namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Geometry variables
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(NORMALIZED_SURFACE_NORMAL);

    // Optimization variables
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(OBJECTIVE_SENSITIVITY);
    KRATOS_DEFINE_VARIABLE(double,OBJECTIVE_SURFACE_SENSITIVITY);
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MAPPED_OBJECTIVE_SENSITIVITY);
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONSTRAINT_SENSITIVITY);
    KRATOS_DEFINE_VARIABLE(double,CONSTRAINT_SURFACE_SENSITIVITY);
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MAPPED_CONSTRAINT_SENSITIVITY);
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(SEARCH_DIRECTION);
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_UPDATE);
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_CHANGE);
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(SHAPE_UPDATE);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(SHAPE_CHANGE);
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_CHANGE);	

	// For edge damping
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DAMPING_FACTOR);

    // For mapping
    KRATOS_DEFINE_VARIABLE(int,MAPPING_ID);

    // For Structure Sensitivity Analysis
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(STRAIN_ENERGY_SHAPE_GRADIENT);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MASS_SHAPE_GRADIENT);
    KRATOS_DEFINE_VARIABLE(int,ACTIVE_NODE_INDEX);
	KRATOS_DEFINE_VARIABLE(Vector,DKDXU);
    KRATOS_DEFINE_VARIABLE(Vector,DKDXU_X);
    KRATOS_DEFINE_VARIABLE(Vector,DKDXU_Y);
    KRATOS_DEFINE_VARIABLE(Vector,DKDXU_Z);

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
	class KratosShapeOptimizationApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosShapeOptimizationApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosShapeOptimizationApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosShapeOptimizationApplication();

		/// Destructor.
		virtual ~KratosShapeOptimizationApplication(){}


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
			return "KratosShapeOptimizationApplication";
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

        // elements

        // for structural optimization
      	const SmallDisplacementAnalyticSensitivityElement mSmallDisplacementAnalyticSensitivityElement3D4N;
      	const SmallDisplacementAnalyticSensitivityElement mSmallDisplacementAnalyticSensitivityElement3D10N;
      	const SmallDisplacementAnalyticSensitivityElement mSmallDisplacementAnalyticSensitivityElement3D8N;
      	const SmallDisplacementAnalyticSensitivityElement mSmallDisplacementAnalyticSensitivityElement3D20N;

        //conditions
        const ShapeOptimizationCondition mShapeOptimizationCondition3D3N;
		const ShapeOptimizationCondition mShapeOptimizationCondition3D4N;
        const ShapeOptimizationCondition mShapeOptimizationCondition2D2N;


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
		KratosShapeOptimizationApplication& operator=(KratosShapeOptimizationApplication const& rOther);

		/// Copy constructor.
		KratosShapeOptimizationApplication(KratosShapeOptimizationApplication const& rOther);


		///@}    

	}; // Class KratosShapeOptimizationApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_SHAPEOPTIMIZATION_APPLICATION_H_INCLUDED  defined 


