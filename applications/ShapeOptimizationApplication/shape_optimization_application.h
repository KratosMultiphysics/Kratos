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
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DF1DX);

	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DF1DX_MAPPED);

	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC1DX);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC2DX);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC3DX);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC4DX);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC5DX);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC6DX);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC7DX);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC8DX);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC9DX);

	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC1DX_MAPPED);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC2DX_MAPPED);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC3DX_MAPPED);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC4DX_MAPPED);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC5DX_MAPPED);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC6DX_MAPPED);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC7DX_MAPPED);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC8DX_MAPPED);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DC9DX_MAPPED);

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

	// For bead optimization
    KRATOS_DEFINE_VARIABLE(double,ALPHA);
    KRATOS_DEFINE_VARIABLE(double,ALPHA_MAPPED);
    KRATOS_DEFINE_VARIABLE(double,DF1DALPHA);
    KRATOS_DEFINE_VARIABLE(double,DF1DALPHA_MAPPED);
    KRATOS_DEFINE_VARIABLE(double,DPDALPHA);
    KRATOS_DEFINE_VARIABLE(double,DPDALPHA_MAPPED);
    KRATOS_DEFINE_VARIABLE(double,DLDALPHA);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BEAD_DIRECTION);

    // For auxiliary operations
    KRATOS_DEFINE_VARIABLE(double,SCALAR_VARIABLE);
    KRATOS_DEFINE_VARIABLE(double,SCALAR_VARIABLE_MAPPED);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(VECTOR_VARIABLE);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(VECTOR_VARIABLE_MAPPED);

	// For in plane mapping operations
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BACKGROUND_COORDINATE);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BACKGROUND_NORMAL);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(OUT_OF_PLANE_DELTA);


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
		~KratosShapeOptimizationApplication() override {}


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
			return "KratosShapeOptimizationApplication";
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


