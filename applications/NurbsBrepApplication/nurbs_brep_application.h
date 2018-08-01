//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


#if !defined(KRATOS_NURBS_BREP_APPLICATION_H_INCLUDED )
#define  KRATOS_NURBS_BREP_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"

#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

#include "includes/model_part.h"


namespace Kratos {

///@name Kratos Globals
///@{

	//Variables definition
	//KRATOS_DEFINE_VARIABLE(Tree::Pointer, SEARCH_TREE)
	//KRATOS_DEFINE_VARIABLE(NurbsBrepModeler::tree::Pointer, SEARCH_TREE)
	//variables for IGA - DEM coupling
	KRATOS_DEFINE_VARIABLE(std::vector<Condition*>, WALL_POINT_CONDITION_POINTERS)
	typedef std::vector<array_1d<double, 3> > std_vector_of_arrays_3d;
  	KRATOS_DEFINE_APPLICATION_VARIABLE(DEM_APPLICATION, std_vector_of_arrays_3d, WALL_POINT_CONDITION_ELASTIC_FORCES)
  	KRATOS_DEFINE_APPLICATION_VARIABLE(DEM_APPLICATION, std_vector_of_arrays_3d, WALL_POINT_CONDITION_TOTAL_FORCES)

	KRATOS_DEFINE_VARIABLE(double, RADIUS)
	KRATOS_DEFINE_VARIABLE(Vector, COORDINATES)

	KRATOS_DEFINE_VARIABLE(double, CONTROL_POINT_WEIGHT)

	//KRATOS_DEFINE_VARIABLE(double, INTEGRATION_WEIGHT)

	KRATOS_DEFINE_VARIABLE(Vector, TANGENTS_BASIS_VECTOR)

	KRATOS_DEFINE_VARIABLE(Vector, LOCAL_PARAMETERS)
	KRATOS_DEFINE_VARIABLE(int,    FACE_BREP_ID)

	KRATOS_DEFINE_VARIABLE(Vector, CONTROL_POINT_IDS)
	KRATOS_DEFINE_VARIABLE(Vector, CONTROL_POINT_IDS_SLAVE)

    KRATOS_DEFINE_VARIABLE(Vector, SHAPE_FUNCTION_VALUES)
    KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES)
    KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES)
	KRATOS_DEFINE_VARIABLE(Vector, NURBS_SHAPE_FUNCTIONS)
	KRATOS_DEFINE_VARIABLE(Matrix, NURBS_SHAPE_FUNCTION_DERIVATIVES)
	KRATOS_DEFINE_VARIABLE(Matrix, NURBS_SHAPE_FUNCTION_SECOND_DERIVATIVES)

	KRATOS_DEFINE_VARIABLE(Vector, TANGENTS_BASIS_VECTOR_SLAVE)

	KRATOS_DEFINE_VARIABLE(Vector, LOCATION_SLAVE)
	KRATOS_DEFINE_VARIABLE(Vector, LOCAL_PARAMETERS_SLAVE)
	KRATOS_DEFINE_VARIABLE(int,    FACE_BREP_ID_SLAVE)

	KRATOS_DEFINE_VARIABLE(Vector, NURBS_SHAPE_FUNCTIONS_SLAVE)
	KRATOS_DEFINE_VARIABLE(Matrix, NURBS_SHAPE_FUNCTION_DERIVATIVES_SLAVE)
	KRATOS_DEFINE_VARIABLE(Matrix, NURBS_SHAPE_FUNCTION_SECOND_DERIVATIVES_SLAVE)

	KRATOS_DEFINE_VARIABLE(double, CURVE_PARAMETER_KNOT_LOCATION_PERCENTAGE)

///@}
class KratosNurbsBrepApplication : public KratosApplication {
public:
	///@name Type Definitions
	///@{
	/// Pointer definition of KratosNurbsBrepApplication
	KRATOS_CLASS_POINTER_DEFINITION(KratosNurbsBrepApplication);
	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosNurbsBrepApplication();

	/// Destructor.
	virtual ~KratosNurbsBrepApplication() {}

	///@}
	///@name Operations
	///@{

	virtual void Register();

	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	virtual std::string Info() const {
		return "KratosNurbsBrepApplication";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const {
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
protected:

private:
	///@name Un accessible methods
	///@{

	/// Assignment operator.
	KratosNurbsBrepApplication& operator=(KratosNurbsBrepApplication const& rOther);

	/// Copy constructor.
	KratosNurbsBrepApplication(KratosNurbsBrepApplication const& rOther);
	///@}

	}; // Class KratosNurbsBrepApplication
}  // namespace Kratos.

#endif // KRATOS_NURBS_BREP_APPLICATION_H_INCLUDED  defined
