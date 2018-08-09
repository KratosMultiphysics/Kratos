//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/NurbsBrepApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes


// Project includes
#include "includes/define_python.h"
#include "includes/define.h"
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

	namespace Python
	{

		using namespace pybind11;


		PYBIND11_MODULE(KratosNurbsBrepApplication, m)
		{

			class_<KratosNurbsBrepApplication,
				KratosNurbsBrepApplication::Pointer,
				KratosApplication >(m, "KratosNurbsBrepApplication")
				.def(init<>())
				;

			AddCustomStrategiesToPython(m);
			AddCustomUtilitiesToPython(m);

			//To enhance weighting to nodes:
			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COORDINATES)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SURFACE_NORMAL)

			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTROL_POINT_WEIGHT)
			//KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTEGRATION_WEIGHT)

			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TANGENTS_BASIS_VECTOR)

			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_PARAMETERS)
			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FACE_BREP_ID)

			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTROL_POINT_IDS)
			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTROL_POINT_IDS_SLAVE)

			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NURBS_SHAPE_FUNCTIONS)
			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NURBS_SHAPE_FUNCTION_DERIVATIVES)
			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NURBS_SHAPE_FUNCTION_SECOND_DERIVATIVES)

			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TANGENTS_BASIS_VECTOR_SLAVE)

			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCATION_SLAVE)
			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_PARAMETERS_SLAVE)
			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FACE_BREP_ID_SLAVE)

			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NURBS_SHAPE_FUNCTIONS_SLAVE)
			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NURBS_SHAPE_FUNCTION_DERIVATIVES_SLAVE)
			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NURBS_SHAPE_FUNCTION_SECOND_DERIVATIVES_SLAVE)

			// internal variable to optimize projections between master and slave curves.
			KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CURVE_PARAMETER_KNOT_LOCATION_PERCENTAGE)
		}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
