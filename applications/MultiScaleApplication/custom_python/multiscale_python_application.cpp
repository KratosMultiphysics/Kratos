/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "multiscale_application.h"
#include "add_linear_solvers_to_python.h"
#include "add_conditions_to_python.h"
#include "add_strategies_to_python.h"
#include "add_utilities_to_python.h"
#include "add_constitutive_laws_to_python.h"
#include "add_custom_io_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

BOOST_PYTHON_MODULE( KratosMultiScaleApplication )
{

    class_ < KratosMultiScaleApplication,
           KratosMultiScaleApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable > ( "KratosMultiScaleApplication" )
           ;

	AddConditionsToPython();
	AddLinearSolversToPython();
	AddStrategiesToPython();
	AddUtilitiesToPython();
	AddConstitutiveLawsToPython();
	AddCustomIOToPython();

    //registering variables in python

	// for lagrange multipliers
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_LAGRANGE )
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ROTATION_LAGRANGE )
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_DOUBLE_LAGRANGE_1 )
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_DOUBLE_LAGRANGE_2 )
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ROTATION_DOUBLE_LAGRANGE_1 )
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ROTATION_DOUBLE_LAGRANGE_2 )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( DOUBLE_LAGRANGE_SCALE_FACTOR )

	// for strategies
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( STRATEGY_SOLUTION_STEP_SOLVED )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONSTITUTIVE_INTAGRATION_ERROR_CODE )

	// for damage constitutive law
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( DAMAGE_T )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( DAMAGE_C )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRACTURE_ENERGY_T )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRACTURE_ENERGY_C )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_STRESS_T )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_STRESS_C )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( DAMAGE_SECANT_MATRIX )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( DAMAGE_MODEL )

	// for custom fracture-energy-based regularization
	KRATOS_REGISTER_VARIABLE( CHARACTERISTIC_LENGTH_MULTIPLIER )

	// for interface constitutive law
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( NORMAL_STIFFNESS )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( TANGENTIAL_STIFFNESS )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( INITIAL_COHESION )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRACTURE_ENERGY_MODE_I )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRACTURE_ENERGY_MODE_II )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_FUNCTION_VALUE )

	// for plots
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_SURFACE_DATA_2D_X )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_SURFACE_DATA_2D_Y )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_SURFACE_DATA_3D_X )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_SURFACE_DATA_3D_Y )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_SURFACE_DATA_3D_Z )

	// for plastic constitutive law
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( ISOTROPIC_HARDENING )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( KINEMATIC_HARDENING )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_STRESS_INFINITY )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( ISOTROPIC_HARDENING_EXPONENT )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( EQUIVALENT_PLASTIC_STRAIN )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( PLASTIC_STRAIN_TENSOR )
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
