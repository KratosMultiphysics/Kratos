//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

#include "custom_processes/apply_component_table_process.hpp"
#include "custom_processes/apply_double_table_process.hpp"
#include "custom_processes/apply_constant_hydrostatic_pressure_process.hpp"
#include "custom_processes/apply_hydrostatic_pressure_table_process.hpp"


namespace Kratos
{
	
namespace Python
{

using namespace boost::python;

void  AddCustomProcessesToPython() 
{
    class_< ApplyComponentTableProcess, bases< Process >, boost::noncopyable > ( "ApplyComponentTableProcess",
        init < ModelPart&, Parameters>());
        
    class_< ApplyDoubleTableProcess, bases< Process >, boost::noncopyable > ( "ApplyDoubleTableProcess",
        init < ModelPart&, Parameters>());

    class_< ApplyConstantHydrostaticPressureProcess, bases< Process >, boost::noncopyable > ( "ApplyConstantHydrostaticPressureProcess",
        init < ModelPart&, Parameters>());

    class_< ApplyHydrostaticPressureTableProcess, bases< Process >, boost::noncopyable > ( "ApplyHydrostaticPressureTableProcess",
        init < ModelPart&, Parameters>());
}

}  // namespace Python.
} // Namespace Kratos
