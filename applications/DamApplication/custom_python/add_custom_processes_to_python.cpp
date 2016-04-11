//   
//   Project Name:           
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//

// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes 
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_python/add_custom_processes_to_python.h"

// Processes
#include "custom_processes/displa_table_interpolation_process.hpp"
#include "custom_processes/temperature_table_interpolation_process.hpp"
#include "custom_processes/point_load_table_interpolation_process.hpp"
#include "custom_processes/line_load_table_interpolation_process.hpp"
#include "custom_processes/surface_load_table_interpolation_process.hpp"
#include "custom_processes/normal_load_table_interpolation_process.hpp"
#include "custom_processes/exact_water_evolution_conditions_load_process.hpp"
#include "custom_processes/exact_bofang_evolution_conditions_temperature_process.hpp"
#include "custom_processes/interpolation_water_evolution_conditions_load_process.hpp"
#include "custom_processes/interpolation_bofang_evolution_conditions_temperature_process.hpp"


namespace Kratos
{
	
namespace Python
{

using namespace boost::python;

void  AddCustomProcessesToPython() 
{
    // Interpolation table for displacements
    class_< DisplaTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "DisplaTableInterpolationProcess",
        init < ModelPart&,double >());

    // Interpolation table for uniform temperature
    class_< TemperatureTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "TemperatureTableInterpolationProcess",
        init < ModelPart&,double >());

    // Interpolation table for point loads
    class_< PointLoadTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "PointLoadTableInterpolationProcess",
        init < ModelPart&,double >());

    // Interpolation table for line loads
    class_< LineLoadTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "LineLoadTableInterpolationProcess",
        init < ModelPart&,double >());

    // Interpolation table for surface loads
    class_< SurfaceLoadTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "SurfaceLoadTableInterpolationProcess",
        init < ModelPart&,double >());  
    
    // Interpolation table for Uniform Normal loads    
    class_< NormalLoadTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "NormalLoadTableInterpolationProcess",
        init < ModelPart&,double >());
   
    // Exact case Evolution for Temperature
    class_< ExactBofangEvolutionConditionsTemperatureProcess, bases< Process >, boost::noncopyable > ( "ExactBofangEvolutionConditionsTemperatureProcess",
        init < ModelPart&,double >());
    
    // Interpolation case Evolution for Temperature
    class_< InterpolationBofangEvolutionConditionsTemperatureProcess, bases< Process >, boost::noncopyable > ( "InterpolationBofangEvolutionConditionsTemperatureProcess",
        init < ModelPart&,double >());
        
    // Exact case Evolution for Water Loads
    class_< ExactWaterEvolutionConditionsLoadProcess, bases< Process >, boost::noncopyable > ( "ExactWaterEvolutionConditionsLoadProcess",
        init < ModelPart&,double >());
        
    // Interpolation case Evolution for Water Loads
    class_< InterpolationWaterEvolutionConditionsLoadProcess, bases< Process >, boost::noncopyable > ( "InterpolationWaterEvolutionConditionsLoadProcess",
        init < ModelPart&,double >());
              
}

}  // namespace Python.
} // Namespace Kratos

