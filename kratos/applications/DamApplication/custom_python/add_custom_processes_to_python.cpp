//   
//   Project Name:           
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//

// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

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

// Processes new interface
#include "custom_processes/bofang_condition_temperature_process.hpp"
#include "custom_processes/dam_hydro_condition_load_process.hpp"
#include "custom_processes/dam_uplift_condition_load_process.hpp"
#include "custom_processes/dam_uplift_circular_condition_load_process.hpp"

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
        
        
    // PROCESSES FOR NEW INTERFACE
    
    // Bofang Process
    class_< BofangConditionTemperatureProcess, bases< Process >, boost::noncopyable > ( "BofangConditionTemperatureProcess",
        init < ModelPart&, Parameters>());
        
    // Hydrostatic condition
    class_< DamHydroConditionLoadProcess, bases< Process >, boost::noncopyable > ( "DamHydroConditionLoadProcess",
        init < ModelPart&, Parameters>());
        
    // Uplift Condition
    class_< DamUpliftConditionLoadProcess, bases< Process >, boost::noncopyable > ( "DamUpliftConditionLoadProcess",
        init < ModelPart&, Parameters>());
    
    // Uplift Condition for arch dams   
    class_< DamUpliftCircularConditionLoadProcess, bases< Process >, boost::noncopyable > ( "DamUpliftCircularConditionLoadProcess",
        init < ModelPart&, Parameters>());

}

}  // namespace Python.
} // Namespace Kratos

