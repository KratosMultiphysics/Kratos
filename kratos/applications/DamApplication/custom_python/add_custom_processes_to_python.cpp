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

#include "custom_processes/exact_evolution_conditions_load_process.hpp"
#include "custom_processes/exact_evolution_conditions_temperature_process.hpp"
#include "custom_processes/interpolation_evolution_conditions_load_process.hpp"
#include "custom_processes/interpolation_evolution_conditions_temperature_process.hpp"
#include "custom_processes/linear_evolution_conditions_temperature_process.hpp"


namespace Kratos
{
	
namespace Python
{

using namespace boost::python;

void  AddCustomProcessesToPython() 
{
    typedef Process                                                 ProcessBaseType;
    typedef ExactEvolutionConditionsLoadProcess                     ExactEvolutionConditionsLoadProcessType;
    typedef ExactEvolutionConditionsTemperatureProcess              ExactEvolutionConditionsTemperatureProcessType;
    typedef InterpolationEvolutionConditionsLoadProcess             InterpolationEvolutionConditionsLoadProcessType;
    typedef InterpolationEvolutionConditionsTemperatureProcess      InterpolationEvolutionConditionsTemperatureProcessType;
    typedef LinearEvolutionConditionsTemperatureProcess             LinearEvolutionConditionsTemperatureProcessType;
   
    // Exact case Evolution for Normal Load
    class_< ExactEvolutionConditionsLoadProcessType, bases< ProcessBaseType >, boost::noncopyable > ( "ExactEvolutionConditionsLoadProcess",
        init < ModelPart&,double >());
        
    // Exact case Evolution for Temperature
    class_< ExactEvolutionConditionsTemperatureProcessType, bases< ProcessBaseType >, boost::noncopyable > ( "ExactEvolutionConditionsTemperatureProcess",
        init < ModelPart&,double >());
        
    // Interpolation case Evolution for Normal Load
    class_< InterpolationEvolutionConditionsLoadProcessType, bases< ProcessBaseType >, boost::noncopyable > ( "InterpolationEvolutionConditionsLoadProcess",
        init < ModelPart&,double >());
        
    // Interpolation case Evolution for Temperature
    class_< InterpolationEvolutionConditionsTemperatureProcessType, bases< ProcessBaseType >, boost::noncopyable > ( "InterpolationEvolutionConditionsTemperatureProcess",
        init < ModelPart&,double >());        

    // Linear case Evolution for Temperature
    class_< LinearEvolutionConditionsTemperatureProcessType, bases< ProcessBaseType >, boost::noncopyable > ( "LinearEvolutionConditionsTemperatureProcess",
        init < ModelPart& >());                
}

}  // namespace Python.
} // Namespace Kratos

