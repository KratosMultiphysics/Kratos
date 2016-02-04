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

#include "custom_processes/evolution_conditions_process.hpp"
#include "custom_processes/interpolated_evolution_conditions_process.hpp"

namespace Kratos
{
	
namespace Python
{

using namespace boost::python;

void  AddCustomProcessesToPython() 
{
    typedef Process                                       ProcessBaseType;
    typedef EvolutionConditionsProcess                    EvolutionConditionsProcessType;
    typedef InterpolatedEvolutionConditionsProcess        InterpolatedEvolutionConditionsProcessType;

    
    // Evolution Conditions Process
    class_< EvolutionConditionsProcessType, bases< ProcessBaseType >, boost::noncopyable > ( "EvolutionConditionsProcess",
        init < ModelPart&,double >());
    
    // Interpolated Evolution Conditions Process
    class_< InterpolatedEvolutionConditionsProcessType, bases< ProcessBaseType >, boost::noncopyable > ( "InterpolatedEvolutionConditionsProcess",
        init < ModelPart&,double >());
}
	
}  // namespace Python.
} // Namespace Kratos

