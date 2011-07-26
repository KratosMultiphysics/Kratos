//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-12-09 20:20:55 $
//   Revision:            $Revision: 1.5 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>

#include "custom_python/add_trilinos_space_to_python.h" 
#include "custom_python/add_trilinos_convergence_criterias_to_python.h" 
#include "custom_python/add_trilinos_schemes_to_python.h" 
#include "custom_python/add_trilinos_linear_solvers_to_python.h" 
#include "custom_python/add_trilinos_strategies_to_python.h" 
#include "custom_python/add_custom_io_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_trilinos_communicator_to_python.h"
#include "custom_python/add_first.h"

////utilities
#include "python/pointer_vector_set_python_interface.h"

// Project includes
#include "includes/define.h"
#include "trilinos_application.h"
//#include "trilinos_space.h"
//#include "spaces/ublas_space.h"
//// #include "add_trilinos_linear_solvers_to_python.h"
//#include "includes/model_part.h"


namespace Kratos {

    namespace Python {


        BOOST_PYTHON_MODULE(KratosTrilinosApplication) {

 	class_<KratosTrilinosApplication,
                    KratosTrilinosApplication::Pointer,
                    bases<KratosApplication>, boost::noncopyable > ("KratosTrilinosApplication")
                    ;

		AddBasicOperations();
		AddConvergenceCriterias();
		AddSchemes();
		AddLinearSolvers();
		AddStrategies();
        AddCustomIOToPython();
        AddCustomUtilitiesToPython();
        AddTrilinosCommunicatorToPython();
	AddFirst();
        }


    } // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
