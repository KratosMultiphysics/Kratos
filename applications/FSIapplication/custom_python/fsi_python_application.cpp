//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-03-25 07:48:22 $
//   Revision:            $Revision: 1.5 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "fsi_application.h"
#include "custom_python/add_convergence_accelerators_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_mappers_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosFSIApplication)
{

    class_<KratosFSIApplication,
           KratosFSIApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosFSIApplication")
           ;

    AddCustomUtilitiesToPython();
    AddMappersToPython();
    AddConvergenceAcceleratorsToPython();

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
