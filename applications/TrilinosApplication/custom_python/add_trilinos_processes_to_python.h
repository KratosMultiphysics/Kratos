

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2011-07-22 15:20:00 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_TRILINOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_TRILINOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED



// System includes


// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

namespace Kratos
{
namespace Python
{
void  AddProcesses(pybind11::module& m);
}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_TRILINOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED  defined
