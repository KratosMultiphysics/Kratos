//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:              March 2017 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_ADD_CUSTOM_MPI_STRATEGIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_MPI_STRATEGIES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos
{

namespace Python
{

void  AddCustomMPIStrategiesToPython(pybind11::module& m);

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_MPI_STRATEGIES_TO_PYTHON_H_INCLUDED  defined 
