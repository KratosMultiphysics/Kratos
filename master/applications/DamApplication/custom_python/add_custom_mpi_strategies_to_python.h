//   
//   Project Name:        
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
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
