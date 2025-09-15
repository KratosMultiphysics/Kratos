//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED  defined 
