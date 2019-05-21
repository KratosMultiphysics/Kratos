//
//  Main authors:    Miguel Angel Celigueta   maceli@cimne.upc.edu
//
//


#if !defined(KRATOS_STRATEGIES_PYTHON_H_INCLUDED )
#define  KRATOS_STRATEGIES_PYTHON_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"


namespace Kratos
{

    namespace Python
    {

      void  AddCustomStrategiesToPython(pybind11::module& m);

    }  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_STRATEGIES_PYTHON_H_INCLUDED  defined
