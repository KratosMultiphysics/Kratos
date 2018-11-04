//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rzorrilla $
//   Date:                $Date: 2016-9-5 19:30:00 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_ADD_CONVERGENCE_ACCELERATORS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CONVERGENCE_ACCELERATORS_TO_PYTHON_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>

// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{

void AddConvergenceAcceleratorsToPython(pybind11::module &m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_CONVERGENCE_ACCELERATORS_TO_PYTHON_H_INCLUDED  defined 
