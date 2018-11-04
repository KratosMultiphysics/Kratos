//
//   Project Name:        Kratos
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-07-23 14:47:45 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_ADD_LINEAR_SOLVERS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_LINEAR_SOLVERS_TO_PYTHON_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos {
namespace Python {

void  AddLinearSolversToPython(pybind11::module& m);

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_ADD_LINEAR_SOLVERS_TO_PYTHON_H_INCLUDED  defined
