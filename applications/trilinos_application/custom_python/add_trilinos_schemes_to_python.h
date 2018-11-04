

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-10-23 12:48:28 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_TRILINOS_ADD_SCHEMES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_TRILINOS_ADD_SCHEMES_TO_PYTHON_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"



namespace Kratos
{
namespace Python
{
void  AddSchemes(pybind11::module& m);
}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_TRILINOS_ADD_SCHEMES_TO_PYTHON_H_INCLUDED  defined 
