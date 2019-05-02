//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_ADD_CUSTOM_DEM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_DEM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED

// System includes

// Project includes
#include "includes/define.h"  

// External includes
#include <pybind11/pybind11.h>

namespace Kratos
{

namespace Python
{

void  AddCustomConstitutiveLawsToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_DEM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED  defined 
