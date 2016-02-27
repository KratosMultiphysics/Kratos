//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_ADD_CUSTOM_DEM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_DEM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED
// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h"  

// External includes
#include "boost/smart_ptr.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomConstitutiveLawsToPython();

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_DEM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED  defined 
