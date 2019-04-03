


#if !defined(KRATOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED

// System includes
// External includes 
//#include "boost/smart_ptr.hpp"
#include <pybind11/pybind11.h>

#include "includes/define_python.h"

namespace Kratos
{

	namespace Python
	{

		void AddCustomProcessesToPython(pybind11::module& m);

	}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED  defined 
