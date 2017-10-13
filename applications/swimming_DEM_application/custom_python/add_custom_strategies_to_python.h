//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Guillermo Casas, gcasas@cimne.upc.edu  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $


#if !defined(KRATOS_STRATEGIES_PYTHON_H_INCLUDED )
#define  KRATOS_STRATEGIES_PYTHON_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"

namespace Kratos
{
    namespace Python
    {
        using namespace boost::python;

        void  AddCustomStrategiesToPython();

    }  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_STRATEGIES_PYTHON_H_INCLUDED  defined 
