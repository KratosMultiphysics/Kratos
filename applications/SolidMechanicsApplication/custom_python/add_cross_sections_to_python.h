//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:       Massimo Petracca $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                     2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(ADD_CROSS_SECTIONS_TO_PYTHON_H_INCLUDED)
#define ADD_CROSS_SECTIONS_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos
{

namespace Python
{

void AddCrossSectionsToPython(pybind11::module& m);

}

}


#endif // ADD_CROSS_SECTIONS_TO_PYTHON_H_INCLUDED
