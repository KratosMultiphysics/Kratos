//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//

#if !defined(KRATOS_STATISTICS_PYTHON_H_INCLUDED )
#define  KRATOS_STATISTICS_PYTHON_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define_python.h"


namespace Kratos {
namespace Python {

void  AddCustomStatisticsToPython(pybind11::module& m);

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_STATISTICS_PYTHON_H_INCLUDED  defined
