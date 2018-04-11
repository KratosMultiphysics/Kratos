//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_ADD_TRILINOS_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_ADD_TRILINOS_COMMUNICATOR_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
{

void  AddTrilinosCommunicatorToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_TRILINOS_COMMUNICATOR_H_INCLUDED  defined
