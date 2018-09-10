// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_ADD_CUSTOM_IO_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_IO_TO_PYTHON_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>

// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{

void  AddCustomIOToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_IO_TO_PYTHON_H_INCLUDED  defined 
