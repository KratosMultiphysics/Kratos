// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ADD_CUSTOM_PROCESSES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_PROCESSES_TO_PYTHON_H_INCLUDED

// System includes 
#include <pybind11/pybind11.h>

// External includes 

// Project includes
#include "includes/define.h"


namespace Kratos
{

  namespace Python
  {

    void  AddCustomProcessesToPython(pybind11::module& m);

  }  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_PROCESSES_TO_PYTHON_H_INCLUDED  defined 
