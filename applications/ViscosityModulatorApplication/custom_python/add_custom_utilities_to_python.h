// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
// | | | (_)        |  \/  |         | |
// | | | |_ ___  ___| .  . | ___   __| |
// | | | | / __|/ __| |\/| |/ _ \ / _` |
// \ \_/ / \__ \ (__| |  | | (_) | (_| |
//  \___/|_|___/\___\_|  |_/\___/ \__,_|  APPLICATION
//                                      
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala
//

#if !defined(KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED  defined 
