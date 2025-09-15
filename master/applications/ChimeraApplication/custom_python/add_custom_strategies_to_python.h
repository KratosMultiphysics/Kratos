//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
//

#if !defined(KRATOS_STRATEGIES_PYTHON_H_INCLUDED )
#define  KRATOS_STRATEGIES_PYTHON_H_INCLUDED


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"


namespace Kratos
{

    namespace Python
    {

      void  AddCustomStrategiesToPython(pybind11::module& m);


    }  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_STRATEGIES_PYTHON_H_INCLUDED  defined
