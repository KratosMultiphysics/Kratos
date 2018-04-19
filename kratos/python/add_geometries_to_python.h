//    |  /           |                 
//    ' /   __| _` | __|  _  \  __| 
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_ADD_GEOMETRIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_GEOMETRIES_TO_PYTHON_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>

// External includes


// Project includes


namespace Kratos
{

namespace Python
{

void  AddGeometriesToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_GEOMETRIES_TO_PYTHON_H_INCLUDED  defined 
