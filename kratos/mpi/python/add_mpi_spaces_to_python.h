
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:     Philipp Bucher (https://github.com/philbucher)
//

#ifndef KRATOS_ADD_MPI_SPACES_TO_PYTHON_H_INCLUDED
#define KRATOS_ADD_MPI_SPACES_TO_PYTHON_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define_python.h"


namespace Kratos {
namespace Python {

void AddMPISpacesToPython(pybind11::module& m);

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_ADD_MPI_SPACES_TO_PYTHON_H_INCLUDED defined
