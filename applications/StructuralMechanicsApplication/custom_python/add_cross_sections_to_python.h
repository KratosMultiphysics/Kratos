// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(ADD_CROSS_SECTIONS_TO_PYTHON_H_INCLUDED)
#define ADD_CROSS_SECTIONS_TO_PYTHON_H_INCLUDED

#include <pybind11/pybind11.h>

namespace Kratos
{

namespace Python
{

void AddCrossSectionsToPython(pybind11::module& m);

}

}


#endif // ADD_CROSS_SECTIONS_TO_PYTHON_H_INCLUDED
