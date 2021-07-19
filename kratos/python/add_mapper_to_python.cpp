//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "mappers/mapper_define.h"
#include "add_mapper_to_python.h"

namespace Kratos {
namespace Python {

void AddMapperToPython(pybind11::module& m)
{
    AddMappingToPython<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>(m);
}

}  // namespace Python.
}  // namespace Kratos.
