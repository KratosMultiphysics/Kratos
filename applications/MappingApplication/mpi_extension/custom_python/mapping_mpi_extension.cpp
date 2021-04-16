//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifdef KRATOS_PYTHON
// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_utilities/mapper_mpi_define.h"
#include "custom_python/add_mapper_to_python.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosMappingMPIExtension,m)
{
    auto py_mapper_factory = pybind11::class_<MapperFactory, MapperFactory::Pointer>(m, "MapperFactory");

    AddMapperToPython<MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType>(m, py_mapper_factory);
}

}
}

#endif // KRATOS_PYTHON defined
