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

#if !defined(KRATOS_MAPPER_DEFINE_H_INCLUDED)
#define KRATOS_MAPPER_DEFINE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "spaces/default_spaces.h"

namespace Kratos {

namespace MapperDefinitions {

    typedef TDefaultDenseSpace<double> DenseSpaceType;

    typedef TDefaultSparseSpace<double> SparseSpaceType;

}  // namespace MapperDefinitions.

}  // namespace Kratos.

#endif // KRATOS_MAPPER_DEFINE_H_INCLUDED defined
