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

#pragma once

// System includes

// External includes
#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"

// Project includes
#include "trilinos_space.h"
#include "spaces/ublas_space.h"

namespace Kratos {

namespace MPIMapperDefinitions {

    typedef TUblasDenseSpace<double> DenseSpaceType;

    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> SparseSpaceType;

}  // namespace MPIMapperDefinitions.

}  // namespace Kratos.
