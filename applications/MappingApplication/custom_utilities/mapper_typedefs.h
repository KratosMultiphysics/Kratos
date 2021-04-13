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
#include "boost/numeric/ublas/vector.hpp"

// External includes

// Project includes
#include "spaces/ublas_space.h"


namespace Kratos
{
    namespace MapperDefinitions
    {
        typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;

        typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    }

}  // namespace Kratos.