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
#include "mpi.h"

// Project includes
#include "mapper_utilities_mpi.h"
#include "mapper_utilities.h"

namespace Kratos
{

namespace MapperUtilitiesMPI
{

void ComputeGlobalBoundingBoxes(ModelPart& rModelPart, std::vector<double>& rGlobalBoundingBoxes)
{
    const auto local_bounding_box = MapperUtilities::ComputeLocalBoundingBox(rModelPart); // TODO is const ok here?

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    rGlobalBoundingBoxes.resize(6*comm_size);

    MPI_Allgather(local_bounding_box.data(), 6, MPI_DOUBLE, rGlobalBoundingBoxes.data(), 6, MPI_DOUBLE, MPI_COMM_WORLD); // TODO working?
}

}  // namespace MapperUtilitiesMPI.


}  // namespace Kratos.
