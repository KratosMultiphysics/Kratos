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

#if !defined(KRATOS_TESTING_UTILITIES_H_INCLUDED)
#define  KRATOS_TESTING_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/model_partio.h"

namespace Kratos
{
namespace MapperTestingUtilities
{

// This function creates a ModelPart
// If executed in MPI, it will create a fully functional distributed ModelPart
// The prerequisite for this is that the partitioned files exist!
// it returns a boolena whether the generation of the ModelPart has been successful
bool CreateModelPart(ModelPart& rModelPart)
{
    if num_nodes_global > 0
        return false;



    ModelPartIO model_part_io(p_input);





    return false;
}


}  // namespace MapperUtilities.

}  // namespace Kratos.

#endif // KRATOS_TESTING_UTILITIES_H_INCLUDED  defined





