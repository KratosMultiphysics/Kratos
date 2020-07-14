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

#if !defined(KRATOS_MAPPING_MATRIX_UTILITIES_H_INCLUDED)
#define  KRATOS_MAPPING_MATRIX_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mapper_local_system.h"

namespace Kratos
{
namespace MappingMatrixUtilities
{

template<class TSparseSpace, class TDenseSpace>
void KRATOS_API(MAPPING_APPLICATION) BuildMappingMatrix(
    Kratos::unique_ptr<typename TSparseSpace::MatrixType>& rpMappingMatrix,
    Kratos::unique_ptr<typename TSparseSpace::VectorType>& rpInterfaceVectorOrigin,
    Kratos::unique_ptr<typename TSparseSpace::VectorType>& rpInterfaceVectorDestination,
    const ModelPart& rModelPartOrigin,
    const ModelPart& rModelPartDestination,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
    const int EchoLevel);

}  // namespace MappinMatrixUtilities.

}  // namespace Kratos.

#endif // KRATOS_MAPPING_MATRIX_UTILITIES_H_INCLUDED  defined
