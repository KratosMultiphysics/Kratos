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
void KRATOS_API(MAPPING_APPLICATION) InitializeSystemVector(
    Kratos::unique_ptr<typename TSparseSpace::VectorType>& rpVector,
    const std::size_t VectorSize);

template<class TSparseSpace, class TDenseSpace>
void KRATOS_API(MAPPING_APPLICATION) BuildMappingMatrix(
    Kratos::unique_ptr<typename TSparseSpace::MatrixType>& rpMappingMatrix,
    Kratos::unique_ptr<typename TSparseSpace::VectorType>& rpInterfaceVectorOrigin,
    Kratos::unique_ptr<typename TSparseSpace::VectorType>& rpInterfaceVectorDestination,
    const ModelPart& rModelPartOrigin,
    const ModelPart& rModelPartDestination,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
    const int EchoLevel);

template<class TSparseSpace, class TDenseSpace>
void KRATOS_API(MAPPING_APPLICATION) CheckRowSum(
    const typename TSparseSpace::MatrixType& rM,
    const std::string& rBaseFileName,
    const bool ThrowError = false);

}  // namespace MappinMatrixUtilities.

}  // namespace Kratos.

#endif // KRATOS_MAPPING_MATRIX_UTILITIES_H_INCLUDED  defined
