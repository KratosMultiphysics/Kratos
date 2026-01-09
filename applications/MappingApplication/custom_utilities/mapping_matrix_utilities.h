//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mapper_local_system.h"

namespace Kratos {

template<class TSparseSpace, class TDenseSpace>
class KRATOS_API(MAPPING_APPLICATION) MappingMatrixUtilities
{
public:
static void InitializeSystemVector(
    Kratos::unique_ptr<typename TSparseSpace::VectorType>& rpVector,
    const std::size_t VectorSize);

static void BuildMappingMatrix(
    Kratos::unique_ptr<typename TSparseSpace::MatrixType>& rpMappingMatrix,
    Kratos::unique_ptr<typename TSparseSpace::VectorType>& rpInterfaceVectorOrigin,
    Kratos::unique_ptr<typename TSparseSpace::VectorType>& rpInterfaceVectorDestination,
    const ModelPart& rModelPartOrigin,
    const ModelPart& rModelPartDestination,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
    const int EchoLevel);

static void BuildMappingMatrixRBFMapper(
    Kratos::unique_ptr<typename TSparseSpace::MatrixType>& rpMappingMatrix,
    Kratos::unique_ptr<typename TSparseSpace::VectorType>& rpInterfaceVectorOrigin,
    Kratos::unique_ptr<typename TSparseSpace::VectorType>& rpInterfaceVectorDestination,
    const ModelPart& rModelPartOrigin,
    const ModelPart& rModelPartDestination,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
    const IndexType NumberOfPolynomialTerms,
    const bool BuildOriginInterpolationMatrix,
    const bool OriginIsIga,
    const int EchoLevel);

static void CheckRowSum(
    const typename TSparseSpace::MatrixType& rM,
    const std::string& rBaseFileName,
    const bool ThrowError = false,
    const double Tolerance = 1e-15);

};

}  // namespace Kratos.
