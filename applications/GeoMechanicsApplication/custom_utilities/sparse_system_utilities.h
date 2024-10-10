// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once

#include <algorithm>
#include <vector>

// External includes
#include "containers/system_vector.h"
#include "geo_aliases.h"
#include "geometries/geometry.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "solving_strategies/schemes/scheme.h"
#include "spaces/ublas_space.h"

#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

namespace Kratos::Geo
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SparseSystemUtilities
{
public:
    /// Vector type definition
    using SystemVectorType = UblasSpace<double, CompressedMatrix, Vector>::VectorType;
    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;


    static void GetUFirstAndSecondDerivativeVector(SystemVectorType&    rFirstDerivativeVector,
                                                   SystemVectorType&    rSecondDerivativeVector,
                                                   const DofsArrayType& rDofSet,
                                                   const ModelPart&     rModelPart,
                                                   const IndexType      bufferIndex);

    static void SetUFirstAndSecondDerivativeVector(const SystemVectorType& rFirstDerivativeVector,
                                                   const SystemVectorType& rSecondDerivativeVector,
                                                   ModelPart&              rModelPart);

private:
    static void GetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                                  const Node&             rNode,
                                                  SystemVectorType&       rFirstDerivativeVector,
                                                  SystemVectorType&       rSecondDerivativeVector,
                                                  const IndexType         bufferIndex);

    static void SetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                                  Node&                   rNode,
                                                  const SystemVectorType& rFirstDerivativeVector,
                                                  const SystemVectorType& rSecondDerivativeVector);

    static void GetDerivativesForVariable(const Variable<double>& rVariable,
                                          const Node&             rNode,
                                          SystemVectorType&       rFirstDerivativeVector,
                                          SystemVectorType&       rSecondDerivativeVector,
                                          const IndexType         bufferIndex);

    static void SetDerivativesForVariable(const Variable<double>& rVariable,
                                          Node&                   rNode,
                                          const SystemVectorType& rFirstDerivativeVector,
                                          const SystemVectorType& rSecondDerivativeVector);
};
} // namespace Kratos::Geo
