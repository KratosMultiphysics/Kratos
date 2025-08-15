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
#include "includes/node.h"
#include "spaces/ublas_space.h"
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
                                                   const IndexType      BufferIndex);

    static void SetUFirstAndSecondDerivativeVector(const SystemVectorType& rFirstDerivativeVector,
                                                   const SystemVectorType& rSecondDerivativeVector,
                                                   ModelPart&              rModelPart);

private:
    static void GetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                                  const Node&             rNode,
                                                  SystemVectorType&       rFirstDerivativeVector,
                                                  SystemVectorType&       rSecondDerivativeVector,
                                                  const IndexType         BufferIndex);

    static void SetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                                  Node&                   rNode,
                                                  const SystemVectorType& rFirstDerivativeVector,
                                                  const SystemVectorType& rSecondDerivativeVector);

    static void GetDerivativesForVariable(const Variable<double>& rVariable,
                                          const Node&             rNode,
                                          SystemVectorType&       rFirstDerivativeVector,
                                          SystemVectorType&       rSecondDerivativeVector,
                                          const IndexType         BufferIndex);

    static void SetDerivativesForVariable(const Variable<double>& rVariable,
                                          Node&                   rNode,
                                          const SystemVectorType& rFirstDerivativeVector,
                                          const SystemVectorType& rSecondDerivativeVector);
};
} // namespace Kratos::Geo
