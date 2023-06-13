//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"

// Application includes
#include "collective_expression.h"

namespace Kratos {

class KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpressionIO
{
public:
    ///@name Public static operations
    ///@{

    template<class TRawDataType>
    static void ReadCArray(
        CollectiveExpression& rCollectiveExpression,
        TRawDataType const* pBegin,
        int const* NumberOfEntities,
        int const** pListShapeBegin,
        int const* ShapeSizes,
        const int NumberOfContainers);

    template<class TRawDataType>
    static void MoveCArray(
        CollectiveExpression& rCollectiveExpression,
        TRawDataType* pBegin,
        int const* NumberOfEntities,
        int const** pListShapeBegin,
        int const* ShapeSizes,
        const int NumberOfContainers);

    template<class TRawDataType>
    static void WriteCArray(
        const CollectiveExpression& rCollectiveExpression,
        TRawDataType* pBegin,
        const int Size);

    ///@}
};

} // namespace Kratos
