//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

#if !defined(GET_SOLUTION_STEP_VALUES_FROM_GEOMETRY_UTILITIES_H_INCLUDED)
#define GET_SOLUTION_STEP_VALUES_FROM_GEOMETRY_UTILITIES_H_INCLUDED

#include "includes/model_part.h"

namespace Kratos
{
    typedef Element::GeometryType GeometryType;

    namespace GetSolutionStepValuesFromGeometryUtilities
    {

        /**
        */
        template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(
            const TVariableType& rThisVariable,
            const GeometryType& rGeometry,
            const IndexType IntegrationPointIndex)
        {
            const SizeType points_number = rShapeFunctionsValues.size2();

            typename TVariableType::Type rResult;
            const auto r_N = rGeometry.ShapeFunctionsValues();

            KRATOS_DEBUG_ERROR_IF(rGeometry.GetGeometry().size() < 1)
                << "Geometry size is smaller 1. FastGetSolutionStepValue cannot be computed." << std::endl;

            rResult = (*this)[0].FastGetSolutionStepValue(rThisVariable) * r_N(IntegrationPointIndex, 0);

            for (IndexType i = 1; i < points_number; ++i) {
                rResult += (*this)[i].FastGetSolutionStepValue(rThisVariable) * r_N(IntegrationPointIndex, i);
            }

            return rResult;
        }

        template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(
            const TVariableType& rThisVariable,
            const GeometryType& rGeometry,
            const IndexType IntegrationPointIndex)
        {
            const SizeType points_number = rShapeFunctionsValues.size2();

            typename TVariableType::Type rResult;
            const auto r_N = rGeometry.ShapeFunctionsValues();

            KRATOS_DEBUG_ERROR_IF(rGeometry.GetGeometry().size() < 1)
                << "Geometry size is smaller 1. FastGetSolutionStepValue cannot be computed." << std::endl;

            rResult = (*this)[0].GetSolutionStepValue(rThisVariable) * r_N(IntegrationPointIndex, 0);

            for (IndexType i = 1; i < points_number; ++i) {
                rResult += (*this)[i].FastGetSolutionStepValue(rThisVariable) * r_N(IntegrationPointIndex, i);
            }

            return rResult;
        }
    }

} // namespace Kratos

#endif // GET_SOLUTION_STEP_VALUES_FROM_GEOMETRY_UTILITIES_H_INCLUDED
