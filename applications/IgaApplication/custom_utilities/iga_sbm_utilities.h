//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_IGA_SBM_UTILITIES_H_INCLUDED)
#define  KRATOS_IGA_SBM_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "iga_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) IgaSbmUtilities
{
public:
    using IndexType = std::size_t;
    using SizeType = std::size_t;
    using GeometryType = Geometry<Node>;

    static void GetDeformedPosition(
        const Condition& rCondition,
        array_1d<double, 3>& rPointDeformedCoordinates);

    static void GetSolutionCoefficientVector(
        const Variable<array_1d<double, 3>>& rVariable,
        const GeometryType& rReferenceGeometry,
        Vector& rValues);
};

///@}
}  // namespace Kratos.

#endif // KRATOS_IGA_SBM_UTILITIES_H_INCLUDED
