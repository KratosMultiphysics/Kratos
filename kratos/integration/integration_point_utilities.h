//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#if !defined(KRATOS_INTEGRATION_POINT_UTILITIES_INCLUDED )
#define  KRATOS_INTEGRATION_POINT_UTILITIES_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "integration/integration_point.h"
#include "geometries/geometry.h"

namespace Kratos
{
class KRATOS_API(KRATOS_CORE) IntegrationPointUtilities
{

public:
    ///@name private static members
    ///@{

    static const std::vector<std::vector<std::array<double, 2>>> s_gauss_legendre;

    ///@}
public:

    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef IntegrationPoint<3> IntegrationPointType;
    typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;

    ///@}
    ///@name Define Integration Points
    ///@{

    static void IntegrationPoints1D(
        typename IntegrationPointsArrayType::iterator& rIntegrationPointsBegin,
        SizeType PointsInU,
        double U0, double U1);

    static void IntegrationPoints2D(
        typename IntegrationPointsArrayType::iterator& rIntegrationPointsBegin,
        SizeType PointsInU, SizeType PointsInV,
        double U0, double U1, double V0, double V1);

    ///@}
};

}
#endif // KRATOS_INTEGRATION_POINT_UTILITIES_INCLUDED  defined