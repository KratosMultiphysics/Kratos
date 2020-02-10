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
class IntegrationPointUtilities
{

private:
    ///@name private static members
    ///@{

    static std::vector<std::vector<std::vector<double>>> s_gauss_legendre;

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
        IntegrationPointsArrayType& rIntegrationPoints,
        SizeType PointsInU,
        double U0, double U1);

    static void IntegrationPoints2D(
        IntegrationPointsArrayType& rIntegrationPoints,
        SizeType PointsInU, SizeType PointsInV,
        double U0, double U1, double V0, double V1);

    ///@}
};
}
#endif // KRATOS_INTEGRATION_POINT_UTILITIES_INCLUDED  defined