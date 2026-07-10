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
    static const std::vector<std::vector<std::array<double, 3>>> s_gauss_triangle;

    ///@}
public:

    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef IntegrationPoint<3> IntegrationPointType;
    typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;

    ///@}
    ///@name Generic generation of Integration Points
    ///@{
    /// 
    static void CreateIntegrationPoints1D(
        IntegrationPointsArrayType & rIntegrationPoints,
        const std::vector<double>&rSpansLocalSpace,
        const IntegrationInfo& rIntegrationInfo);
    
    static void CreateIntegrationPoints2D(
        IntegrationPointsArrayType & rIntegrationPoints,
        const std::vector<double>&rSpansLocalSpaceU,
        const std::vector<double>&rSpansLocalSpaceV,
        const IntegrationInfo& rIntegrationInfo);

    ///@}
    ///@name Create Integration Points
    ///@{

    static void CreateIntegrationPoints1DGauss(
        IntegrationPointsArrayType& rIntegrationPoints,
        const std::vector<double>& rSpanIntervals,
        const SizeType IntegrationPointsPerSpan);
    
    static void CreateIntegrationPoints2DGauss(
        IntegrationPointsArrayType& rIntegrationPoints,
        const std::vector<double>& rSpanIntervalsU,
        const std::vector<double>& rSpanIntervalsV,
        const SizeType IntegrationPointsPerSpanU,
        const SizeType IntegrationPointsPerSpanV);

    static void CreateIntegrationPoints1DGrid(
        IntegrationPointsArrayType & rIntegrationPoints,
        const std::vector<double>& rSpanIntervals,
        const SizeType IntegrationPointsPerSpan);

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

    static void IntegrationPoints3D(
        typename IntegrationPointsArrayType::iterator& rIntegrationPointsBegin,
        SizeType PointsInU, SizeType PointsInV, SizeType PointsInW,
        double U0, double U1, double V0, double V1, double W0, double W1);

    /// Triangular shape
    static void IntegrationPointsTriangle2D(
        typename IntegrationPointsArrayType::iterator & rIntegrationPointsBegin,
        SizeType PointsIndex,
        double U0, double U1, double U2, double V0, double V1, double V2);

    ///@}
};

}
#endif // KRATOS_INTEGRATION_POINT_UTILITIES_INCLUDED  defined
