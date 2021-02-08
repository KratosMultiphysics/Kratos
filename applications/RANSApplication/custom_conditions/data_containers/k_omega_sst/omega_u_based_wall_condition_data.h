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

#if !defined(KRATOS_K_OMEGA_SST_WALL_CONDITION_DATA_OMEGA_U_BASED_CONDITION_DATA_H_INCLUDED)
#define KRATOS_K_OMEGA_SST_WALL_CONDITION_DATA_OMEGA_U_BASED_CONDITION_DATA_H_INCLUDED

// System includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_conditions/scalar_wall_flux_condition_data.h"

namespace Kratos
{
///@name Kratos Classes
///@{

namespace KOmegaSSTWallConditionData
{
template<unsigned int TDim>
class OmegaUBasedWallConditionData : public ScalarWallFluxConditionData
{
public:
    using BaseType = ScalarWallFluxConditionData;
    using NodeType = Node<3>;
    using GeometryType = typename BaseType::GeometryType;
    using ShapeFunctionDerivativesArrayType = typename GeometryType::ShapeFunctionsGradientsType;

    static const Variable<double>& GetScalarVariable();

    static void Check(
        const Condition& rCondition,
        const ProcessInfo& rCurrentProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    static const std::string GetName()
    {
        return "KOmegaOmegaUBasedConditionData";
    }

    OmegaUBasedWallConditionData(
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo)
    : BaseType(rGeometry, rProperties, rProcessInfo)
    {
    }

    void CalculateConstants(
        const ProcessInfo& rCurrentProcessInfo);

    bool IsWallFluxComputable() const;

    double CalculateWallFlux(
        const Vector& rShapeFunctions);

protected:
    double mKappa;
    double mInvKappa;
    double mBeta;
    double mYPlus;
    double mCmu25;
    double mDensity;
    double mTurbulentViscosity;
    double mSigmaOmega1;
    double mSigmaOmega2;
    double mBetaStar;
    double mWallHeight;
};

///@}

} // namespace KOmegaSSTWallConditionData

} // namespace Kratos

#endif