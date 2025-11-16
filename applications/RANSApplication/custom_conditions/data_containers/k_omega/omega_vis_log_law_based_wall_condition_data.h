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

#if !defined(KRATOS_K_OMEGA_WALL_CONDITION_DATA_OMEGA_VIS_LOG_WALL_BASED_CONDITION_DATA_H_INCLUDED)
#define KRATOS_K_OMEGA_WALL_CONDITION_DATA_OMEGA_VIS_LOG_WALL_BASED_CONDITION_DATA_H_INCLUDED

// System includes

// Project includes
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_conditions/scalar_wall_flux_condition_data.h"

namespace Kratos
{
///@name Kratos Classes
///@{

namespace KOmegaWallConditionData
{
class OmegaVisLogBasedWallConditionData : public ScalarWallFluxConditionData
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ScalarWallFluxConditionData;

    using GeometryType = typename BaseType::GeometryType;

    ///@}
    ///@name Static Operations
    ///@{

    static const Variable<double>& GetScalarVariable();

    static void Check(
        const Condition& rCondition,
        const ProcessInfo& rCurrentProcessInfo);

    static const std::string GetName() { return "KOmegaOmegaVisLogWallBasedConditionData"; }

    ///@}
    ///@name Life Cycle
    ///@{

    OmegaVisLogBasedWallConditionData(
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo)
    : BaseType(rGeometry, rProperties, rProcessInfo)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    const static bool IsWallFluxComputed(const GeometryType& rGeometry) {
        return (rGeometry.GetValue(RANS_IS_STRUCTURE) == 1);
    }

    void CalculateConstants(
        const ProcessInfo& rCurrentProcessInfo);

    double CalculateWallFlux(
        const Vector& rShapeFunctions,
        const ScalarWallFluxConditionData::Parameters& rParameters);

    ///@}

protected:
    ///@name Protected Members
    ///@{

    Vector mGaussPointWeights;
    Matrix mShapeFunctions;
    Geometry<Node>::ShapeFunctionsGradientsType mShapeFunctionDerivatives;

    double mWallDistance;
    double mOmegaViscous;
    double mOmegaLog;
    double mOmegaBlended;
    double mUTau;
    double mCmu25;
    double mKappa;

    ///@}
};

///@}

} // namespace KOmegaWallConditionData

} // namespace Kratos

#endif