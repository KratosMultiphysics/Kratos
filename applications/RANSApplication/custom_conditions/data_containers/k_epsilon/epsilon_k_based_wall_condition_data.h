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

namespace KEpsilonWallConditionData
{
class EpsilonKBasedWallConditionData : public ScalarWallFluxConditionData
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ScalarWallFluxConditionData;

    using GeometryType = BaseType::GeometryType;

    ///@}
    ///@name Static Operations
    ///@{

    static const Variable<double>& GetScalarVariable();

    static void Check(
        const Condition& rCondition,
        const ProcessInfo& rCurrentProcessInfo);

    static const std::string GetName() { return "KEpsilonEpsilonKBasedConditionData"; }

    ///@}
    ///@name Life Cycle
    ///@{

    EpsilonKBasedWallConditionData(
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo)
    : BaseType(rGeometry, rProperties, rProcessInfo)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void CalculateConstants(
        const ProcessInfo& rCurrentProcessInfo);

    double CalculateWallFlux(
        const Vector& rShapeFunctions,
        const ScalarWallFluxConditionData::Parameters& rParameters);

    ///@}

protected:
    ///@name Protected Members
    ///@{

    double mEpsilonSigma;
    double mCmu25;

    ///@}
};

///@}

} // namespace KEpsilonWallConditionData

} // namespace Kratos