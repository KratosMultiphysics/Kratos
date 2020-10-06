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

#if !defined(KRATOS_K_OMEGA_WALL_CONDITION_DATA_OMEGA_U_BASED_CONDITION_DATA_H_INCLUDED)
#define KRATOS_K_OMEGA_WALL_CONDITION_DATA_OMEGA_U_BASED_CONDITION_DATA_H_INCLUDED

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

namespace KOmegaWallConditionData
{
class OmegaUBasedWallConditionData : public ScalarWallFluxConditionData
{
public:
    using BaseType = ScalarWallFluxConditionData;
    using NodeType = Node<3>;
    using GeomtryType = BaseType::GeometryType;

    static const Variable<double>& GetScalarVariable();

    static void Check(
        const GeometryType& rGeometry,
        const ProcessInfo& rCurrentProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    static const std::string GetName()
    {
        return "KOmegaOmegaUBasedConditionData";
    }

    OmegaUBasedWallConditionData(
        const GeomtryType& rGeometry)
    : BaseType(rGeometry)
    {
    }

    void CalculateConstants(
        const ProcessInfo& rCurrentProcessInfo) override;

    bool IsWallFluxComputable() const override;

    double CalculateWallFlux(
        const Vector& rShapeFunctions) const override;

protected:
    double mOmegaSigma;
    double mKappa;
    double mInvKappa;
    double mBeta;
    double mYPlus;
    double mCmu25;
};

///@}

} // namespace KOmegaWallConditionData

} // namespace Kratos

#endif