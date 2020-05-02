//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_EVM_K_EPSILON_WALL_CONDITION_DATA_UTILITIES_EPSILON_K_BASED_CONDITION_DATA_H_INCLUDED)
#define KRATOS_RANS_EVM_K_EPSILON_WALL_CONDITION_DATA_UTILITIES_EPSILON_K_BASED_CONDITION_DATA_H_INCLUDED

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
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

namespace EvmKEpsilonWallConditionDataUtilities
{
class EpsilonKBasedWallConditionData : public ScalarWallFluxConditionData
{
public:
    using BaseType = ScalarWallFluxConditionData;
    using NodeType = Node<3>;
    using GeomtryType = BaseType::GeometryType;

    static const Variable<double>& GetScalarVariable();
    static const Variable<double>& GetScalarRateVariable();

    static void Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    static const std::string GetName() {return "KEpsilonEpsilonKBasedConditionData";}

    EpsilonKBasedWallConditionData(const GeomtryType& rGeometry) : BaseType(rGeometry)
    {
    }

    void CalculateConstants(const ProcessInfo& rCurrentProcessInfo) override;

    bool IsWallFluxComputable() const override;

    double CalculateWallFlux(const Vector& rShapeFunctions) const override;

protected:
    double mEpsilonSigma;
    double mKappa;
    double mYPlus;
    double mCmu25;
};

///@}

} // namespace EvmKEpsilonWallConditionDataUtilities

} // namespace Kratos

#endif