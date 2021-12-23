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

#if !defined(KRATOS_K_OMEGA_SST_ELEMENT_DATA_K_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_K_OMEGA_SST_ELEMENT_DATA_K_ELEMENT_DATA_H_INCLUDED

// System includes

// Project includes
#include "containers/array_1d.h"
#include "containers/variable.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_element_data.h"

namespace Kratos
{
///@name  Functions
///@{

namespace KOmegaSSTElementData
{
template <unsigned int TDim>
class KElementData : public ConvectionDiffusionReactionElementData
{
public:
    using BaseType = ConvectionDiffusionReactionElementData;
    using NodeType = Node<3>;
    using GeomtryType = BaseType::GeometryType;

    static const Variable<double>& GetScalarVariable();

    static void Check(
        const GeometryType& rGeometry,
        const ProcessInfo& rCurrentProcessInfo);

    static const std::string GetName()
    {
        return "KOmegaSSTKElementData";
    }

    KElementData(
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo)
        : BaseType(rGeometry, rProperties, rProcessInfo)
    {
    }

    void Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateConstants(
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateGaussPointData(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const int Step = 0);

    array_1d<double, 3> CalculateEffectiveVelocity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const;

    double CalculateEffectiveKinematicViscosity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const;

    double CalculateReactionTerm(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const;

    double CalculateSourceTerm(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const;

    double CalculateEffectiveViscosity(
        const ProcessInfo& rCurrentProcessInfo);

protected:
    BoundedMatrix<double, TDim, TDim> mVelocityGradient;
    array_1d<double, 3> mEffectiveVelocity;
    array_1d<double, 3> mTurbulentKineticEnergyGradient;
    array_1d<double, 3> mTurbulentSpecificEnergyDissipationRateGradient;

    double mSigmaK1;
    double mSigmaK2;
    double mSigmaOmega2;
    double mBetaStar;
    double mTurbulentKineticEnergy;
    double mTurbulentSpecificEnergyDissipationRate;
    double mKinematicViscosity;
    double mTurbulentKinematicViscosity;
    double mWallDistance;
    double mCrossDiffusion;
    double mBlendedSimgaK;
    double mVelocityDivergence;
    double mDensity;
};

///@}

} // namespace KOmegaSSTElementData

} // namespace Kratos

#endif