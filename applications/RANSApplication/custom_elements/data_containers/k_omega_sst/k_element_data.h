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
class KElementData : public ConvectionDiffusionReactionElementData<TDim>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ConvectionDiffusionReactionElementData<TDim>;

    using NodeType = Node;

    using GeometryType = typename BaseType::GeometryType;

    using ArrayD = typename BaseType::ArrayD;

    ///@}
    ///@name Life Cycle
    ///@{

    KElementData(
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo)
        : BaseType(rGeometry, rProperties, rProcessInfo)
    {
    }

    ~KElementData() override = default;

    ///@}
    ///@name Static Operations
    ///@{

    static const Variable<double>& GetScalarVariable();

    static void Check(
        const Element& rElement,
        const ProcessInfo& rCurrentProcessInfo);

    static const std::string GetName()
    {
        return "KOmegaSSTKElementData";
    }

    ///@}
    ///@name Operations
    ///@{

    void CalculateConstants(
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateGaussPointData(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const int Step = 0);

    ///@}

protected:
    ///@name Protected Members
    ///@{

    using BaseType::mEffectiveVelocity;
    using BaseType::mEffectiveKinematicViscosity;
    using BaseType::mReactionTerm;
    using BaseType::mSourceTerm;

    BoundedMatrix<double, TDim, TDim> mVelocityGradient;
    ArrayD mTurbulentKineticEnergyGradient;
    ArrayD mTurbulentSpecificEnergyDissipationRateGradient;

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
    double mA1;

    ///@}
    ///@name Protected Operations
    ///@{

    double CalculateEffectiveViscosity(
        const ProcessInfo& rCurrentProcessInfo);

    ///@}
};

///@}

} // namespace KOmegaSSTElementData

} // namespace Kratos

#endif