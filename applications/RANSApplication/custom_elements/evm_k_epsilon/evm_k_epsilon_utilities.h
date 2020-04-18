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

#if !defined(KRATOS_RANS_EVM_K_EPSILON_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_EVM_K_EPSILON_UTILITIES_H_INCLUDED

// System includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_element_data.h"

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

namespace EvmKepsilonModelUtilities
{
/// Node type
using NodeType = Node<3>;
using GeometryType = Geometry<NodeType>;

template <unsigned int TDim>
class EpsilonElementData : public ScalarConvectionDiffusionReactionElementData
{
public:
    static const Variable<double>& GetScalarVariable();
    static const Variable<double>& GetScalarRateVariable();
    static const Variable<double>& GetScalarRelaxedRateVariable();

    static void Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    void CalculateElementData(const GeometryType& rGeometry,
                              const Vector& rShapeFunctions,
                              const Matrix& rShapeFunctionDerivatives,
                              const ProcessInfo& rCurrentProcessInfo,
                              const int Step = 0) override;

    double CalculateEffectiveKinematicViscosity(const Vector& rShapeFunctions,
                                                const Matrix& rShapeFunctionDerivatives,
                                                const ProcessInfo& rCurrentProcessInfo) const override;

    double CalculateReactionTerm(const Vector& rShapeFunctions,
                                 const Matrix& rShapeFunctionDerivatives,
                                 const ProcessInfo& rCurrentProcessInfo) const override;

    double CalculateSourceTerm(const Vector& rShapeFunctions,
                               const Matrix& rShapeFunctionDerivatives,
                               const ProcessInfo& rCurrentProcessInfo) const override;

private:
    BoundedMatrix<double, TDim, TDim> mVelocityGradient;

    double mC1;
    double mC2;
    double mGamma;
    double mTurbulentKineticEnergy;
    double mTurbulentKinematicViscosity;
    double mKinematicViscosity;
    double mVelocityDivergence;
};

template <unsigned int TDim>
class KElementData : public ScalarConvectionDiffusionReactionElementData
{
public:
    static const Variable<double>& GetScalarVariable();
    static const Variable<double>& GetScalarRateVariable();
    static const Variable<double>& GetScalarRelaxedRateVariable();

    static void Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    void CalculateElementData(const GeometryType& rGeometry,
                              const Vector& rShapeFunctions,
                              const Matrix& rShapeFunctionDerivatives,
                              const ProcessInfo& rCurrentProcessInfo,
                              const int Step = 0) override;

    double CalculateEffectiveKinematicViscosity(const Vector& rShapeFunctions,
                                                const Matrix& rShapeFunctionDerivatives,
                                                const ProcessInfo& rCurrentProcessInfo) const override;

    double CalculateReactionTerm(const Vector& rShapeFunctions,
                                 const Matrix& rShapeFunctionDerivatives,
                                 const ProcessInfo& rCurrentProcessInfo) const override;

    double CalculateSourceTerm(const Vector& rShapeFunctions,
                               const Matrix& rShapeFunctionDerivatives,
                               const ProcessInfo& rCurrentProcessInfo) const override;

private:
    BoundedMatrix<double, TDim, TDim> mVelocityGradient;

    double mGamma;
    double mTurbulentKineticEnergy;
    double mTurbulentKinematicViscosity;
    double mKinematicViscosity;
    double mVelocityDivergence;
};

double CalculateTurbulentViscosity(const double C_mu,
                                   const double turbulent_kinetic_energy,
                                   const double turbulent_energy_dissipation_rate,
                                   const double f_mu);

double CalculateFmu(const double y_plus);

double CalculateF2(const double turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double turbulent_energy_dissipation_rate);

template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                           const double turbulent_kinematic_viscosity,
                           const double turbulent_kinetic_energy);

double CalculateGamma(const double C_mu,
                      const double f_mu,
                      const double turbulent_kinetic_energy,
                      const double turbulent_kinematic_viscosity);

} // namespace EvmKepsilonModelUtilities

///@}

} // namespace Kratos

#endif