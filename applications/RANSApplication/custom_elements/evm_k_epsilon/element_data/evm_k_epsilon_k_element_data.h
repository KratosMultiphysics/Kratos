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

#if !defined(KRATOS_RANS_EVM_K_EPSILON_ELEMENT_DATA_UTILITIES_K_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_RANS_EVM_K_EPSILON_ELEMENT_DATA_UTILITIES_K_ELEMENT_DATA_H_INCLUDED

// System includes

// Project includes
#include "containers/variable.h"
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

namespace EvmKEpsilonElementDataUtilities
{
template <unsigned int TDim>
class KElementData : public ScalarConvectionDiffusionReactionElementData
{
public:
    using BaseType = ScalarConvectionDiffusionReactionElementData;
    using NodeType = Node<3>;
    using GeomtryType = BaseType::GeometryType;

    static const Variable<double>& GetScalarVariable();
    static const Variable<double>& GetScalarRateVariable();
    static const Variable<double>& GetScalarRelaxedRateVariable();

    static void Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    KElementData(const GeomtryType& rGeometry) : BaseType(rGeometry)
    {
    }

    void CalculateConstants(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateGaussPointData(const Vector& rShapeFunctions,
                                 const Matrix& rShapeFunctionDerivatives,
                                 const int Step = 0) override;

    double CalculateEffectiveKinematicViscosity(const Vector& rShapeFunctions,
                                                const Matrix& rShapeFunctionDerivatives) const override;

    double CalculateReactionTerm(const Vector& rShapeFunctions,
                                 const Matrix& rShapeFunctionDerivatives) const override;

    double CalculateSourceTerm(const Vector& rShapeFunctions,
                               const Matrix& rShapeFunctionDerivatives) const override;

protected:
    BoundedMatrix<double, TDim, TDim> mVelocityGradient;

    double mGamma;
    double mTurbulentKinematicViscosity;
    double mKinematicViscosity;
    double mVelocityDivergence;
    double mInvTkeSigma;
    double mCmu;
};

///@}

} // namespace EvmKEpsilonElementDataUtilities

} // namespace Kratos

#endif