//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah
//                   Bence Rochlitz
//
//  Supervised by:   Jordi Cotela
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_K_OMEGA_ELEMENT_DATA_OMEGA_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_K_OMEGA_ELEMENT_DATA_OMEGA_ELEMENT_DATA_H_INCLUDED

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
///@name  Functions
///@{

namespace KOmegaElementData
{
template <unsigned int TDim>
class OmegaElementData : public ConvectionDiffusionReactionElementData<TDim>
{
public:
    using BaseType = ConvectionDiffusionReactionElementData<TDim>;
    using NodeType = Node<3>;
    using GeometryType = typename BaseType::GeometryType;

    static const Variable<double>& GetScalarVariable();

    static void Check(
        const GeometryType& rGeometry,
        const ProcessInfo& rCurrentProcessInfo);

    static const std::string GetName()
    {
        return "KOmegaOmegaElementData";
    }

    OmegaElementData(
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo)
        : BaseType(rGeometry, rProperties, rProcessInfo)
    {
    }

    void CalculateConstants(
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateGaussPointData(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const int Step = 0);

protected:
    BoundedMatrix<double, TDim, TDim> mVelocityGradient;

    double mTurbulentKineticEnergy;
    double mTurbulentKinematicViscosity;
    double mKinematicViscosity;
    double mVelocityDivergence;
    double mSigmaOmega;
    double mBeta;
    double mGamma;
    double mDensity;
};

///@}

} // namespace KOmegaElementData

} // namespace Kratos

#endif