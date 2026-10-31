//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

// System includes

// External includes

// Project includes
#include "custom_elements/gsbm_laplacian_element.h"
#include "iga_application_variables.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{

namespace
{

using IndexType = std::size_t;

/// Consistent mass \sum_{gp} w N \otimes N, integrated as LaplacianElement does
Matrix ComputeMassMatrix(const Element& rElement)
{
    const auto& r_geometry = rElement.GetGeometry();
    const auto integration_method = rElement.GetIntegrationMethod();
    const auto& r_integration_points = r_geometry.IntegrationPoints(integration_method);
    const Matrix& r_N = r_geometry.ShapeFunctionsValues(integration_method);

    Matrix mass_matrix = ZeroMatrix(r_geometry.size(), r_geometry.size());
    for (IndexType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
        const auto N = row(r_N, i_point);
        noalias(mass_matrix) += r_integration_points[i_point].Weight() * outer_prod(N, N);
    }
    return mass_matrix;
}

/// Heat capacity \rho c from the density and specific heat of the settings
double ComputeHeatCapacity(const Element& rElement, const ConvectionDiffusionSettings& rSettings)
{
    const auto& r_properties = rElement.GetProperties();
    return r_properties[rSettings.GetDensityVariable()]
        * r_properties[rSettings.GetSpecificHeatVariable()];
}

} // namespace

// The volume fraction is constant over the element's knot span, so weighting
// the test functions (v_alpha = alpha v) scales the whole local system

void GsbmLaplacianElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    // (v, \rho c \partial_t u_h) with \partial_t u_h = \sum_k c_k u_h^{n+1-k}
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        const auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
        const double c_0 = rCurrentProcessInfo[BDF_COEFFICIENTS][0];
        noalias(rLeftHandSideMatrix) +=
            c_0 * ComputeHeatCapacity(*this, r_settings) * ComputeMassMatrix(*this);
    }

    rLeftHandSideMatrix *= GetValue(VOLUME_FRACTION);

    KRATOS_CATCH("")
}

void GsbmLaplacianElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // -\rho c M \sum_k c_k u_h^{n+1-k}, with u_h^{n+1} the current iterate
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        const auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
        const auto& r_unknown_variable = r_settings.GetUnknownVariable();
        const Vector& r_coefficients = rCurrentProcessInfo[BDF_COEFFICIENTS];
        const auto& r_geometry = GetGeometry();

        Vector time_derivative = ZeroVector(r_geometry.size());
        for (IndexType k = 0; k < r_coefficients.size(); ++k) {
            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                time_derivative[i] += r_coefficients[k]
                    * r_geometry[i].FastGetSolutionStepValue(r_unknown_variable, k);
            }
        }
        noalias(rRightHandSideVector) -= ComputeHeatCapacity(*this, r_settings)
            * prod(ComputeMassMatrix(*this), time_derivative);
    }

    rRightHandSideVector *= GetValue(VOLUME_FRACTION);

    KRATOS_CATCH("")
}

void GsbmLaplacianElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);
    rMassMatrix *= GetValue(VOLUME_FRACTION);

    KRATOS_CATCH("")
}

int GsbmLaplacianElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(Has(VOLUME_FRACTION)) << "GsbmLaplacianElement #" << Id()
        << ": VOLUME_FRACTION is not set. Run VolumeFractionGsbmProcess first."
        << std::endl;

    // The time term needs \rho c and as many steps as BDF coefficients
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        const auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
        KRATOS_ERROR_IF_NOT(
            r_settings.IsDefinedDensityVariable() && r_settings.IsDefinedSpecificHeatVariable()
        ) << "GsbmLaplacianElement #" << Id() << ": the transient problem needs the density "
            << "and specific heat variables in CONVECTION_DIFFUSION_SETTINGS." << std::endl;
        KRATOS_ERROR_IF_NOT(
            GetProperties().Has(r_settings.GetDensityVariable()) &&
            GetProperties().Has(r_settings.GetSpecificHeatVariable())
        ) << "GsbmLaplacianElement #" << Id() << ": no density or specific heat in the "
            << "properties." << std::endl;
        KRATOS_ERROR_IF(
            GetGeometry()[0].GetBufferSize() < rCurrentProcessInfo[BDF_COEFFICIENTS].size()
        ) << "GsbmLaplacianElement #" << Id() << ": the buffer size is smaller than the "
            << "number of BDF coefficients." << std::endl;
    }

    return BaseType::Check(rCurrentProcessInfo);
}

} // namespace Kratos
