// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "custom_elements/laplacian_shifted_boundary_element.h"


namespace Kratos
{

LaplacianShiftedBoundaryElement::LaplacianShiftedBoundaryElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : LaplacianElement(
        NewId,
        pGeometry)
{
}

LaplacianShiftedBoundaryElement::LaplacianShiftedBoundaryElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : LaplacianElement(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer LaplacianShiftedBoundaryElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundaryElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer LaplacianShiftedBoundaryElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundaryElement>(NewId, pGeom, pProperties);
}

LaplacianShiftedBoundaryElement::~LaplacianShiftedBoundaryElement()
{
}

void LaplacianShiftedBoundaryElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Add base Laplacian contribution
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    // Get convection-diffusion data container
    auto p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto &r_settings = *p_settings;
    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();

    // Find the surrogate face local id
    std::size_t sur_bd_id;
    const auto& r_geom = GetGeometry();
    const auto r_boundaries = r_geom.GenerateBoundariesEntities();
    const std::size_t n_bds = r_boundaries.size();
    const std::size_t n_bd_nodes = r_boundaries[0].PointsNumber();
    for (std::size_t i_bd = 0; i_bd < n_bds; ++i_bd) {
        const auto& r_boundary = r_boundaries[i_bd];
        std::size_t n_flag_nodes = 0;
        for (const auto& r_node : r_boundary) {
            if (r_node.Is(BOUNDARY)) {
                n_flag_nodes++;
            }
        }
        if (n_flag_nodes == n_bd_nodes) {
            sur_bd_id = i_bd;
            break;
        }
    }

    // Get the surrogate boundary nodal data
    DenseVector<double> nodal_unknown(n_bd_points);
    DenseVector<double> nodal_conductivity(n_bd_points);
    for (std::size_t i_bd_node = 0; i_bd_node < n_bd_points; i_bd_node) {
        nodal_unknown[i_bd_node] = r_sur_bd_geom[i_bd_node].FastGetSolutionStepValue(r_unknown_var);
        nodal_conductivity[i_bd_node] = r_sur_bd_geom[i_bd_node].FastGetSolutionStepValue(r_diffusivity_var);
    }

    // Get the surrogate boundary geometry data
    const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
    const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
    const unsigned int dim = r_sur_bd_geom.WorkingSpaceDimension();

    const auto& r_integration_points = r_sur_bd_geom.IntegrationPoints(BaseType::GetIntegrationMethod());
    const auto& DN_De = r_sur_bd_geom.ShapeFunctionsLocalGradients(BaseType::GetIntegrationMethod());
    const Matrix& N_g = r_sur_bd_geom.ShapeFunctionsValues(BaseType::GetIntegrationMethod());
    Element::GeometryType::JacobiansType J0;
    r_sur_bd_geom.Jacobian(J0, BaseType::GetIntegrationMethod());

    // Integrate the surrogate boundary flux
    double det_J0;
    Matrix DN_DX(n_bd_points, dim);
    Matrix inv_J0(dim, dim);
    DenseVector temp(n_bd_points);
    DenseVector grad_proj(n_bd_points);
    Matrix aux_LHS = ZeroMatrix(n_bd_points);
    DenseVector aux_RHS = ZeroVector(n_bd_points);
    const std::size_t n_gauss = r_integration_points.size();
    for(std::size_t i_point = 0; i_point < n_gauss; ++i_point) {
        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[i_point], inv_J0, det_J0);

        // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[i_point], inv_J0);

        const auto N_g = row(N_g, i_point);
        const double w_g = r_integration_points[i_point].Weight() * det_J0;
        const double k_g = inner_prod(N_g, nodal_conductivity);
        const array_1d<double,3> normal_g = r_sur_bd_geom.Normal(r_integration_points[i_point].Coordinates());

        //TODO: THIS CAN BE INCLUDED IN THE ONE BELOW
        grad_proj = ZeroVector(n_bd_points);
        for (std::size_t d = 0; d < dim; ++d) {
            for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
                grad_proj(d) += DN_DX(i_node, dim) * normal_g(dim);
            }
        }

        for (std::size_t i_node = 0; i_node = n_bd_points; ++i_node) {
            for (std::size_t j_node = 0; j_node = n_bd_points; ++j_node) {
                aux_LHS(i_node, j_node) += w_g * k_g * N_g(i_node) * grad_proj(j_node);
            }
        }
        // noalias(rLeftHandSideMatrix) += w_g * k_g * prod(DN_DX, trans(DN_DX)); //

        // // Calculating the local RHS
        // const double qgauss = inner_prod(N, heat_flux_local);

        // noalias(rRightHandSideVector) += IntToReferenceWeight * qgauss * N;
    }

    // Do the assembly with the NodesInFace

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Add base Laplacian contribution
    BaseType::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    // Add the surrogate boundary face contribution

}

void LaplacianShiftedBoundaryElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Add base Laplacian contribution
    BaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Add the surrogate boundary face contribution

}

int LaplacianShiftedBoundaryElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Surrogate boundary checks

    // Base Laplacian element check
    return BaseType::Check(rCurrentProcessInfo);
}

} // Namespace Kratos


