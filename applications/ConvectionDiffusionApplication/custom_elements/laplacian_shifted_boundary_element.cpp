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
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/laplacian_shifted_boundary_element.h"


namespace Kratos
{

template<std::size_t TDim>
LaplacianShiftedBoundaryElement<TDim>::LaplacianShiftedBoundaryElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : LaplacianElement(
        NewId,
        pGeometry)
{
}

template<std::size_t TDim>
LaplacianShiftedBoundaryElement<TDim>::LaplacianShiftedBoundaryElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : LaplacianElement(
        NewId,
        pGeometry,
        pProperties)
{
}

template<std::size_t TDim>
Element::Pointer LaplacianShiftedBoundaryElement<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundaryElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TDim>
Element::Pointer LaplacianShiftedBoundaryElement<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundaryElement>(NewId, pGeom, pProperties);
}

template<std::size_t TDim>
LaplacianShiftedBoundaryElement<TDim>::~LaplacianShiftedBoundaryElement()
{
}

template<std::size_t TDim>
void LaplacianShiftedBoundaryElement<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Add base Laplacian contribution
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    // Check if the element belongs to the surrogate interface
    // Note that the INTERFACE flag is assumed to be set in the layer of elements attached to the surrogate interface
    if (Is(INTERFACE)) {
        // Get convection-diffusion data container
        auto p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
        auto &r_settings = *p_settings;
        const auto& r_unknown_var = r_settings.GetUnknownVariable();
        const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();

        // Find the surrogate face local id
        // Note that the BOUNDARY flag is assumed to be set in the surrogate interface nodes
        std::size_t sur_bd_id;
        const auto& r_geom = GetGeometry();
        const std::size_t n_nodes = r_geom.PointsNumber();
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

        // Get the parent geometry data
        double dom_size_parent;
        array_1d<double, TNumNodes> N_parent;
        BoundedMatrix<double, TNumNodes, TDim> DN_DX_parent;
        GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, dom_size_parent);

        // Get the surrogate boundary geometry data
        const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
        const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
        const unsigned int n_dim = r_sur_bd_geom.WorkingSpaceDimension();

        const auto& r_integration_points = r_sur_bd_geom.IntegrationPoints(BaseType::GetIntegrationMethod());
        const Matrix& N_g = r_sur_bd_geom.ShapeFunctionsValues(BaseType::GetIntegrationMethod());

        // Get the surrogate boundary NodesInFaces in order to perform the assembly
        // Note that the first node in the NodesInFace is the node contrary to the face
        DenseMatrix<unsigned int> nodes_in_faces;
        r_geom.NodesInFaces(nodes_in_faces);
        const DenseVector<std::size_t> sur_bd_local_ids = row(nodes_in_faces, sur_bd_id);

        // Get the surrogate boundary nodal data
        DenseVector<double> nodal_unknown(n_bd_points);
        DenseVector<double> nodal_conductivity(n_bd_points);
        for (std::size_t i_bd_node = 0; i_bd_node < n_bd_points; ++i_bd_node) {
            nodal_unknown[i_bd_node] = r_sur_bd_geom[i_bd_node].FastGetSolutionStepValue(r_unknown_var);
            nodal_conductivity[i_bd_node] = r_sur_bd_geom[i_bd_node].FastGetSolutionStepValue(r_diffusivity_var);
        }

        // Integrate the surrogate boundary flux
        const std::size_t n_gauss = r_integration_points.size();
        for(std::size_t i_g = 0; i_g < n_gauss; ++i_g) {
            // Get surrogate boundary Gauss pt. geometry data
            const DenseVector<double> N_i_g = row(N_g, i_g);
            const double w_i_g = r_integration_points[i_g].Weight() * r_sur_bd_geom.DeterminantOfJacobian(i_g, BaseType::GetIntegrationMethod());
            const double k_i_g = inner_prod(N_i_g, nodal_conductivity);
            const array_1d<double,3> normal_g = r_sur_bd_geom.Normal(r_integration_points[i_g].Coordinates());

            // Add the surrogate boundary flux contribution
            // Note that the local face ids. are already taken into account in the assembly
            double aux_1;
            double aux_2;
            std::size_t i_loc_id;
            for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
                i_loc_id = sur_bd_local_ids[i_node + 1];
                aux_1 = w_i_g * k_i_g * N_i_g(i_node);
                for (std::size_t j_node = 0; j_node < n_bd_points; ++j_node) {
                    for (std::size_t d = 0; d < n_dim; ++d) {
                        aux_2 = aux_1 * DN_DX_parent(j_node, d) * normal_g(d);
                        rLeftHandSideMatrix(i_loc_id, j_node) -= aux_2;
                        rRightHandSideVector(i_loc_id) += aux_2 * nodal_unknown(j_node);
                    }
                }
            }
        }
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim>
void LaplacianShiftedBoundaryElement<TDim>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Add base Laplacian contribution
    BaseType::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    // Add the surrogate boundary face contribution
    //TODO:
}

template<std::size_t TDim>
void LaplacianShiftedBoundaryElement<TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Add base Laplacian contribution
    BaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Add the surrogate boundary face contribution
    //TODO:
}

template<std::size_t TDim>
int LaplacianShiftedBoundaryElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Surrogate boundary checks

    // Base Laplacian element check
    return BaseType::Check(rCurrentProcessInfo);
}

template class LaplacianShiftedBoundaryElement<2>;
template class LaplacianShiftedBoundaryElement<3>;

} // Namespace Kratos
