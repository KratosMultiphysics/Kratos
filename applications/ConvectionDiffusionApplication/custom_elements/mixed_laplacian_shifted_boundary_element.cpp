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
#include "custom_elements/mixed_laplacian_shifted_boundary_element.h"


namespace Kratos
{

template<std::size_t TDim>
MixedLaplacianShiftedBoundaryElement<TDim>::MixedLaplacianShiftedBoundaryElement(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry)
    : MixedLaplacianElement<TDim, NumNodes>(NewId, pGeometry)
{
}

template<std::size_t TDim>
MixedLaplacianShiftedBoundaryElement<TDim>::MixedLaplacianShiftedBoundaryElement(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    typename PropertiesType::Pointer pProperties)
    : MixedLaplacianElement<TDim, NumNodes>(NewId, pGeometry, pProperties)
{
}

template<std::size_t TDim>
Element::Pointer MixedLaplacianShiftedBoundaryElement<TDim>::Create(
    IndexType NewId,
    typename BaseType::NodesArrayType const& ThisNodes,
    typename PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MixedLaplacianShiftedBoundaryElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TDim>
Element::Pointer MixedLaplacianShiftedBoundaryElement<TDim>::Create(
    IndexType NewId,
    typename GeometryType::Pointer pGeom,
    typename PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MixedLaplacianShiftedBoundaryElement>(NewId, pGeom, pProperties);
}

template<std::size_t TDim>
void MixedLaplacianShiftedBoundaryElement<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Add base Laplacian contribution
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    // Check if the element belongs to the surrogate interface
    // Note that the INTERFACE flag is assumed to be set in the layer of elements attached to the surrogate interface
    if (this->Is(INTERFACE)) {
        // Get convection-diffusion data container
        auto &r_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
        const auto& r_unknown_grad_var = r_settings.GetGradientVariable();
        const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();

        // Find the surrogate face local id
        // Note that it might happen that an interface element has no surrogate face (i.e. a unique node in the surrogate skin)
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();
        if (sur_bd_ids_vect.size() != 0) {
            // Get the parent geometry data
            double dom_size_parent;
            const auto& r_geom = this->GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, TDim> DN_DX_parent;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, dom_size_parent);
            const auto& r_boundaries = r_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_geom.NodesInFaces(nodes_in_faces);

            // Get the nodal values
            BoundedVector<double, NumNodes> nodal_diffusivity;
            BoundedVector<array_1d<double,TDim>, NumNodes> nodal_unknown_gradient;
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                nodal_diffusivity[i_node] = r_geom[i_node].FastGetSolutionStepValue(r_diffusivity_var);
                const auto& r_i_node_grad = r_geom[i_node].FastGetSolutionStepValue(r_unknown_grad_var);
                auto& r_nodal_unknown_gradient_i = nodal_unknown_gradient[i_node];
                for (std::size_t d = 0; d < TDim; ++d) {
                    r_nodal_unknown_gradient_i[d] = r_i_node_grad[d];
                }
            }

            // Loop the surrogate faces
            // Note that there is the chance that the surrogate face is not unique
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the current surrogate face geometry information
                const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
                const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
                const DenseVector<std::size_t> sur_bd_local_ids = row(nodes_in_faces, sur_bd_id);
                const auto& r_integration_points = r_sur_bd_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
                const std::size_t n_sur_bd_gauss = r_integration_points.size();
                Vector DetJ = ZeroVector(n_sur_bd_gauss);
                r_sur_bd_geom.DeterminantOfJacobian(DetJ, GeometryData::IntegrationMethod::GI_GAUSS_2);
                const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

                // Get the gradient of the node contrary to the surrogate face
                // Note that this is used to calculate the normal as n = - DN_DX_cont_node / norm_2(DN_DX_cont_node)
                BoundedVector<double,TDim> n_sur_bd = row(DN_DX_parent, sur_bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(n_sur_bd);
                n_sur_bd *= -h_sur_bd;

                // Allocate Gauss pt. auxiliary arrays
                double w_g;
                double aux;
                Vector aux_N(TDim);
                std::size_t i_loc_id;
                std::size_t j_loc_id;

                // Add the surrogate boundary flux contribution
                // Note that the local face ids. are already taken into account in the assembly
                for (std::size_t g = 0; g < n_sur_bd_gauss; ++g) {
                    // Calculate Gauss pt. geometry data
                    noalias(aux_N) = row(r_sur_bd_N, g);
                    w_g = DetJ[g] * r_integration_points[g].Weight();

                    // Add the nodal gradient flux contribution
                    for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
                        i_loc_id = sur_bd_local_ids[i_node + 1];
                        for (std::size_t j_node = 0; j_node < n_bd_points; ++j_node) {
                            j_loc_id = sur_bd_local_ids[j_node + 1];
                            aux = w_g * aux_N[i_node] * aux_N[j_node] * nodal_diffusivity[j_loc_id];
                            for (std::size_t d = 0; d < TDim; ++d) {
                                rLeftHandSideMatrix(i_loc_id*BlockSize, j_loc_id*BlockSize + d + 1) -= aux * n_sur_bd[d];
                            }
                            rRightHandSideVector(i_loc_id*BlockSize) += aux * inner_prod(nodal_unknown_gradient[j_loc_id], n_sur_bd);
                        }
                    }
                }
            }
        }
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim>
void MixedLaplacianShiftedBoundaryElement<TDim>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Add base Laplacian contribution
    BaseType::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    // Check if the element belongs to the surrogate interface
    // Note that the INTERFACE flag is assumed to be set in the layer of elements attached to the surrogate interface
    if (this->Is(INTERFACE)) {
        // Get convection-diffusion data container
        auto &r_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
        const auto& r_unknown_grad_var = r_settings.GetGradientVariable();
        const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();

        // Find the surrogate face local id
        // Note that it might happen that an interface element has no surrogate face (i.e. a unique node in the surrogate skin)
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();
        if (sur_bd_ids_vect.size() != 0) {
            // Get the parent geometry data
            double dom_size_parent;
            const auto& r_geom = this->GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, TDim> DN_DX_parent;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, dom_size_parent);
            const auto& r_boundaries = r_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_geom.NodesInFaces(nodes_in_faces);

            // Get the nodal values
            BoundedVector<double, NumNodes> nodal_diffusivity;
            BoundedVector<array_1d<double,TDim>, NumNodes> nodal_unknown_gradient;
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                nodal_diffusivity[i_node] = r_geom[i_node].FastGetSolutionStepValue(r_diffusivity_var);
                const auto& r_i_node_grad = r_geom[i_node].FastGetSolutionStepValue(r_unknown_grad_var);
                auto& r_nodal_unknown_gradient_i = nodal_unknown_gradient[i_node];
                for (std::size_t d = 0; d < TDim; ++d) {
                    r_nodal_unknown_gradient_i[d] = r_i_node_grad[d];
                }
            }

            // Loop the surrogate faces
            // Note that there is the chance that the surrogate face is not unique
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the current surrogate face geometry information
                const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
                const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
                const DenseVector<std::size_t> sur_bd_local_ids = row(nodes_in_faces, sur_bd_id);
                const auto& r_integration_points = r_sur_bd_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
                const std::size_t n_sur_bd_gauss = r_integration_points.size();
                Vector DetJ = ZeroVector(n_sur_bd_gauss);
                r_sur_bd_geom.DeterminantOfJacobian(DetJ, GeometryData::IntegrationMethod::GI_GAUSS_2);
                const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

                // Get the gradient of the node contrary to the surrogate face
                // Note that this is used to calculate the normal as n = - DN_DX_cont_node / norm_2(DN_DX_cont_node)
                BoundedVector<double,TDim> n_sur_bd = row(DN_DX_parent, sur_bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(n_sur_bd);
                n_sur_bd *= -h_sur_bd;

                // Allocate Gauss pt. auxiliary arrays
                double w_g;
                double aux;
                Vector aux_N(TDim);
                std::size_t i_loc_id;
                std::size_t j_loc_id;

                // Add the surrogate boundary flux contribution
                // Note that the local face ids. are already taken into account in the assembly
                for (std::size_t g = 0; g < n_sur_bd_gauss; ++g) {
                    // Calculate Gauss pt. geometry data
                    noalias(aux_N) = row(r_sur_bd_N, g);
                    w_g = DetJ[g] * r_integration_points[g].Weight();

                    // Add the nodal gradient flux contribution
                    for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
                        i_loc_id = sur_bd_local_ids[i_node + 1];
                        for (std::size_t j_node = 0; j_node < n_bd_points; ++j_node) {
                            j_loc_id = sur_bd_local_ids[j_node + 1];
                            aux = w_g * aux_N[i_node] * aux_N[j_node] * nodal_diffusivity[j_loc_id];
                            for (std::size_t d = 0; d < TDim; ++d) {
                                rLeftHandSideMatrix(i_loc_id*BlockSize, j_loc_id*BlockSize + d + 1) -= aux * n_sur_bd[d];
                            }
                        }
                    }
                }
            }
        }
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim>
void MixedLaplacianShiftedBoundaryElement<TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
KRATOS_TRY

    // Add base Laplacian contribution
    BaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Check if the element belongs to the surrogate interface
    // Note that the INTERFACE flag is assumed to be set in the layer of elements attached to the surrogate interface
    if (this->Is(INTERFACE)) {
        // Get convection-diffusion data container
        auto &r_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
        const auto& r_unknown_grad_var = r_settings.GetGradientVariable();
        const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();

        // Find the surrogate face local id
        // Note that it might happen that an interface element has no surrogate face (i.e. a unique node in the surrogate skin)
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();
        if (sur_bd_ids_vect.size() != 0) {
            // Get the parent geometry data
            double dom_size_parent;
            const auto& r_geom = this->GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, TDim> DN_DX_parent;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, dom_size_parent);
            const auto& r_boundaries = r_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_geom.NodesInFaces(nodes_in_faces);

            // Get the nodal values
            BoundedVector<double, NumNodes> nodal_diffusivity;
            BoundedVector<array_1d<double,TDim>, NumNodes> nodal_unknown_gradient;
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                nodal_diffusivity[i_node] = r_geom[i_node].FastGetSolutionStepValue(r_diffusivity_var);
                const auto& r_i_node_grad = r_geom[i_node].FastGetSolutionStepValue(r_unknown_grad_var);
                auto& r_nodal_unknown_gradient_i = nodal_unknown_gradient[i_node];
                for (std::size_t d = 0; d < TDim; ++d) {
                    r_nodal_unknown_gradient_i[d] = r_i_node_grad[d];
                }
            }

            // Loop the surrogate faces
            // Note that there is the chance that the surrogate face is not unique
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the current surrogate face geometry information
                const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
                const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
                const DenseVector<std::size_t> sur_bd_local_ids = row(nodes_in_faces, sur_bd_id);
                const auto& r_integration_points = r_sur_bd_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
                const std::size_t n_sur_bd_gauss = r_integration_points.size();
                Vector DetJ = ZeroVector(n_sur_bd_gauss);
                r_sur_bd_geom.DeterminantOfJacobian(DetJ, GeometryData::IntegrationMethod::GI_GAUSS_2);
                const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

                // Get the gradient of the node contrary to the surrogate face
                // Note that this is used to calculate the normal as n = - DN_DX_cont_node / norm_2(DN_DX_cont_node)
                BoundedVector<double,TDim> n_sur_bd = row(DN_DX_parent, sur_bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(n_sur_bd);
                n_sur_bd *= -h_sur_bd;

                // Allocate Gauss pt. auxiliary arrays
                double w_g;
                double aux;
                Vector aux_N(TDim);
                std::size_t i_loc_id;
                std::size_t j_loc_id;

                // Add the surrogate boundary flux contribution
                // Note that the local face ids. are already taken into account in the assembly
                for (std::size_t g = 0; g < n_sur_bd_gauss; ++g) {
                    // Calculate Gauss pt. geometry data
                    noalias(aux_N) = row(r_sur_bd_N, g);
                    w_g = DetJ[g] * r_integration_points[g].Weight();

                    // Add the nodal gradient flux contribution
                    for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
                        i_loc_id = sur_bd_local_ids[i_node + 1];
                        for (std::size_t j_node = 0; j_node < n_bd_points; ++j_node) {
                            j_loc_id = sur_bd_local_ids[j_node + 1];
                            aux = w_g * aux_N[i_node] * aux_N[j_node] * nodal_diffusivity[j_loc_id];
                            rRightHandSideVector(i_loc_id*BlockSize) += aux * inner_prod(nodal_unknown_gradient[j_loc_id], n_sur_bd);
                        }
                    }
                }
            }
        }
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim>
int MixedLaplacianShiftedBoundaryElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Base Laplacian element check
    return BaseType::Check(rCurrentProcessInfo);
}

template<std::size_t TDim>
std::vector<std::size_t> MixedLaplacianShiftedBoundaryElement<TDim>::GetSurrogateFacesIds()
{
    const std::size_t n_faces = TDim + 1;
    auto& r_neigh_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

    // Check the current element faces
    // Note that we relly on the fact that the neighbours are sorted according to the faces
    std::vector<std::size_t> surrogate_faces_ids;
    for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
        auto p_neigh_elem = r_neigh_elems(i_face).get();
        if (p_neigh_elem != nullptr && p_neigh_elem->Is(BOUNDARY)) {
            surrogate_faces_ids.push_back(i_face);
        }
    }

    return surrogate_faces_ids;
}

template class MixedLaplacianShiftedBoundaryElement<2>;
template class MixedLaplacianShiftedBoundaryElement<3>;

} // Namespace Kratos
