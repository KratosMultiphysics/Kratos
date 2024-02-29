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
    : LaplacianElement(NewId, pGeometry)
{
}

template<std::size_t TDim>
LaplacianShiftedBoundaryElement<TDim>::LaplacianShiftedBoundaryElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : LaplacianElement(NewId, pGeometry, pProperties)
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
        // Note that it might happen that an interface element has no surrogate face (i.e. a unique node in the surrogate skin)
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();
        if (sur_bd_ids_vect.size() != 0) {
            // Get the parent geometry data
            double dom_size_parent;
            const auto& r_geom = GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, TDim> DN_DX_parent;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, dom_size_parent);
            const auto& r_boundaries = r_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_geom.NodesInFaces(nodes_in_faces);

            // Get the unknowns vector
            BoundedVector<double, NumNodes> nodal_unknown;
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                nodal_unknown[i_node] = r_geom[i_node].FastGetSolutionStepValue(r_unknown_var);
            }

            // Loop the surrogate faces
            // Note that there is the chance that the surrogate face is not unique
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the current surrogate face geometry information
                const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
                const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
                const DenseVector<std::size_t> sur_bd_local_ids = row(nodes_in_faces, sur_bd_id);
                const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

                // Get the surrogate boundary average conductivity
                double k_avg = 0.0;
                for (std::size_t i_bd_node = 0; i_bd_node < n_bd_points; ++i_bd_node) {
                    k_avg += r_sur_bd_geom[i_bd_node].FastGetSolutionStepValue(r_diffusivity_var);
                }
                k_avg /= n_bd_points;

                // Get the gradient of the node contrary to the surrogate face
                // Note that this is used to calculate the normal as n = - DN_DX_cont_node / norm_2(DN_DX_cont_node)
                // const BoundedVector<double,TDim> DN_DX_cont_node = row(DN_DX_parent, sur_bd_local_ids[0]);
                BoundedVector<double,TDim> n_sur_bd = row(DN_DX_parent, sur_bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(n_sur_bd);
                n_sur_bd *= -h_sur_bd;

                // Calculate the gradient projection
                const BoundedVector<double,NumNodes> DN_DX_proj_n = prod(DN_DX_parent, n_sur_bd);

                // Add the surrogate boundary flux contribution
                // Note that the local face ids. are already taken into account in the assembly
                // Note that the integration weight is calculated as TDim * Parent domain size * norm(DN_DX_cont_node)
                double aux_1;
                double aux_2;
                std::size_t i_loc_id;
                BoundedVector<double,TDim> j_node_grad;
                const double aux_w_k = TDim * dom_size_parent * k_avg / h_sur_bd;
                for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
                    aux_1 = aux_w_k * r_sur_bd_N(0,i_node);
                    i_loc_id = sur_bd_local_ids[i_node + 1];
                    for (std::size_t j_node = 0; j_node < NumNodes; ++j_node) {
                        aux_2 = aux_1 * DN_DX_proj_n(j_node);
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

    // Check if the element belongs to the surrogate interface
    // Note that the INTERFACE flag is assumed to be set in the layer of elements attached to the surrogate interface
    if (Is(INTERFACE)) {
        // Get convection-diffusion data container
        auto p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
        auto &r_settings = *p_settings;
        const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();

        // Find the surrogate face local id
        // Note that it might happen that an interface element has no surrogate face (i.e. a unique node in the surrogate skin)
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();
        if (sur_bd_ids_vect.size() != 0) {
            // Get the parent geometry data
            double dom_size_parent;
            const auto& r_geom = GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, TDim> DN_DX_parent;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, dom_size_parent);
            const auto& r_boundaries = r_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_geom.NodesInFaces(nodes_in_faces);

            // Loop the surrogate faces
            // Note that there is the chance that the surrogate face is not unique
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the current surrogate face geometry information
                const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
                const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
                const DenseVector<std::size_t> sur_bd_local_ids = row(nodes_in_faces, sur_bd_id);
                const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

                // Get the surrogate boundary average conductivity
                double k_avg = 0.0;
                for (std::size_t i_bd_node = 0; i_bd_node < n_bd_points; ++i_bd_node) {
                    k_avg += r_sur_bd_geom[i_bd_node].FastGetSolutionStepValue(r_diffusivity_var);
                }
                k_avg /= n_bd_points;

                // Get the gradient of the node contrary to the surrogate face
                // Note that this is used to calculate the normal as n = - DN_DX_cont_node / norm_2(DN_DX_cont_node)
                // const BoundedVector<double,TDim> DN_DX_cont_node = row(DN_DX_parent, sur_bd_local_ids[0]);
                BoundedVector<double,TDim> n_sur_bd = row(DN_DX_parent, sur_bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(n_sur_bd);
                n_sur_bd *= -h_sur_bd;

                // Calculate the gradient projection
                const BoundedVector<double,NumNodes> DN_DX_proj_n = prod(DN_DX_parent, n_sur_bd);

                // Add the surrogate boundary flux contribution
                // Note that the local face ids. are already taken into account in the assembly
                // Note that the integration weight is calculated as TDim * Parent domain size * norm(DN_DX_cont_node)
                double aux_1;
                double aux_2;
                std::size_t i_loc_id;
                BoundedVector<double,TDim> j_node_grad;
                const double aux_w_k = TDim * dom_size_parent * k_avg / h_sur_bd;
                for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
                    aux_1 = aux_w_k * r_sur_bd_N(0,i_node);
                    i_loc_id = sur_bd_local_ids[i_node + 1];
                    for (std::size_t j_node = 0; j_node < NumNodes; ++j_node) {
                        aux_2 = aux_1 * DN_DX_proj_n(j_node);
                        rLeftHandSideMatrix(i_loc_id, j_node) -= aux_2;
                    }
                }
            }
        }
    }
}

template<std::size_t TDim>
void LaplacianShiftedBoundaryElement<TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Add base Laplacian contribution
    BaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Check if the element belongs to the surrogate interface
    // Note that the INTERFACE flag is assumed to be set in the layer of elements attached to the surrogate interface
    if (Is(INTERFACE)) {
        // Get convection-diffusion data container
        auto p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
        auto &r_settings = *p_settings;
        const auto& r_unknown_var = r_settings.GetUnknownVariable();
        const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();

        // Find the surrogate face local id
        // Note that it might happen that an interface element has no surrogate face (i.e. a unique node in the surrogate skin)
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();
        if (sur_bd_ids_vect.size() != 0) {
            // Get the parent geometry data
            double dom_size_parent;
            const auto& r_geom = GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, TDim> DN_DX_parent;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, dom_size_parent);
            const auto& r_boundaries = r_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_geom.NodesInFaces(nodes_in_faces);

            // Get the unknowns vector
            BoundedVector<double, NumNodes> nodal_unknown;
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                nodal_unknown[i_node] = r_geom[i_node].FastGetSolutionStepValue(r_unknown_var);
            }

            // Loop the surrogate faces
            // Note that there is the chance that the surrogate face is not unique
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the current surrogate face geometry information
                const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
                const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
                const DenseVector<std::size_t> sur_bd_local_ids = row(nodes_in_faces, sur_bd_id);
                const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

                // Get the surrogate boundary average conductivity
                double k_avg = 0.0;
                for (std::size_t i_bd_node = 0; i_bd_node < n_bd_points; ++i_bd_node) {
                    k_avg += r_sur_bd_geom[i_bd_node].FastGetSolutionStepValue(r_diffusivity_var);
                }
                k_avg /= n_bd_points;

                // Get the gradient of the node contrary to the surrogate face
                // Note that this is used to calculate the normal as n = - DN_DX_cont_node / norm_2(DN_DX_cont_node)
                // const BoundedVector<double,TDim> DN_DX_cont_node = row(DN_DX_parent, sur_bd_local_ids[0]);
                BoundedVector<double,TDim> n_sur_bd = row(DN_DX_parent, sur_bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(n_sur_bd);
                n_sur_bd *= -h_sur_bd;

                // Calculate the gradient projection
                const BoundedVector<double,NumNodes> DN_DX_proj_n = prod(DN_DX_parent, n_sur_bd);

                // Add the surrogate boundary flux contribution
                // Note that the local face ids. are already taken into account in the assembly
                // Note that the integration weight is calculated as TDim * Parent domain size * norm(DN_DX_cont_node)
                double aux_1;
                double aux_2;
                std::size_t i_loc_id;
                BoundedVector<double,TDim> j_node_grad;
                const double aux_w_k = TDim * dom_size_parent * k_avg / h_sur_bd;
                for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
                    aux_1 = aux_w_k * r_sur_bd_N(0,i_node);
                    i_loc_id = sur_bd_local_ids[i_node + 1];
                    for (std::size_t j_node = 0; j_node < NumNodes; ++j_node) {
                        aux_2 = aux_1 * DN_DX_proj_n(j_node);
                        rRightHandSideVector(i_loc_id) += aux_2 * nodal_unknown(j_node);
                    }
                }
            }
        }
    }
}

template<std::size_t TDim>
int LaplacianShiftedBoundaryElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Base Laplacian element check
    return BaseType::Check(rCurrentProcessInfo);
}

template<std::size_t TDim>
std::vector<std::size_t> LaplacianShiftedBoundaryElement<TDim>::GetSurrogateFacesIds()
{
    const std::size_t n_faces = TDim + 1;
    auto& r_neigh_elems = GetValue(NEIGHBOUR_ELEMENTS);

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

template class LaplacianShiftedBoundaryElement<2>;
template class LaplacianShiftedBoundaryElement<3>;

} // Namespace Kratos
