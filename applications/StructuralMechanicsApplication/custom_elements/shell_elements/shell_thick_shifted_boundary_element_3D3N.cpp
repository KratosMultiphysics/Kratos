// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
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
#include "utilities/geometry_utilities.h"

// Application includes
#include "shell_thick_shifted_boundary_element_3D3N.hpp"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{

template <ShellKinematics TKinematics>
ShellThickShiftedBoundaryElement3D3N<TKinematics>::ShellThickShiftedBoundaryElement3D3N(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

template <ShellKinematics TKinematics>
ShellThickShiftedBoundaryElement3D3N<TKinematics>::ShellThickShiftedBoundaryElement3D3N(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

template <ShellKinematics TKinematics>
void ShellThickShiftedBoundaryElement3D3N<TKinematics>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // CutFEM-style formulation
    const int sbm_formulation_type_dispatch = this->GetProperties().Has(SBM_FORMULATION_TYPE)
        ? this->GetProperties()[SBM_FORMULATION_TYPE]
        : 0;
    if (sbm_formulation_type_dispatch == 4) {
        CalculateLocalSystemCutFEM(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        return;
    }

    // Add base Laplacian contribution
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    // Check if the element belongs to the surrogate interface
    // Note that the INTERFACE flag is assumed to be set in the layer of elements attached to the surrogate interface
    if (this->Is(INTERFACE)) {
        // Find the surrogate face local id
        // Note that it might happen that an interface element has no surrogate face (i.e. a unique node in the surrogate skin)
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();

        Matrix left_hand_side = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
        Vector right_hand_side = ZeroVector(rRightHandSideVector.size());

        Matrix left_hand_side_taylor = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
        Matrix left_hand_side_taylor_penalty = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
        Matrix left_hand_side_taylor_gap = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
        Matrix left_hand_side_taylor_penalty_gap = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());

        Matrix left_hand_side_gap_SBM = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());

        // Gap SBM Neumann-consistency term
        Matrix left_hand_side_gap_neumann = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
        Vector right_hand_side_shear_neumann = ZeroVector(rRightHandSideVector.size());

        // Per-DOF constrained mask (ux,uy,uz,rx,ry,rz)
        array_1d<double, 6> constrained_dofs;
        if (this->GetProperties().Has(SBM_CONSTRAINED_DOFS)) {
            const auto& r_constrained_dofs = this->GetProperties()[SBM_CONSTRAINED_DOFS];
            KRATOS_ERROR_IF(r_constrained_dofs.size() != 6)
                << "SBM_CONSTRAINED_DOFS must have size 6 (ux,uy,uz,rx,ry,rz) in Properties " << this->GetProperties().Id() << std::endl;
            for (std::size_t d = 0; d < 6; ++d) {
                constrained_dofs[d] = r_constrained_dofs[d];
            }
        } else {
            for (std::size_t d = 0; d < 6; ++d) {
                constrained_dofs[d] = 1.0;
            }
        }

        // Surrogate-boundary flux formulation choice: 0 = MLS (default), 1 = Taylor
        const int sbm_formulation_type = this->GetProperties().Has(SBM_FORMULATION_TYPE)
            ? this->GetProperties()[SBM_FORMULATION_TYPE]
            : 0;

        // Current nodal state, needed for the Taylor/Gap SBM penalty residual correction below.
        Vector unknown_values = ZeroVector(LocalSize);
        {
            const auto& r_geom_outer = this->GetGeometry();
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                const auto& r_disp = r_geom_outer[i_node].FastGetSolutionStepValue(DISPLACEMENT);
                const auto& r_rot = r_geom_outer[i_node].FastGetSolutionStepValue(ROTATION);
                for (std::size_t d = 0; d < 3; ++d) {
                    unknown_values(i_node*BlockSize + d) = r_disp[d];
                    unknown_values(i_node*BlockSize + d + 3) = r_rot[d];
                }
            }
        }

        if (sur_bd_ids_vect.size() != 0) {
            typename BaseType::CalculationData data(this->mpCoordinateTransformation, rCurrentProcessInfo);
            data.CalculateLHS = true;
            data.CalculateRHS = true;

            BaseType::InitializeCalculationData(data);

            array_1d<double,3>& gp0 = data.gpLocations[0];
            gp0[0] = 1.0/3.0;
            gp0[1] = 1.0/3.0;
            gp0[2] = 1.0/3.0;

            data.gpIndex = 0;

            noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

            BaseType::CalculateSectionResponse(data);

            // Get the parent geometry data
            const double dom_size_parent = data.TotalArea;
            const double lcs_x12 = data.LCS0.X1() - data.LCS0.X2();
            const double lcs_x23 = data.LCS0.X2() - data.LCS0.X3();
            const double lcs_x31 = data.LCS0.X3() - data.LCS0.X1();
            const double lcs_y12 = data.LCS0.Y1() - data.LCS0.Y2();
            const double lcs_y23 = data.LCS0.Y2() - data.LCS0.Y3();
            const double lcs_y31 = data.LCS0.Y3() - data.LCS0.Y1();
            const double lcs_A2 = 2.0 * data.TotalArea;

            const auto& r_geom = this->GetGeometry();

            array_1d<double, NumNodes> N_parent;
            N_parent[0] = N_parent[1] = N_parent[2] = 1.0 / 3.0;
            BoundedMatrix<double, NumNodes, 2> DN_DX_parent;
            DN_DX_parent(0,0) = -lcs_y23 / lcs_A2;  DN_DX_parent(0,1) =  lcs_x23 / lcs_A2;
            DN_DX_parent(1,0) = -lcs_y31 / lcs_A2;  DN_DX_parent(1,1) =  lcs_x31 / lcs_A2;
            DN_DX_parent(2,0) = -lcs_y12 / lcs_A2;  DN_DX_parent(2,1) =  lcs_x12 / lcs_A2;
            BoundedMatrix<double,8,LocalSize> B;
            const auto &r_boundaries = r_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_geom.NodesInFaces(nodes_in_faces);

            // Prescribed Dirichlet target (u0, theta0)
            const auto& r_penalty_bc_val = this->GetValue(DISPLACEMENT);
            const auto& r_penalty_bc_val_rot = this->GetValue(ROTATION);
            array_1d<double,6> penalty_bc_values;
            for (std::size_t d = 0; d < 3; ++d) {
                penalty_bc_values[d] = r_penalty_bc_val[d];
                penalty_bc_values[d + 3] = r_penalty_bc_val_rot[d];
            }

            // Loop the surrogate faces
            // Note that there is the chance that the surrogate face is not unique
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                    // Get the current surrogate face geometry information
                    const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
                    const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
                    const DenseVector<std::size_t> sur_bd_local_ids = column(nodes_in_faces, sur_bd_id);
                    const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

                    // Get the gradient of the node contrary to the surrogate face
                    // Note that this is used to calculate the normal as n = - DN_DX_cont_node / norm_2(DN_DX_cont_node)
                    BoundedVector<double,2> n_sur_bd = row(DN_DX_parent, sur_bd_local_ids[0]);

                    const double h_sur_bd = 1.0 / norm_2(n_sur_bd);
                    n_sur_bd *= -h_sur_bd;

                    // Calculate the required projections
                    array_1d<double,6> cauchy_traction;
                    BoundedMatrix<double,6,LocalSize> aux_CB_projection;
                    aux_CB_projection = ZeroMatrix(6, LocalSize);

                    const auto &r_stress = data.generalizedStresses; //cons_law_values.GetStressVector();
                    const auto &r_C = data.D; //cons_law_values.GetConstitutiveMatrix();

                    // TODO!! Modify these functions to work in 3D
                    B = ZeroMatrix(8, 18);

                    for ( IndexType i = 0; i < 3; ++i ) {
                        const IndexType initial_index = i*6;
                        B(0, initial_index    ) = DN_DX_parent(i, 0);
                        B(1, initial_index + 1) = DN_DX_parent(i, 1);
                        B(2, initial_index    ) = DN_DX_parent(i, 1);
                        B(2, initial_index + 1) = DN_DX_parent(i, 0);

                        B(3, initial_index + 4) = DN_DX_parent(i, 0);
                        B(4, initial_index + 3) = -DN_DX_parent(i, 1);
                        B(5, initial_index + 3) = -DN_DX_parent(i, 0);
                        B(5, initial_index + 4) = DN_DX_parent(i, 1);

                        B(6, initial_index + 2) = DN_DX_parent(i, 0);
                        B(6, initial_index + 4) = N_parent[i];
                        B(7, initial_index + 2) = DN_DX_parent(i, 1);
                        B(7, initial_index + 3) = -N_parent[i];
                    }

                    Matrix D(8, 8);
                    noalias(D) = ZeroMatrix(8, 8);

                    for(std::size_t i = 0; i < r_C.size1(); ++i)
                    {
                        for(std::size_t j = 0; j < r_C.size1(); ++j)
                        {
                            D(i,j) = r_C(i,j);
                        }
                    }

                    BoundedMatrix<double,8,LocalSize> B_temp = ZeroMatrix(8, 18);
                    BoundedMatrix<double,8,LocalSize> B_parent = ZeroMatrix(8, 18);

                    MatrixType T(18, 18);
                    data.LCS.ComputeTotalRotationMatrix(T);
                    B_temp = prod(data.B, T); 

                    Matrix R(8, 8);
                    this->mSections[0]->GetRotationMatrixForGeneralizedStrains(-(this->mSections[0]->GetOrientationAngle()), R);

                    R(6,7) = -1.0*R(6,7); 
                    R(7,6) = -1.0*R(7,6);

                    B_parent = prod(R, B_temp);

                    const array_1d<double,3> local_e1 = -data.LCS.Vx();
                    const array_1d<double,3> local_e2 = -data.LCS.Vy();
                    const array_1d<double,3> local_e3 = data.LCS.Vz();
                    CalculateCBProjectionLinearisation(D, B_parent, n_sur_bd, local_e1, local_e2, local_e3, aux_CB_projection);
                    CalculateCauchyTractionVector(r_stress, n_sur_bd, local_e1, local_e2, local_e3, cauchy_traction);

                    // Add the surrogate boundary flux contribution
                    // Note that the local face ids. are already taken into account in the assembly
                    // Note that the integration weight is calculated as 2 * Parent domain size * norm(DN_DX_cont_node)
                    double aux_val;
                    std::size_t i_loc_id;
                    const double aux_w = 2 * dom_size_parent / h_sur_bd;

                    for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
                        aux_val = aux_w * r_sur_bd_N(0,i_node);
                        i_loc_id = sur_bd_local_ids[i_node + 1];
                        for (std::size_t d = 0; d < 6; ++d) {
                            if (constrained_dofs[d] == 0.0) {
                                continue;
                            }
                            right_hand_side(i_loc_id*BlockSize+d) += aux_val * cauchy_traction[d];
                            for (std::size_t j_node = 0; j_node < NumNodes * 6; ++j_node) {
                                left_hand_side(i_loc_id*BlockSize+d, j_node) -= aux_val * aux_CB_projection(d,j_node);
                            }
                        }
                    }

                    ///////// Taylor expansion from the surrogate boundary to the true boundary
                    const double distance = -0.5 * (
                        r_geom[sur_bd_local_ids[1]].FastGetSolutionStepValue(DISTANCE)
                        + r_geom[sur_bd_local_ids[2]].FastGetSolutionStepValue(DISTANCE));

                    // Shift direction: the level-set's own (CST-reconstructed) gradient
                    array_1d<double,2> grad_distance = ZeroVector(2);
                    for (std::size_t i_node = 0; i_node < 3; ++i_node) {
                        const double dist_i = r_geom[i_node].FastGetSolutionStepValue(DISTANCE);
                        grad_distance[0] += dist_i * DN_DX_parent(i_node, 0);
                        grad_distance[1] += dist_i * DN_DX_parent(i_node, 1);
                    }
                    const double grad_distance_norm = norm_2(grad_distance);
                    const array_1d<double,2> shift_dir = grad_distance / grad_distance_norm;

                    Matrix aux_N_taylor = ZeroMatrix(6, 18);
                    array_1d<double,3> r_sur_bd_N_taylor = ZeroVector(3);
                    r_sur_bd_N_taylor[0] = 0.0 + distance * (DN_DX_parent(0,0)*shift_dir[0] + DN_DX_parent(0,1)*shift_dir[1]);
                    r_sur_bd_N_taylor[1] = r_sur_bd_N(0,0) + distance * (DN_DX_parent(1,0)*shift_dir[0] + DN_DX_parent(1,1)*shift_dir[1]);
                    r_sur_bd_N_taylor[2] = r_sur_bd_N(0,1) + distance * (DN_DX_parent(2,0)*shift_dir[0] + DN_DX_parent(2,1)*shift_dir[1]);

                    // j_e := |e^ext|/|e~| 
                    const std::size_t je_n1 = sur_bd_local_ids[1];
                    const std::size_t je_n2 = sur_bd_local_ids[2];
                    const double je_dist_n1 = -r_geom[je_n1].FastGetSolutionStepValue(DISTANCE);
                    const double je_dist_n2 = -r_geom[je_n2].FastGetSolutionStepValue(DISTANCE);
                    const array_1d<double,3> je_shift_3d_dir = shift_dir[0]*local_e1 + shift_dir[1]*local_e2;
                    const array_1d<double,3> je_p1 = r_geom[je_n1].Coordinates();
                    const array_1d<double,3> je_p2 = r_geom[je_n2].Coordinates();
                    const array_1d<double,3> je_q1 = je_p1 + je_dist_n1 * je_shift_3d_dir;
                    const array_1d<double,3> je_q2 = je_p2 + je_dist_n2 * je_shift_3d_dir;
                    const double e_true_len = norm_2(je_q1 - je_q2);
                    const double e_surrogate_len = norm_2(je_p1 - je_p2);
                    const double j_e = e_true_len / e_surrogate_len;
                    const double aux_w_taylor = aux_w * j_e;


                    for (std::size_t i_node = 0; i_node < 3; ++i_node) {
                        aux_val = aux_w * r_sur_bd_N_taylor[i_node];
                        const double aux_val_gap = aux_w_taylor * r_sur_bd_N_taylor[i_node];
                        i_loc_id = sur_bd_local_ids[i_node];
                        for (std::size_t d = 0; d < 6; ++d) {
                            if (constrained_dofs[d] == 0.0) {
                                continue;
                            }
                            for (std::size_t j_node = 0; j_node < NumNodes * 6; ++j_node) {
                                left_hand_side_taylor(i_loc_id*BlockSize+d, j_node) -= aux_val * aux_CB_projection(d,j_node);
                                left_hand_side_taylor_gap(i_loc_id*BlockSize+d, j_node) -= aux_val_gap * aux_CB_projection(d,j_node);
                            }
                        }
                    }

                    // Curvature-scaled Neumann-consistency term
                    BoundedMatrix<double,6,LocalSize> aux_CB_projection_neumann = ZeroMatrix(6, LocalSize);
                    {
                        const Matrix aux_CB_local = prod(D, B_parent);
                        const array_1d<double,2> n_true = shift_dir;
                        array_1d<double,2> tau;
                        tau[0] = -n_true[1];
                        tau[1] = n_true[0];
                        const double n_tilde_dot_tau = n_sur_bd[0]*tau[0] + n_sur_bd[1]*tau[1];
                        for (std::size_t j = 0; j < LocalSize; ++j) {
                            const double m1_tau = tau[0]*aux_CB_local(3,j) + tau[1]*aux_CB_local(5,j);
                            const double m2_tau = tau[0]*aux_CB_local(5,j) + tau[1]*aux_CB_local(4,j);
                            aux_CB_projection_neumann(4,j) = n_tilde_dot_tau * m1_tau;
                            aux_CB_projection_neumann(3,j) = n_tilde_dot_tau * m2_tau;
                        }
                    }

                    // Penalty
                    const double penalty_parameter = this->GetProperties().Has(PENALTY_COEFFICIENT)
                        ? this->GetProperties()[PENALTY_COEFFICIENT]
                        : 0.0;

                    Matrix rAuxMat = ZeroMatrix(6, 18); //Taylor
                    for (IndexType i_node = 0; i_node < 3; ++i_node) {
                        for (IndexType d = 0; d < 5; ++d) {
                            if (constrained_dofs[d] == 0.0) {
                                continue;
                            }
                            rAuxMat(d, i_node*6 + d) = r_sur_bd_N_taylor[i_node];
                        }
                    }

                    // Gap SBM Neumann-consistency term
                    Matrix rAuxMat_neumann_ungated = ZeroMatrix(6, 18);
                    for (IndexType i_node = 0; i_node < 3; ++i_node) {
                        for (IndexType d = 0; d < 5; ++d) {
                            rAuxMat_neumann_ungated(d, i_node*6 + d) = r_sur_bd_N_taylor[i_node];
                        }
                    }
                    left_hand_side_gap_neumann -= prod(trans(rAuxMat_neumann_ungated), aux_CB_projection_neumann) * aux_w_taylor;

                    // material-stiffness-norm scaling
                    const double rho_C = norm_frobenius(D);
                    double rho_membrane_sq = 0.0, rho_bending_sq = 0.0, rho_shear_sq = 0.0;
                    for (std::size_t i = 0; i < 3; ++i) {
                        for (std::size_t j = 0; j < 3; ++j) {
                            rho_membrane_sq += D(i,j)*D(i,j);
                        }
                    }
                    for (std::size_t i = 3; i < 6; ++i) {
                        for (std::size_t j = 3; j < 6; ++j) {
                            rho_bending_sq += D(i,j)*D(i,j);
                        }
                    }
                    for (std::size_t i = 6; i < 8; ++i) {
                        for (std::size_t j = 6; j < 8; ++j) {
                            rho_shear_sq += D(i,j)*D(i,j);
                        }
                    }
                    const double rho_membrane = std::sqrt(rho_membrane_sq);
                    const double rho_bending = std::sqrt(rho_bending_sq);
                    const double rho_shear = std::sqrt(rho_shear_sq);
                    Matrix rho_diag = ZeroMatrix(6, 6);
                    rho_diag(0,0) = rho_membrane;
                    rho_diag(1,1) = rho_membrane;
                    rho_diag(2,2) = rho_shear;
                    rho_diag(3,3) = rho_bending;
                    rho_diag(4,4) = rho_bending;
                    rho_diag(5,5) = rho_C; 

                    const Matrix aux_taylor_penalty_base = prod(trans(rAuxMat), Matrix(prod(rho_diag, rAuxMat)));
                    left_hand_side_taylor_penalty = aux_taylor_penalty_base * penalty_parameter*aux_w/ h_sur_bd;
                    left_hand_side_taylor_penalty_gap = aux_taylor_penalty_base * penalty_parameter*aux_w_taylor/ h_sur_bd;

                    // Taylor penalty forcing/residual
                    if (sbm_formulation_type == 1) {
                        array_1d<double,6> penalty_bc_values_5 = penalty_bc_values;
                        penalty_bc_values_5[5] = 0.0; // rz row of rAuxMat is always zero, keep it inert
                        right_hand_side += (penalty_parameter*aux_w/h_sur_bd) * prod(trans(rAuxMat), Vector(prod(rho_diag, penalty_bc_values_5)));
                        right_hand_side -= prod(left_hand_side_taylor_penalty, unknown_values);
                    }

                    // Gap SBM domain term
                    const double local_x[3] = {data.LCS0.X1(), data.LCS0.X2(), data.LCS0.X3()};
                    const double local_y[3] = {data.LCS0.Y1(), data.LCS0.Y2(), data.LCS0.Y3()};
                    const std::size_t edge_n1 = sur_bd_local_ids[1];
                    const std::size_t edge_n2 = sur_bd_local_ids[2];
                    const double dist_n1 = -r_geom[edge_n1].FastGetSolutionStepValue(DISTANCE);
                    const double dist_n2 = -r_geom[edge_n2].FastGetSolutionStepValue(DISTANCE);
                    const double p1x = local_x[edge_n1], p1y = local_y[edge_n1];
                    const double p2x = local_x[edge_n2], p2y = local_y[edge_n2];
                    const double q1x = p1x + dist_n1 * shift_dir[0];
                    const double q1y = p1y + dist_n1 * shift_dir[1];
                    const double q2x = p2x + dist_n2 * shift_dir[0];
                    const double q2y = p2y + dist_n2 * shift_dir[1];
                    const double quad_area = 0.5 * std::abs(
                        p1x*q1y - q1x*p1y +
                        q1x*q2y - q2x*q1y +
                        q2x*p2y - p2x*q2y +
                        p2x*p1y - p1x*p2y);
                    const double edge_len = std::sqrt((p2x-p1x)*(p2x-p1x) + (p2y-p1y)*(p2y-p1y));
                    const double H_e = quad_area / edge_len;

                    Matrix BTD = Matrix(18, 8, 0.0);
                    BTD = prod(trans(B_parent), D);
                    Matrix ghost_block = prod(BTD, B_parent);
                    left_hand_side_gap_SBM += ghost_block * H_e * aux_w;
                }
        }

        // Assemble contributions according to the selected surrogate-boundary flux formulation.
        if (sbm_formulation_type == 0) {
            // MLS SBM: primal-consistency flux term
            rLeftHandSideMatrix += left_hand_side;
            rRightHandSideVector += right_hand_side;
        } else if (sbm_formulation_type == 1) {
            // Taylor SBM
            rLeftHandSideMatrix += left_hand_side;
            rLeftHandSideMatrix -= trans(left_hand_side_taylor);
            rLeftHandSideMatrix += left_hand_side_taylor_penalty;
            rRightHandSideVector += right_hand_side;
        } else {
            // Gap SBM
            Matrix left_hand_side_gap_sbm_total = left_hand_side_gap_SBM + left_hand_side_taylor_gap
                + trans(left_hand_side_taylor_gap) + left_hand_side_taylor_penalty_gap + left_hand_side_gap_neumann;
            rLeftHandSideMatrix += left_hand_side_gap_sbm_total;
            rRightHandSideVector -= prod(left_hand_side_gap_sbm_total, unknown_values);
            rRightHandSideVector += right_hand_side_shear_neumann;
        }
    }

    KRATOS_CATCH("")
}

template <ShellKinematics TKinematics>
void ShellThickShiftedBoundaryElement3D3N<TKinematics>::CalculateLocalSystemCutFEM(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geom = this->GetGeometry();
    array_1d<double, 3> dist;
    for (std::size_t i = 0; i < 3; ++i) {
        dist[i] = r_geom[i].FastGetSolutionStepValue(DISTANCE);
    }
    const std::size_t n_pos = (dist[0] > 0.0 ? 1 : 0) + (dist[1] > 0.0 ? 1 : 0) + (dist[2] > 0.0 ? 1 : 0);

    // Full-triangle domain stiffness. 
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    if (n_pos == 0) {
        // Fully outside: zero everything 
        rLeftHandSideMatrix = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
        rRightHandSideVector = ZeroVector(rRightHandSideVector.size());
        return;
    }
    if (n_pos == 3) {
        // Fully inside: full domain contribution, 
        return;
    }

    double full_area;
    array_1d<double, NumNodes> N_parent;
    BoundedMatrix<double, NumNodes, 2> DN_DX_parent;
    GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, full_area);

    // Exact zero-crossing on an edge whose two endpoints have opposite-sign distance.
    auto edge_intersection = [&](std::size_t a, std::size_t b) -> array_1d<double, 2> {
        const double t = dist[a] / (dist[a] - dist[b]);
        array_1d<double, 2> p;
        p[0] = r_geom[a].X() + t * (r_geom[b].X() - r_geom[a].X());
        p[1] = r_geom[a].Y() + t * (r_geom[b].Y() - r_geom[a].Y());
        return p;
    };

    std::vector<array_1d<double, 2>> cut_points;
    const std::size_t edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (const auto& e : edges) {
        if ((dist[e[0]] > 0.0) != (dist[e[1]] > 0.0)) {
            cut_points.push_back(edge_intersection(e[0], e[1]));
        }
    }
    KRATOS_ERROR_IF(cut_points.size() != 2)
        << "CutFEM SBM: expected exactly 2 cut points, got " << cut_points.size()
        << " in element " << this->Id() << std::endl;

    const array_1d<double, 2>& P1 = cut_points[0];
    const array_1d<double, 2>& P2 = cut_points[1];

    auto tri_area = [](double px, double py, double qx, double qy, double rx, double ry) {
        return 0.5 * std::abs(px * (qy - ry) + qx * (ry - py) + rx * (py - qy));
    };

    double trimmed_area;
    if (n_pos == 1) {
        // Physical region = the small triangle at the single positive node.
        const std::size_t pos_node = (dist[0] > 0.0) ? 0 : ((dist[1] > 0.0) ? 1 : 2);
        trimmed_area = tri_area(r_geom[pos_node].X(), r_geom[pos_node].Y(), P1[0], P1[1], P2[0], P2[1]);
    } else {
        // n_pos == 2: physical region = full triangle minus the small triangle at the
        // single negative node.
        const std::size_t neg_node = (dist[0] < 0.0) ? 0 : ((dist[1] < 0.0) ? 1 : 2);
        const double cut_off_area = tri_area(r_geom[neg_node].X(), r_geom[neg_node].Y(), P1[0], P1[1], P2[0], P2[1]);
        trimmed_area = full_area - cut_off_area;
    }

    const double area_fraction = trimmed_area / full_area;
    rRightHandSideVector *= area_fraction;

    typename BaseType::CalculationData data(this->mpCoordinateTransformation, rCurrentProcessInfo);
    data.CalculateLHS = true;
    data.CalculateRHS = true;
    BaseType::InitializeCalculationData(data);

    // Plain trimmed-area rescale
    rLeftHandSideMatrix *= area_fraction;

    // Same shear-consistency override the other formulations use for the boundary/flux
    // term: DSG's shear rows are an assumed-natural-strain construction 
    if (data.shearFormulation == 0) {
        const double x12 = data.LCS0.X1() - data.LCS0.X2();
        const double x23 = data.LCS0.X2() - data.LCS0.X3();
        const double x31 = data.LCS0.X3() - data.LCS0.X1();
        const double x21 = -x12, x32 = -x23, x13 = -x31;
        const double y12 = data.LCS0.Y1() - data.LCS0.Y2();
        const double y23 = data.LCS0.Y2() - data.LCS0.Y3();
        const double y31 = data.LCS0.Y3() - data.LCS0.Y1();
        const double A2 = 2.0 * data.TotalArea;
        for (std::size_t j = 0; j < LocalSize; ++j) {
            data.B(6, j) = 0.0;
            data.B(7, j) = 0.0;
        }
        data.B(6, 2) = y23 / A2;  data.B(6, 4) = 1.0 / 3.0;
        data.B(7, 2) = x32 / A2;  data.B(7, 3) = -1.0 / 3.0;
        data.B(6, 8) = y31 / A2;  data.B(6, 10) = 1.0 / 3.0;
        data.B(7, 8) = x13 / A2;  data.B(7, 9) = -1.0 / 3.0;
        data.B(6, 14) = y12 / A2; data.B(6, 16) = 1.0 / 3.0;
        data.B(7, 14) = x21 / A2; data.B(7, 15) = -1.0 / 3.0;
    }

    array_1d<double, 3>& gp0 = data.gpLocations[0];
    gp0[0] = 1.0 / 3.0; gp0[1] = 1.0 / 3.0; gp0[2] = 1.0 / 3.0;
    data.gpIndex = 0;
    noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);
    BaseType::CalculateSectionResponse(data);

    // Exact outward normal: distance is linear over this flat CST triangle
    array_1d<double, 2> grad_d;
    grad_d[0] = dist[0] * DN_DX_parent(0, 0) + dist[1] * DN_DX_parent(1, 0) + dist[2] * DN_DX_parent(2, 0);
    grad_d[1] = dist[0] * DN_DX_parent(0, 1) + dist[1] * DN_DX_parent(1, 1) + dist[2] * DN_DX_parent(2, 1);
    const double grad_norm = norm_2(grad_d);
    array_1d<double, 2> n_cut;
    n_cut[0] = -grad_d[0] / grad_norm;
    n_cut[1] = -grad_d[1] / grad_norm;

    const double seg_dx = P2[0] - P1[0];
    const double seg_dy = P2[1] - P1[1];
    const double seg_len = std::sqrt(seg_dx * seg_dx + seg_dy * seg_dy);
    const double mid_x = 0.5 * (P1[0] + P2[0]);
    const double mid_y = 0.5 * (P1[1] + P2[1]);

    // Shape function values of the ORIGINAL (uncut) triangle at the segment midpoint
    array_1d<double, 3> N_mid;
    N_mid[0] = tri_area(mid_x, mid_y, r_geom[1].X(), r_geom[1].Y(), r_geom[2].X(), r_geom[2].Y()) / full_area;
    N_mid[1] = tri_area(mid_x, mid_y, r_geom[2].X(), r_geom[2].Y(), r_geom[0].X(), r_geom[0].Y()) / full_area;
    N_mid[2] = 1.0 - N_mid[0] - N_mid[1];

    // aux_CB_projection/cauchy_traction: constant over the whole flat CST triangle 
    array_1d<double, 6> cauchy_traction;
    BoundedMatrix<double, 6, LocalSize> aux_CB_projection = ZeroMatrix(6, LocalSize);
    const auto& r_stress = data.generalizedStresses;
    const auto& r_C = data.D;

    BoundedMatrix<double, 8, LocalSize> B = ZeroMatrix(8, 18);
    for (IndexType i = 0; i < 3; ++i) {
        const IndexType initial_index = i * 6;
        B(0, initial_index) = DN_DX_parent(i, 0);
        B(1, initial_index + 1) = DN_DX_parent(i, 1);
        B(2, initial_index) = DN_DX_parent(i, 1);
        B(2, initial_index + 1) = DN_DX_parent(i, 0);
        B(3, initial_index + 4) = DN_DX_parent(i, 0);
        B(4, initial_index + 3) = -DN_DX_parent(i, 1);
        B(5, initial_index + 3) = -DN_DX_parent(i, 0);
        B(5, initial_index + 4) = DN_DX_parent(i, 1);
        B(6, initial_index + 2) = DN_DX_parent(i, 0);
        B(6, initial_index + 4) = N_parent[i];
        B(7, initial_index + 2) = DN_DX_parent(i, 1);
        B(7, initial_index + 3) = -N_parent[i];
    }

    Matrix D(8, 8);
    noalias(D) = ZeroMatrix(8, 8);
    for (std::size_t i = 0; i < r_C.size1(); ++i) {
        for (std::size_t j = 0; j < r_C.size1(); ++j) {
            D(i, j) = r_C(i, j);
        }
    }

    BoundedMatrix<double, 8, LocalSize> B_temp = ZeroMatrix(8, 18);
    BoundedMatrix<double, 8, LocalSize> B_parent = ZeroMatrix(8, 18);
    MatrixType T(18, 18);
    data.LCS.ComputeTotalRotationMatrix(T);
    B_temp = prod(data.B, T);
    Matrix R(8, 8);
    this->mSections[0]->GetRotationMatrixForGeneralizedStrains(-(this->mSections[0]->GetOrientationAngle()), R);
    R(6, 7) = -1.0 * R(6, 7);
    R(7, 6) = -1.0 * R(7, 6);
    B_parent = prod(R, B_temp);

    const array_1d<double,3> local_e1 = -data.LCS.Vx();
    const array_1d<double,3> local_e2 = -data.LCS.Vy();
    const array_1d<double,3> local_e3 = data.LCS.Vz();
    CalculateCBProjectionLinearisation(D, B_parent, n_cut, local_e1, local_e2, local_e3, aux_CB_projection);
    CalculateCauchyTractionVector(r_stress, n_cut, local_e1, local_e2, local_e3, cauchy_traction);

    // Per-DOF constrained mask (ux,uy,uz,rx,ry,rz)
    array_1d<double, 6> constrained_dofs;
    if (this->GetProperties().Has(SBM_CONSTRAINED_DOFS)) {
        const auto& r_constrained_dofs = this->GetProperties()[SBM_CONSTRAINED_DOFS];
        KRATOS_ERROR_IF(r_constrained_dofs.size() != 6)
            << "SBM_CONSTRAINED_DOFS must have size 6 (ux,uy,uz,rx,ry,rz) in Properties " << this->GetProperties().Id() << std::endl;
        for (std::size_t d = 0; d < 6; ++d) {
            constrained_dofs[d] = r_constrained_dofs[d];
        }
    } else {
        for (std::size_t d = 0; d < 6; ++d) {
            constrained_dofs[d] = 1.0;
        }
    }

    const double penalty_parameter = this->GetProperties().Has(PENALTY_COEFFICIENT)
        ? this->GetProperties()[PENALTY_COEFFICIENT]
        : 0.0;
    const double rho_C = norm_frobenius(D);
    const double h_e = std::sqrt(full_area);
    const double aux_weight_penalty = seg_len * penalty_parameter * rho_C / h_e;

    const auto& r_penalty_bc_val = this->GetValue(DISPLACEMENT);
    const auto& r_penalty_bc_val_rot = this->GetValue(ROTATION);
    array_1d<double, 6> bc_values;
    for (std::size_t d = 0; d < 3; ++d) {
        bc_values[d] = r_penalty_bc_val[d];
        bc_values[d + 3] = r_penalty_bc_val_rot[d];
    }

    Vector unknown_values = ZeroVector(LocalSize);
    for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
        const auto& r_disp = r_geom[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        const auto& r_rot = r_geom[i_node].FastGetSolutionStepValue(ROTATION);
        for (std::size_t d = 0; d < 3; ++d) {
            unknown_values(i_node * BlockSize + d) = r_disp[d];
            unknown_values(i_node * BlockSize + d + 3) = r_rot[d];
        }
    }

    Matrix left_hand_side_cut = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
    Vector right_hand_side_cut = ZeroVector(rRightHandSideVector.size());

    for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
        const double aux_val = seg_len * N_mid[i_node];
        for (std::size_t d = 0; d < 6; ++d) {
            if (constrained_dofs[d] == 0.0) {
                continue;
            }
            for (std::size_t j_node = 0; j_node < LocalSize; ++j_node) {
                // Term B (primal consistency)
                left_hand_side_cut(i_node * BlockSize + d, j_node) -= aux_val * aux_CB_projection(d, j_node);
                // Term C (adjoint consistency)
                left_hand_side_cut(j_node, i_node * BlockSize + d) -= aux_val * aux_CB_projection(d, j_node);
                // Term C forcing
                right_hand_side_cut(j_node) -= aux_val * aux_CB_projection(d, j_node) * bc_values[d];
            }
        }
        // Term D (penalty)
        for (std::size_t j_node = 0; j_node < NumNodes; ++j_node) {
            const double aux_2 = aux_weight_penalty * N_mid[i_node] * N_mid[j_node];
            for (std::size_t d = 0; d < 6; ++d) {
                if (constrained_dofs[d] == 0.0) {
                    continue;
                }
                left_hand_side_cut(i_node * BlockSize + d, j_node * BlockSize + d) += aux_2;
            }
        }
        for (std::size_t d = 0; d < 6; ++d) {
            if (constrained_dofs[d] == 0.0) {
                continue;
            }
            right_hand_side_cut(i_node * BlockSize + d) += aux_weight_penalty * N_mid[i_node] * bc_values[d];
        }
    }
    right_hand_side_cut -= prod(left_hand_side_cut, unknown_values);

    rLeftHandSideMatrix += left_hand_side_cut;
    rRightHandSideVector += right_hand_side_cut;

    KRATOS_CATCH("")
}

template <ShellKinematics TKinematics>
void ShellThickShiftedBoundaryElement3D3N<TKinematics>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // Delegate to CalculateLocalSystem so this can never drift out of sync with it.
    VectorType dummy_rhs;
    CalculateLocalSystem(rLeftHandSideMatrix, dummy_rhs, rCurrentProcessInfo);
    KRATOS_CATCH("")
}

template <ShellKinematics TKinematics>
void ShellThickShiftedBoundaryElement3D3N<TKinematics>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // Delegate to CalculateLocalSystem so this can never drift out of sync with it.
    MatrixType dummy_lhs;
    CalculateLocalSystem(dummy_lhs, rRightHandSideVector, rCurrentProcessInfo);
    KRATOS_CATCH("")
}

template <ShellKinematics TKinematics>
int ShellThickShiftedBoundaryElement3D3N<TKinematics>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    return BaseType::Check(rCurrentProcessInfo);
}

template <ShellKinematics TKinematics>
std::vector<std::size_t> ShellThickShiftedBoundaryElement3D3N<TKinematics>::GetSurrogateFacesIds()
{
    const std::size_t n_faces = 2 + 1;
    auto& r_neigh_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

    // Check the current element faces
    // Note that we rely on the fact that the neighbours are sorted according to the faces
    std::vector<std::size_t> surrogate_faces_ids;
    for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
        auto p_neigh_elem = r_neigh_elems(i_face).get();
        if (p_neigh_elem != nullptr && p_neigh_elem->Is(BOUNDARY)) {
            surrogate_faces_ids.push_back(i_face);
        }
    }

    return surrogate_faces_ids;
}

template <ShellKinematics TKinematics>
void ShellThickShiftedBoundaryElement3D3N<TKinematics>::CalculateCauchyTractionVector(
    const Vector& rVoigtStress,
    const array_1d<double,2>& rUnitNormal,
    const array_1d<double,3>& rLocalE1,
    const array_1d<double,3>& rLocalE2,
    const array_1d<double,3>& rLocalE3,
    array_1d<double,6>& rCauchyTraction)
{
    // rUnitNormal is the LOCAL (data.LCS0 in-plane) conormal.
    const double t1 = rVoigtStress[0]*rUnitNormal[0] + rVoigtStress[2]*rUnitNormal[1];
    const double t2 = rVoigtStress[2]*rUnitNormal[0] + rVoigtStress[1]*rUnitNormal[1];
    const double m1 = rVoigtStress[3]*rUnitNormal[0] + rVoigtStress[5]*rUnitNormal[1];
    const double m2 = rVoigtStress[5]*rUnitNormal[0] + rVoigtStress[4]*rUnitNormal[1];
    const double q  = rVoigtStress[6]*rUnitNormal[0] + rVoigtStress[7]*rUnitNormal[1];

    for (std::size_t k = 0; k < 3; ++k) {
        rCauchyTraction[k]     = t1*rLocalE1[k] + t2*rLocalE2[k] + q*rLocalE3[k];
        rCauchyTraction[3 + k] = m2*rLocalE1[k] + m1*rLocalE2[k];
    }
}

template <ShellKinematics TKinematics>
void ShellThickShiftedBoundaryElement3D3N<TKinematics>::CalculateCBProjectionLinearisation(
    const Matrix& rC,
    const BoundedMatrix<double,8,LocalSize>& rB,
    const array_1d<double,2>& rUnitNormal,
    const array_1d<double,3>& rLocalE1,
    const array_1d<double,3>& rLocalE2,
    const array_1d<double,3>& rLocalE3,
    BoundedMatrix<double,6,LocalSize>& rAuxMat)
{
    Matrix aux_CB = ZeroMatrix(8, LocalSize);
    aux_CB = prod(rC, rB);

    for (std::size_t j = 0; j < LocalSize; ++j) {
        const double t1 = rUnitNormal[0]*aux_CB(0,j) + rUnitNormal[1]*aux_CB(2,j);
        const double t2 = rUnitNormal[0]*aux_CB(2,j) + rUnitNormal[1]*aux_CB(1,j);
        const double m1 = rUnitNormal[0]*aux_CB(3,j) + rUnitNormal[1]*aux_CB(5,j);
        const double m2 = rUnitNormal[0]*aux_CB(5,j) + rUnitNormal[1]*aux_CB(4,j);
        const double q  = rUnitNormal[0]*aux_CB(6,j) + rUnitNormal[1]*aux_CB(7,j);

        for (std::size_t k = 0; k < 3; ++k) {
            rAuxMat(k, j)     = t1*rLocalE1[k] + t2*rLocalE2[k] + q*rLocalE3[k];
            rAuxMat(3 + k, j) = m2*rLocalE1[k] + m1*rLocalE2[k];
        }
    }
}

// template class ShellThickShiftedBoundaryElement3D3N<ShellKinematics::LINEAR>;
template class ShellThickShiftedBoundaryElement3D3N<ShellKinematics::NONLINEAR_COROTATIONAL>;

} // Namespace Kratos
