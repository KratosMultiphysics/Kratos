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

    // KRATOS_WATCH("HOW MANY?")

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

        Matrix left_hand_side_gap_SBM = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());

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

        // KRATOS_WATCH("Is interface");

        if (sur_bd_ids_vect.size() != 0) {

            // #define OPT_1_POINT_INTEGRATION
            typename BaseType::CalculationData data(this->mpCoordinateTransformation, rCurrentProcessInfo);
            data.CalculateLHS = true; //TO DO
            data.CalculateRHS = true; //TO DO

            BaseType::InitializeCalculationData(data);

            // for (SizeType i = 0; i < 3; ++i) {
            //     data.gpIndex = i;
                
                array_1d<double,3>& gp0 = data.gpLocations[0]; //TO DO
                gp0[0] = 1.0/3.0;
                gp0[1] = 1.0/3.0;
                gp0[2] = 1.0/3.0;

                data.gpIndex = 0;

                // calculate the total strain displ. matrix
                // BaseType::CalculateBMatrix(data);

                // // compute generalized strains
                // noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

                // calculate section response
                BaseType::CalculateSectionResponse(data);
                // #undef OPT_1_POINT_INTEGRATION
                // Get the parent geometry data
                double dom_size_parent;
                const auto& r_geom = this->GetGeometry();

                // KRATOS_WATCH(r_geom)

                array_1d<double, NumNodes> N_parent;
                BoundedMatrix<double, NumNodes, 2> DN_DX_parent;
                GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, dom_size_parent);
                BoundedMatrix<double,8,LocalSize> B;
                // StructuralMechanicsElementUtilities::CalculateB(*this, DN_DX_parent, B); //I obtained B in the parent class
                const auto &r_boundaries = r_geom.GenerateBoundariesEntities();
                DenseMatrix<unsigned int> nodes_in_faces;
                r_geom.NodesInFaces(nodes_in_faces);

                // Calculate the stress at the element midpoint
                // Note that in here we are assuming constant strain kinematics
                // KinematicVariables kinematic_variables(StrainSize, 2, NumNodes);
                // ConstitutiveVariables constitutive_variables(StrainSize);
                // // const auto &r_int_pts = this->IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_1);
                // ConstitutiveLaw::Parameters cons_law_values(r_geom, this->GetProperties(), rCurrentProcessInfo);
                // auto& r_cons_law_options = cons_law_values.GetOptions();
                // // r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, this->UseElementProvidedStrain());
                // // r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
                // // r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
                // cons_law_values.SetStrainVector(constitutive_variables.StrainVector);
                // cons_law_values.SetStressVector(constitutive_variables.StressVector); //additional
                // cons_law_values.SetConstitutiveMatrix(constitutive_variables.D);

                // CalculateKinematicVariables(kinematic_variables, 0, GeometryData::IntegrationMethod::GI_GAUSS_1);
                // CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, 0, r_int_pts, this->GetStressMeasure(), this->IsElementRotated());

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
                    const auto& r_sur_bd_dN = r_sur_bd_geom.ShapeFunctionsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_1);

                    // Get the gradient of the node contrary to the surrogate face
                    // Note that this is used to calculate the normal as n = - DN_DX_cont_node / norm_2(DN_DX_cont_node)
                    // const BoundedVector<double,2> DN_DX_cont_node = row(DN_DX_parent, sur_bd_local_ids[0]);
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

                        // B(3, initial_index + 3) = DN_DX_parent(i, 0);
                        // B(4, initial_index + 4) = DN_DX_parent(i, 1);
                        // B(5, initial_index + 3) = DN_DX_parent(i, 1);
                        // B(5, initial_index + 4) = DN_DX_parent(i, 0);

                        B(3, initial_index + 4) = DN_DX_parent(i, 0);
                        B(4, initial_index + 3) = -DN_DX_parent(i, 1);
                        B(5, initial_index + 3) = -DN_DX_parent(i, 0);
                        B(5, initial_index + 4) = DN_DX_parent(i, 1);

                        B(6, initial_index + 2) = DN_DX_parent(i, 0);
                        B(6, initial_index + 4) = N_parent[i];
                        B(7, initial_index + 2) = DN_DX_parent(i, 1);
                        B(7, initial_index + 3) = -N_parent[i];

                        // B(6, initial_index + 2) = -DN_DX_parent(i, 1);
                        // B(6, initial_index + 3) = N_parent[i];
                        // B(7, initial_index + 2) = DN_DX_parent(i, 0);
                        // B(7, initial_index + 4) = N_parent[i];
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

                    // KRATOS_WATCH(R)
                    // KRATOS_WATCH(B_parent)
                    // KRATOS_WATCH(B)
                    // KRATOS_WATCH(-(this->mSections[0]->GetOrientationAngle()))
                    // exit(0);

                    CalculateCBProjectionLinearisation(D, B_parent, n_sur_bd, aux_CB_projection);
                    CalculateCauchyTractionVector(r_stress, n_sur_bd, cauchy_traction);
    
                    // Add the surrogate boundary flux contribution
                    // Note that the local face ids. are already taken into account in the assembly
                    // Note that the integration weight is calculated as 2 * Parent domain size * norm(DN_DX_cont_node)
                    double aux_val;
                    std::size_t i_loc_id;
                    const double aux_w = 2 * dom_size_parent / h_sur_bd;// * data.dA; //TO DO Check


                    //This loop is the problem
                    for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
                        aux_val = aux_w * r_sur_bd_N(0,i_node);
                        i_loc_id = sur_bd_local_ids[i_node + 1];
                        for (std::size_t d = 0; d < 6; ++d) { //3 or 6
                            if (constrained_dofs[d] == 0.0) {
                                continue;
                            }
                            right_hand_side(i_loc_id*BlockSize+d) += aux_val * cauchy_traction[d];
                            // for (std::size_t j_node = 0; j_node < NumNodes; ++j_node) {
                            for (std::size_t j_node = 0; j_node < NumNodes * 6; ++j_node) { //TO DO
                                left_hand_side(i_loc_id*BlockSize+d, j_node) -= aux_val * aux_CB_projection(d,j_node);
                                // rLeftHandSideMatrix(i_loc_id*BlockSize+d, j_node*BlockSize+d) -= aux_val * aux_CB_projection(d,j_node*BlockSize + d);
                            }
                        }
                    }

                    ///////// Taylor expansion from the surrogate boundary to the true boundary
                    const double distance = -0.5 * (
                        r_geom[sur_bd_local_ids[1]].FastGetSolutionStepValue(DISTANCE)
                        + r_geom[sur_bd_local_ids[2]].FastGetSolutionStepValue(DISTANCE));

                    Matrix aux_N_taylor = ZeroMatrix(6, 18);
                    array_1d<double,3> r_sur_bd_N_taylor = ZeroVector(3);
                    r_sur_bd_N_taylor[0] = 0.0 + (distance * DN_DX_parent(0,0));
                    r_sur_bd_N_taylor[1] = r_sur_bd_N(0,0)+ (distance * DN_DX_parent(1,0));
                    r_sur_bd_N_taylor[2] = r_sur_bd_N(0,1)+ (distance * DN_DX_parent(2,0));

                    // KRATOS_WATCH(r_sur_bd_N)
                    // KRATOS_WATCH(N_parent)
                    // KRATOS_WATCH(DN_DX_parent)
                    // // KRATOS_WATCH(r_geom)
                    // KRATOS_WATCH(sur_bd_local_ids)
                    // KRATOS_WATCH(r_sur_bd_N_taylor)
                    // KRATOS_WATCH(h_sur_bd)
                    // KRATOS_WATCH(distance)
                    // KRATOS_WATCH(distance_2)

                    // exit(0);

                    //This loop is for the sbm with taylor expansion.
                    for (std::size_t i_node = 0; i_node < 3; ++i_node) {
                        aux_val = aux_w * r_sur_bd_N_taylor[i_node];
                        i_loc_id = sur_bd_local_ids[i_node];
                        for (std::size_t d = 0; d < 6; ++d) { //3 or 6
                            if (constrained_dofs[d] == 0.0) {
                                continue;
                            }
                            for (std::size_t j_node = 0; j_node < NumNodes * 6; ++j_node) {
                                left_hand_side_taylor(i_loc_id*BlockSize+d, j_node) -= aux_val * aux_CB_projection(d,j_node);
                            }
                        }
                    }

                    // KRATOS_WATCH(aux_val)
                    // KRATOS_WATCH(aux_CB_projection)
                    // KRATOS_WATCH(left_hand_side_taylor)
                    // KRATOS_WATCH(left_hand_side)
                    // exit(0);    


                     // Penalty
                    double aux_val_penalty;
                    const double penalty_parameter = this->GetProperties()[PENALTY_COEFFICIENT];

                    // KRATOS_WATCH(aux_w)
                    // exit(0);
                    // array_1d<double,3> r_sur_bd_N_penalty = ZeroVector(3);
                    // r_sur_bd_N_penalty[0] = 0.0 ;
                    // r_sur_bd_N_penalty[1] = r_sur_bd_N(0,0);
                    // r_sur_bd_N_penalty[2] = r_sur_bd_N(0,1);

                    // for (std::size_t i_node = 0; i_node < 3; ++i_node) {
                    //     aux_val_penalty = aux_w * r_sur_bd_N_penalty[i_node] / h_sur_bd; // different than Nitsche
                    //     i_loc_id = sur_bd_local_ids[i_node ];
                    //     for (std::size_t d = 0; d < 5; ++d) { //HACK! if I put 6, it will affect the solution. 5: fixed, 3: SS
                    //         for (std::size_t j_node = 0; j_node < NumNodes; ++j_node) { //HACK!
                    //             left_hand_side_taylor_penalty(i_loc_id*BlockSize+d, j_node*BlockSize+d) += r_sur_bd_N_penalty[i_node] * aux_val_penalty * penalty_parameter; //Penalty
                    //         }
                    //     }
                    // }

                    Matrix rAuxMat = ZeroMatrix(6, 18); //Taylor
                    Matrix trans_rAuxMat = ZeroMatrix(6, 18);

                    array_1d<double,3> r_sur_bd_N_= ZeroVector(3);
                    r_sur_bd_N_[0] = 0.0;;
                    r_sur_bd_N_[1] = r_sur_bd_N(0,0);
                    r_sur_bd_N_[2] = r_sur_bd_N(0,1);

                    for (IndexType i_node = 0; i_node < 3; ++i_node) {
                        for (IndexType d = 0; d < 5; ++d) {  //3 or 6
                            if (constrained_dofs[d] == 0.0) {
                                continue;
                            }
                            rAuxMat(d, i_node*6 + d) = r_sur_bd_N_taylor[i_node];
                            trans_rAuxMat(d, i_node*6 + d) = r_sur_bd_N_[i_node];
                        }
                    }

                    // KRATOS_WATCH(rAuxMat)
                    // KRATOS_WATCH(trans_rAuxMat)

                    left_hand_side_taylor_penalty = prod(trans(rAuxMat), rAuxMat)* penalty_parameter*aux_w/ h_sur_bd;

                    // Taylor penalty forcing/residual 
                    if (sbm_formulation_type == 1) {
                        array_1d<double,6> penalty_bc_values_5 = penalty_bc_values;
                        penalty_bc_values_5[5] = 0.0; // rz row of rAuxMat is always zero, keep it inert
                        right_hand_side += (penalty_parameter*aux_w/h_sur_bd) * prod(trans(rAuxMat), penalty_bc_values_5);
                        right_hand_side -= prod(left_hand_side_taylor_penalty, unknown_values);
                    }

                    // exit(0);
                    
                    // Gap SBM domain term
                    Matrix BTD = Matrix(18, 8, 0.0);
                    BTD = prod(trans(B_parent), D);
                    left_hand_side_gap_SBM += prod(BTD, B_parent) * std::abs(distance) * aux_w;

                    // KRATOS_WATCH(left_hand_side)
                    // KRATOS_WATCH(aux_CB_projection)
                    // exit(0);
                                
                    // // KRATOS_WATCH(n_sur_bd)
                    // KRATOS_WATCH(aux_CB_projection)        
                    // // KRATOS_WATCH(sur_bd_ids_vect.size())
                    // KRATOS_WATCH(aux_w)
                    // KRATOS_WATCH(left_hand_side)
                    // // KRATOS_WATCH(r_sur_bd_N)
                    // // KRATOS_WATCH(n_bd_points)
                    // exit(0);
                }
            // }
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
            Matrix left_hand_side_gap_sbm_total = left_hand_side_gap_SBM + left_hand_side_taylor
                + trans(left_hand_side_taylor) + left_hand_side_taylor_penalty;
            rLeftHandSideMatrix += left_hand_side_gap_sbm_total;
            rRightHandSideVector -= prod(left_hand_side_gap_sbm_total, unknown_values);
        }
    }

    KRATOS_CATCH("")
}

template <ShellKinematics TKinematics>
void ShellThickShiftedBoundaryElement3D3N<TKinematics>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Add base Laplacian contribution
    BaseType::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    // // Check if the element belongs to the surrogate interface
    // // Note that the INTERFACE flag is assumed to be set in the layer of elements attached to the surrogate interface
    // if (this->Is(INTERFACE)) {
    //     // Find the surrogate face local id
    //     // Note that it might happen that an interface element has no surrogate face (i.e. a unique node in the surrogate skin)
    //     const auto sur_bd_ids_vect = GetSurrogateFacesIds();
    //     if (sur_bd_ids_vect.size() != 0) {
    //         // Get the parent geometry data
    //         double dom_size_parent;
    //         const auto& r_geom = this->GetGeometry();
    //         array_1d<double, NumNodes> N_parent;
    //         BoundedMatrix<double, NumNodes, 2> DN_DX_parent;
    //         GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, dom_size_parent);
    //         BoundedMatrix<double,StrainSize,LocalSize> B;
    //         StructuralMechanicsElementUtilities::CalculateB(*this, DN_DX_parent, B);
    //         const auto &r_boundaries = r_geom.GenerateBoundariesEntities();
    //         DenseMatrix<unsigned int> nodes_in_faces;
    //         r_geom.NodesInFaces(nodes_in_faces);

    //         // Calculate the stress at the element midpoint
    //         // Note that in here we are assuming constant strain kinematics
    //         KinematicVariables kinematic_variables(StrainSize, 2, NumNodes);
    //         ConstitutiveVariables constitutive_variables(StrainSize);
    //         // const auto &r_int_pts = this->IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_1);
    //         ConstitutiveLaw::Parameters cons_law_values(r_geom, this->GetProperties(), rCurrentProcessInfo);
    //         auto& r_cons_law_options = cons_law_values.GetOptions();
    //         // r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, this->UseElementProvidedStrain());
    //         r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    //         r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    //         cons_law_values.SetStrainVector(constitutive_variables.StrainVector);
    //         CalculateKinematicVariables(kinematic_variables, 0, GeometryData::IntegrationMethod::GI_GAUSS_1);

    //         // exit(0); 

    //         // CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, 0, r_int_pts, this->GetStressMeasure(), this->IsElementRotated());

    //         // Loop the surrogate faces
    //         // Note that there is the chance that the surrogate face is not unique
    //         for (std::size_t sur_bd_id : sur_bd_ids_vect) {
    //             // Get the current surrogate face geometry information
    //             const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
    //             const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
    //             const DenseVector<std::size_t> sur_bd_local_ids = column(nodes_in_faces, sur_bd_id);
    //             const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

    //             // Get the gradient of the node contrary to the surrogate face
    //             // Note that this is used to calculate the normal as n = - DN_DX_cont_node / norm_2(DN_DX_cont_node)
    //             // const BoundedVector<double,2> DN_DX_cont_node = row(DN_DX_parent, sur_bd_local_ids[0]);
    //             BoundedVector<double,2> n_sur_bd = row(DN_DX_parent, sur_bd_local_ids[0]);
    //             const double h_sur_bd = 1.0 / norm_2(n_sur_bd);
    //             n_sur_bd *= -h_sur_bd;

    //             // Calculate the required projections
    //             array_1d<double,2> cauchy_traction;
    //             BoundedMatrix<double,2,LocalSize> aux_CB_projection;
    //             const auto &r_stress = cons_law_values.GetStressVector();
    //             const auto &r_C = cons_law_values.GetConstitutiveMatrix();
    //             CalculateCauchyTractionVector(r_stress, n_sur_bd, cauchy_traction);
    //             CalculateCBProjectionLinearisation(r_C, B, n_sur_bd, aux_CB_projection);

    //             // Add the surrogate boundary flux contribution
    //             // Note that the local face ids. are already taken into account in the assembly
    //             // Note that the integration weight is calculated as 2 * Parent domain size * norm(DN_DX_cont_node)
    //             double aux_val;
    //             std::size_t i_loc_id;
    //             const double aux_w = 2 * dom_size_parent / h_sur_bd;
    //             for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
    //                 aux_val = aux_w * r_sur_bd_N(0,i_node);
    //                 i_loc_id = sur_bd_local_ids[i_node + 1];
    //                 for (std::size_t d = 0; d < 2; ++d) {
    //                     for (std::size_t j_node = 0; j_node < NumNodes; ++j_node) {
    //                         rLeftHandSideMatrix(i_loc_id*BlockSize+d, j_node*BlockSize+d) -= aux_val * aux_CB_projection(d,j_node*BlockSize + d);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    KRATOS_CATCH("")
}

template <ShellKinematics TKinematics>
void ShellThickShiftedBoundaryElement3D3N<TKinematics>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Add base Laplacian contribution
    BaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Check if the element belongs to the surrogate interface
    // Note that the INTERFACE flag is assumed to be set in the layer of elements attached to the surrogate interface

    // TODO 

    // if (this->Is(INTERFACE)) {
    //     // Find the surrogate face local id
    //     // Note that it might happen that an interface element has no surrogate face (i.e. a unique node in the surrogate skin)
    //     const auto sur_bd_ids_vect = GetSurrogateFacesIds();
    //     if (sur_bd_ids_vect.size() != 0) {
    //         // Get the parent geometry data
    //         double dom_size_parent;
    //         const auto& r_geom = this->GetGeometry();
    //         array_1d<double, NumNodes> N_parent;
    //         BoundedMatrix<double, NumNodes, 2> DN_DX_parent;
    //         GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, dom_size_parent);
    //         BoundedMatrix<double,StrainSize,LocalSize> B;
    //         StructuralMechanicsElementUtilities::CalculateB(*this, DN_DX_parent, B);
    //         const auto &r_boundaries = r_geom.GenerateBoundariesEntities();
    //         DenseMatrix<unsigned int> nodes_in_faces;
    //         r_geom.NodesInFaces(nodes_in_faces);

    //         // Calculate the stress at the element midpoint
    //         // Note that in here we are assuming constant strain kinematics
    //         KinematicVariables kinematic_variables(StrainSize, 2, NumNodes);
    //         ConstitutiveVariables constitutive_variables(StrainSize);
    //         // const auto &r_int_pts = this->IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_1);
    //         ConstitutiveLaw::Parameters cons_law_values(r_geom, this->GetProperties(), rCurrentProcessInfo);
    //         auto& r_cons_law_options = cons_law_values.GetOptions();
    //         // r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, this->UseElementProvidedStrain());
    //         // r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    //         // r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    //         cons_law_values.SetStrainVector(constitutive_variables.StrainVector);
    //         // CalculateKinematicVariables(kinematic_variables, 0, GeometryData::IntegrationMethod::GI_GAUSS_1);
    //         // CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, 0, r_int_pts, this->GetStressMeasure(), this->IsElementRotated());

    //         // Loop the surrogate faces
    //         // Note that there is the chance that the surrogate face is not unique
    //         for (std::size_t sur_bd_id : sur_bd_ids_vect) {
    //             // Get the current surrogate face geometry information
    //             const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
    //             const unsigned int n_bd_points = r_sur_bd_geom.PointsNumber();
    //             const DenseVector<std::size_t> sur_bd_local_ids = column(nodes_in_faces, sur_bd_id);
    //             const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

    //             // Get the gradient of the node contrary to the surrogate face
    //             // Note that this is used to calculate the normal as n = - DN_DX_cont_node / norm_2(DN_DX_cont_node)
    //             // const BoundedVector<double,2> DN_DX_cont_node = row(DN_DX_parent, sur_bd_local_ids[0]);
    //             BoundedVector<double,2> n_sur_bd = row(DN_DX_parent, sur_bd_local_ids[0]);
    //             const double h_sur_bd = 1.0 / norm_2(n_sur_bd);
    //             n_sur_bd *= -h_sur_bd;

    //             // Calculate the required projections
    //             array_1d<double,2> cauchy_traction;
    //             BoundedMatrix<double,2,LocalSize> aux_CB_projection;
    //             const auto &r_stress = cons_law_values.GetStressVector();
    //             const auto &r_C = cons_law_values.GetConstitutiveMatrix();
    //             CalculateCauchyTractionVector(r_stress, n_sur_bd, cauchy_traction);
    //             CalculateCBProjectionLinearisation(r_C, B, n_sur_bd, aux_CB_projection);

    //             // Add the surrogate boundary flux contribution
    //             // Note that the local face ids. are already taken into account in the assembly
    //             // Note that the integration weight is calculated as 2 * Parent domain size * norm(DN_DX_cont_node)
    //             double aux_val;
    //             std::size_t i_loc_id;
    //             const double aux_w = 2 * dom_size_parent / h_sur_bd;
    //             for (std::size_t i_node = 0; i_node < n_bd_points; ++i_node) {
    //                 aux_val = aux_w * r_sur_bd_N(0,i_node);
    //                 i_loc_id = sur_bd_local_ids[i_node + 1];
    //                 for (std::size_t d = 0; d < 2; ++d) {
    //                     rRightHandSideVector(i_loc_id*BlockSize+d) += aux_val * cauchy_traction[d];
    //                 }
    //             }
    //         }
    //     }
    // }

    KRATOS_CATCH("")
}

template <ShellKinematics TKinematics>
void ShellThickShiftedBoundaryElement3D3N<TKinematics>::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const typename GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    // const auto& r_geometry = GetGeometry();

    // // const GeometryType::IntegrationPointsArrayType& r_integration_points = this->IntegrationPoints(rIntegrationMethod);
    // // Shape functions
    // rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    // // rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);

    // // KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    // // Compute B
    // // CalculateB( rThisKinematicVariables.B, rThisKinematicVariables.DN_DX, r_integration_points, PointNumber );

    // // Compute equivalent F
    // // GetValuesVector(rThisKinematicVariables.Displacements);
    // Vector strain_vector(mConstitutiveLawVector[0]->GetStrainSize());
    // noalias(strain_vector) = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);
    // // ComputeEquivalentF(rThisKinematicVariables.F, strain_vector);
    // rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

template <ShellKinematics TKinematics>
int ShellThickShiftedBoundaryElement3D3N<TKinematics>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Check that only simplicial geometries are used
    const auto geom_type = this->GetGeometry().GetGeometryType();
    // KRATOS_ERROR_IF_NOT(geom_type == GeometryData::KratosGeometryType::Kratos_Triangle2D3 || geom_type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) <<
    //     "ShellThickShiftedBoundaryElement3D3N only supports simplicial geometries (linear triangles and tetrahedra)." << std::endl;

    // Base SmallDisplacement element check
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
    array_1d<double,6>& rCauchyTraction)
{
    // if constexpr (TDim == 2) {
        rCauchyTraction[0] = rVoigtStress[0]*rUnitNormal[0] + rVoigtStress[2]*rUnitNormal[1];
        rCauchyTraction[1] = rVoigtStress[2]*rUnitNormal[0] + rVoigtStress[1]*rUnitNormal[1];
    // } else {
    //     rCauchyTraction[0] = rVoigtStress[0]*rUnitNormal[0] + rVoigtStress[3]*rUnitNormal[1] + rVoigtStress[5]*rUnitNormal[2];
    //     rCauchyTraction[1] = rVoigtStress[3]*rUnitNormal[0] + rVoigtStress[1]*rUnitNormal[1] + rVoigtStress[4]*rUnitNormal[2];
    //     rCauchyTraction[2] = rVoigtStress[5]*rUnitNormal[0] + rVoigtStress[4]*rUnitNormal[1] + rVoigtStress[2]*rUnitNormal[2];
    // }
}

template <ShellKinematics TKinematics>
void ShellThickShiftedBoundaryElement3D3N<TKinematics>::CalculateCBProjectionLinearisation(
    const Matrix& rC,
    const BoundedMatrix<double,8,LocalSize>& rB,
    const array_1d<double,2>& rUnitNormal,
    BoundedMatrix<double,6,LocalSize>& rAuxMat)
{

    Matrix aux_CB = ZeroMatrix(8, LocalSize);
    aux_CB = prod(rC, rB);

    // if constexpr (TDim == 2) {
    
        // membrane part
        for (std::size_t j = 0; j < LocalSize; ++j) {
            rAuxMat(0,j) = rUnitNormal[0]*aux_CB(0,j) + rUnitNormal[1]*aux_CB(2,j)*1;
        }
        for (std::size_t j = 0; j < LocalSize; ++j) {
            rAuxMat(1,j) = rUnitNormal[0]*aux_CB(2,j) + rUnitNormal[1]*aux_CB(1,j)*1;
        }

        // bending part    
        for (std::size_t j = 0; j < LocalSize; ++j) {
            rAuxMat(4,j) = (rUnitNormal[0]*aux_CB(3,j) + rUnitNormal[1]*aux_CB(5,j))*1;
        }
        for (std::size_t j = 0; j < LocalSize; ++j) {
            rAuxMat(3,j) = (rUnitNormal[0]*aux_CB(5,j) + rUnitNormal[1]*aux_CB(4,j))*1;
        }

        // // TO DO: shear and drilling part
        for (std::size_t j = 0; j < LocalSize; ++j) {
            rAuxMat(2,j) = rUnitNormal[0]*aux_CB(6,j) + rUnitNormal[1]*aux_CB(7,j);
        }


    // } else {
    //     for (std::size_t j = 0; j < LocalSize; ++j) {
    //         rAuxMat(0,j) = rUnitNormal[0]*aux_CB(0,j) + rUnitNormal[1]*aux_CB(3,j) + rUnitNormal[2]*aux_CB(5,j);
    //     }
    //     for (std::size_t j = 0; j < LocalSize; ++j) {
    //         rAuxMat(1,j) = rUnitNormal[0]*aux_CB(3,j) + rUnitNormal[1]*aux_CB(1,j) + rUnitNormal[2]*aux_CB(4,j);
    //     }
    //     for (std::size_t j = 0; j < LocalSize; ++j) {
    //         rAuxMat(2,j) = rUnitNormal[0]*aux_CB(5,j) + rUnitNormal[1]*aux_CB(4,j) + rUnitNormal[2]*aux_CB(2,j);
    //     }
    // }
}

// template class ShellThickShiftedBoundaryElement3D3N<ShellKinematics::LINEAR>;
template class ShellThickShiftedBoundaryElement3D3N<ShellKinematics::NONLINEAR_COROTATIONAL>;

} // Namespace Kratos
