//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andea Gorgi
//                  
//

// System includes

// External includes

// Project includes
#include "custom_conditions/sbm_contact_2D_condition.h"
#include "includes/global_pointer_variables.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <vector>

namespace Kratos
{

void SbmContact2DCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    InitializeMemberVariables();
    InitializeSbmMemberVariables();

    if (mpProjectionNodeMaster != nullptr) {
        const auto& r_geometry_master = GetMasterGeometry();
        std::ofstream output_file("txt_files/Contact_Projection_Coordinates.txt", std::ios::app);
        output_file << mpProjectionNodeMaster->X() << " " << mpProjectionNodeMaster->Y() << " "
                    << r_geometry_master.Center().X() << " " << r_geometry_master.Center().Y() << "\n";
    }
}


    void SbmContact2DCondition::InitializeMaterial()
    {
        KRATOS_TRY
        if ((*mpPropMaster)[CONSTITUTIVE_LAW] != nullptr ) {
            const GeometryType& r_geometry_master = GetMasterGeometry();

            const auto& N_values = r_geometry_master.ShapeFunctionsValues(this->GetIntegrationMethod());

            mpConstitutiveLawMaster = (*mpPropMaster)[CONSTITUTIVE_LAW]->Clone();
            mpConstitutiveLawMaster->InitializeMaterial( *mpPropMaster, r_geometry_master, row(N_values , 0 ));

        } else
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;


        if ( (*mpPropSlave)[CONSTITUTIVE_LAW] != nullptr ) {
            const GeometryType& r_geometry_slave = GetSlaveGeometry();

            const auto& N_values = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());

            mpConstitutiveLawSlave = (*mpPropSlave)[CONSTITUTIVE_LAW]->Clone();
            mpConstitutiveLawSlave->InitializeMaterial( *mpPropSlave, r_geometry_slave, row(N_values , 0 ));

        } else
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

        KRATOS_CATCH( "" );

    }

void SbmContact2DCondition::InitializeMemberVariables()
{
    const auto& r_geometry_master = GetMasterGeometry();
    const auto& r_geometry_slave = GetSlaveGeometry();

    const auto& r_DN_De_master = r_geometry_master.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const auto& r_DN_De_slave = r_geometry_slave.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

    mMasterDim = r_DN_De_master[0].size2();
    mSlaveDim = r_DN_De_slave[0].size2();

    mMasterBasisFunctionsOrder = static_cast<SizeType>(std::lround(std::sqrt(static_cast<double>(r_DN_De_master[0].size1())))) - 1;
    mSlaveBasisFunctionsOrder = static_cast<SizeType>(std::lround(std::sqrt(static_cast<double>(r_DN_De_slave[0].size1())))) - 1;

    if (mMasterBasisFunctionsOrder < 1) {
        mMasterBasisFunctionsOrder = 1;
    }
    if (mSlaveBasisFunctionsOrder < 1) {
        mSlaveBasisFunctionsOrder = 1;
    }

    if (this->Has(KNOT_SPAN_SIZES)) {
        const Vector& mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
        KRATOS_ERROR_IF(mesh_size_uv.size() == 0) << ":::[SbmContact2DCondition]::: Empty mesh-size vector on condition with Id "
                                                 << this->Id() << std::endl;

        double h = mesh_size_uv[0];
        for (IndexType i = 1; i < mesh_size_uv.size(); ++i) {
            h = std::min(h, mesh_size_uv[i]);
        }

        mMasterCharacteristicLength = h;
        mSlaveCharacteristicLength = h;
    } else {
        KRATOS_ERROR << ":::[SbmContact2DCondition]::: Missing mesh-size information (KNOT_SPAN_SIZES) on condition with Id "
                     << this->Id() << std::endl;
    }

    array_1d<double, 3> tangent;
    array_1d<double, 3> normal;

    r_geometry_master.Calculate(LOCAL_TANGENT, tangent);
    const double master_tangent_norm = norm_2(tangent);
    KRATOS_ERROR_IF(master_tangent_norm < std::numeric_limits<double>::epsilon())
        << ":::[SbmContact2DCondition]::: Master tangent norm is zero on condition with Id " << this->Id() << std::endl;
    tangent /= master_tangent_norm;

    normal[0] = tangent[1];
    normal[1] = -tangent[0];
    normal[2] = 0.0;

    mNormalMaster.resize(3);
    mNormalMaster[0] = normal[0];
    mNormalMaster[1] = normal[1];
    mNormalMaster[2] = normal[2];

    r_geometry_slave.Calculate(LOCAL_TANGENT, tangent);
    const double slave_tangent_norm = norm_2(tangent);
    KRATOS_ERROR_IF(slave_tangent_norm < std::numeric_limits<double>::epsilon())
        << ":::[SbmContact2DCondition]::: Slave tangent norm is zero on condition with Id " << this->Id() << std::endl;
    tangent /= slave_tangent_norm;

    normal[0] = tangent[1];
    normal[1] = -tangent[0];
    normal[2] = 0.0;

    mNormalSlave.resize(3);
    mNormalSlave[0] = normal[0];
    mNormalSlave[1] = normal[1];
    mNormalSlave[2] = normal[2];

    this->SetValue(NORMAL, mNormalMaster);
    this->SetValue(NORMAL_MASTER, mNormalMaster);
    this->SetValue(NORMAL_SLAVE, mNormalSlave);                                       
}

void SbmContact2DCondition::InitializeSbmMemberVariables()
{
    auto& r_geometry_master = GetMasterGeometry();
    auto& r_geometry_slave = GetSlaveGeometry();

    auto& master_neighbour_nodes = r_geometry_master.GetValue(NEIGHBOUR_NODES);
    KRATOS_ERROR_IF(master_neighbour_nodes.size() == 0) << ":::[SbmContact2DCondition]::: Missing master NEIGHBOUR_NODES on condition with Id "
                                                        << this->Id() << std::endl;

    mpProjectionNodeMaster = &master_neighbour_nodes[0];

    auto& slave_neighbour_nodes = mpProjectionNodeMaster->GetValue(NEIGHBOUR_NODES);
    KRATOS_ERROR_IF(slave_neighbour_nodes.size() == 0) << ":::[SbmContact2DCondition]::: Missing slave NEIGHBOUR_NODES on projection node "
                                                       << mpProjectionNodeMaster->Id() << " (condition Id " << this->Id() << ")." << std::endl;

    mpProjectionNodeSlave = &slave_neighbour_nodes[0];

    auto& slave_neighbour_geometries = mpProjectionNodeSlave->GetValue(NEIGHBOUR_GEOMETRIES);
    KRATOS_ERROR_IF(slave_neighbour_geometries.empty()) << ":::[SbmContact2DCondition]::: Missing neighbour geometries on projection node "
                                                        << mpProjectionNodeSlave->Id() << " (condition Id " << this->Id() << ")." << std::endl;
    slave_neighbour_geometries[0]->SetValue(IS_ASSEMBLED, false);
    mpProjectionNodeMaster->SetValue(IS_ASSEMBLED, false);

    if (mDistanceMaster.size() != mMasterDim) {
        mDistanceMaster.resize(mMasterDim, false);
    }
    if (mDistanceSlave.size() != mSlaveDim) {
        mDistanceSlave.resize(mSlaveDim, false);
    }

    noalias(mDistanceMaster) = ZeroVector(mMasterDim);
    noalias(mDistanceSlave) = ZeroVector(mSlaveDim);

    const auto& master_center = r_geometry_master.Center();
    const auto& slave_center = r_geometry_slave.Center();

    mDistanceMaster[0] = mpProjectionNodeMaster->X() - master_center.X();
    mDistanceMaster[1] = mpProjectionNodeMaster->Y() - master_center.Y();

    mDistanceSlave[0] = mpProjectionNodeSlave->X() - slave_center.X();
    mDistanceSlave[1] = mpProjectionNodeSlave->Y() - slave_center.Y();

    mTrueNormalMaster = mpProjectionNodeMaster->GetValue(NORMAL);
    mTrueNormalSlave = mpProjectionNodeSlave->GetValue(NORMAL);

    this->SetValue(NORMAL_SKIN_MASTER, mTrueNormalMaster);
    this->SetValue(NORMAL_SKIN_SLAVE, mTrueNormalSlave);

    mMasterTrueDotSurrogateNormal = inner_prod(mNormalMaster, mTrueNormalMaster);
    mSlaveTrueDotSurrogateNormal = inner_prod(mNormalSlave, mTrueNormalSlave);

    KRATOS_ERROR_IF(mMasterTrueDotSurrogateNormal <= 0.0) << ":::[SbmContact2DCondition]::: Error. Master: Negative n_ntilde_master"
                                                          << " (condition Id " << this->Id() << ")." << std::endl;
    KRATOS_ERROR_IF(mSlaveTrueDotSurrogateNormal <= 0.0) << ":::[SbmContact2DCondition]::: Error. Slave: Negative n_ntilde_slave"
                                                         << " (condition Id " << this->Id() << ")." << std::endl;

    //FIXME: debug contact
    this->SetValue(PROJECTION_NODE_COORDINATES, mpProjectionNodeMaster->Coordinates());
    this->SetValue(PROJECTION_NODE_MASTER_ID, mpProjectionNodeMaster->Id());
    this->SetValue(PROJECTION_NODE_SLAVE_ID, mpProjectionNodeSlave->Id());

    double curvature_true_master = mpProjectionNodeMaster->GetValue(CURVATURE);
    double curvature_true_slave = mpProjectionNodeSlave->GetValue(CURVATURE);

    this->SetValue(CURVATURE_PROJECTION_MASTER, curvature_true_master);
    this->SetValue(CURVATURE_PROJECTION_SLAVE, curvature_true_slave);



    // //------------------- REFERENCE WEIGHT

    const double thickness_master = (*mpPropMaster).Has(THICKNESS) ? (*mpPropMaster)[THICKNESS] : 1.0;
    mIntegrationWeightMaster = r_geometry_master.IntegrationPoints()[0].Weight() * thickness_master;

    const double thickness_slave = (*mpPropSlave).Has(THICKNESS) ? (*mpPropSlave)[THICKNESS] : 1.0;

    // weight on the master skin
    const double curvature_master = mpProjectionNodeMaster->GetValue(CURVATURE);
    const double rho_surrogate_skin_master = inner_prod(mDistanceMaster, mTrueNormalMaster)/std::abs(inner_prod(mDistanceMaster, mTrueNormalMaster)) 
                                            * norm_2(mDistanceMaster);
    double det_J_surrogate_skin_master = inner_prod(mNormalMaster, mTrueNormalMaster)/(1-curvature_master*rho_surrogate_skin_master);

    if (norm_2(mDistanceMaster) < 1e-12) {
        det_J_surrogate_skin_master = 1.0;
    }

    // weight on the slave skin
    const double curvature_slave = mpProjectionNodeSlave->GetValue(CURVATURE);
    const auto distance_slave_master = mpProjectionNodeSlave->Coordinates()-mpProjectionNodeMaster->Coordinates();

    const double rho_skin_master_skin_slave = inner_prod(distance_slave_master, mTrueNormalSlave)/std::abs(inner_prod(distance_slave_master, mTrueNormalSlave)) 
                                            * norm_2(distance_slave_master);
    double det_J_skin_master_skin_slave = -inner_prod(mTrueNormalMaster, mTrueNormalSlave)/(1-curvature_slave*rho_skin_master_skin_slave);

    if (norm_2(distance_slave_master) < 1e-12) {
        det_J_skin_master_skin_slave = 1.0;
    }

    // weight on surrogate slave
    const double rho_surrogate_skin_slave = inner_prod(mDistanceSlave, mTrueNormalSlave)/std::abs(inner_prod(mDistanceSlave, mTrueNormalSlave)) 
                                                * norm_2(mDistanceSlave);
    double det_J_surrogate_skin_slave = inner_prod(mNormalSlave, mTrueNormalSlave)/(1-curvature_slave*rho_surrogate_skin_slave);

    if (norm_2(mDistanceSlave) < 1e-12) {
        det_J_surrogate_skin_slave = 1.0;
    }

    mIntegrationWeightSlave = r_geometry_master.IntegrationPoints()[0].Weight() * det_J_surrogate_skin_master * det_J_skin_master_skin_slave/det_J_surrogate_skin_slave * thickness_slave;
    

    this->SetValue(INTEGRATION_WEIGHT, mIntegrationWeightMaster); 
}

void SbmContact2DCondition::ComputeGradientTaylorExpansionContribution(
    const GeometryType& rGeometry,
    const Vector& rDistanceVector,
    const SizeType BasisFunctionsOrder,
    Matrix& rGradHSum) const
{
    const SizeType number_of_control_points = rGeometry.PointsNumber();
    const auto& r_DN_De = rGeometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(BasisFunctionsOrder);
    for (IndexType n = 1; n <= BasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    if (rGradHSum.size1() != 3 || rGradHSum.size2() != number_of_control_points)
    {
        rGradHSum.resize(3, number_of_control_points);
    }

    // Neumann (Taylor expansion of the gradient)
    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // Reset for each control point
        double H_taylor_term_X = 0.0; 
        double H_taylor_term_Y = 0.0; 
        double H_taylor_term_Z = 0.0; 

        if (mMasterDim == 2 && mSlaveDim == 2) {
            for (IndexType n = 2; n <= BasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_X += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k+1); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            // 3D
            // TODO:
            KRATOS_ERROR << "3D not implemented yet in SbmContact2DCondition::ComputeGradientTaylorExpansionContribution" << std::endl;
        }
        
        rGradHSum(0,i) = H_taylor_term_X + r_DN_De[0](i, 0);
        rGradHSum(1,i) = H_taylor_term_Y + r_DN_De[0](i, 1);
        // if (mDim == 3)
        //     rGradHSum(2,i) = H_taylor_term_Z + r_DN_De[0](i, 2);
        // else 
            rGradHSum(2,i) = 0;
    }
}

void SbmContact2DCondition::ComputeTaylorSumContribution(
    const GeometryType& rGeometry,
    const Vector& rDistanceVector,
    const SizeType BasisFunctionsOrder,
    Matrix& rHSum) const
{
    const SizeType number_of_nodes = rGeometry.PointsNumber();
    const Matrix& rN = rGeometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    if (rHSum.size1() != 1 || rHSum.size2() != number_of_nodes) {
        rHSum.resize(1, number_of_nodes, false);
    }

    noalias(rHSum) = ZeroMatrix(1, number_of_nodes);

    if (BasisFunctionsOrder == 0) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rHSum(0, i) = rN(0, i);
        }
        return;
    }

    std::vector<Matrix> shape_function_derivatives(BasisFunctionsOrder);
    for (IndexType n = 1; n <= BasisFunctionsOrder; ++n) {
        shape_function_derivatives[n - 1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        double taylor_term = 0.0;
        for (IndexType n = 1; n <= BasisFunctionsOrder; ++n) {
            Matrix& r_shape_function_derivatives = shape_function_derivatives[n - 1];
            for (IndexType k = 0; k <= n; ++k) {
                const IndexType n_k = n - k;
                const double derivative = r_shape_function_derivatives(i, k);
                taylor_term += ComputeTaylorTerm(
                    derivative,
                    rDistanceVector[0],
                    static_cast<int>(n_k),
                    rDistanceVector[1],
                    static_cast<int>(k));
            }
        }

        rHSum(0, i) = rN(0, i) + taylor_term;
    }
}

void SbmContact2DCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType dim = 2;

    const auto& r_geometry_master = GetMasterGeometry();
    const SizeType number_of_nodes_master = r_geometry_master.size();

    const auto& r_geometry_slave = GetSlaveGeometry();
    const SizeType number_of_nodes_slave = r_geometry_slave.size();

    const SizeType mat_size = (number_of_nodes_master+number_of_nodes_slave) * dim;

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}

    void SbmContact2DCondition::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        const int activation_level = this->GetValue(ACTIVATION_LEVEL);

        // INITIALIZE AND RESIZE 

        double penalty_master = (*mpPropMaster)[PENALTY_FACTOR];
        double penalty_slave  = (*mpPropSlave)[PENALTY_FACTOR];

        Vector penalty_vector = ZeroVector(2);
        penalty_vector[0] = penalty_master; penalty_vector[1] = penalty_slave;

        this->SetValue(PENALTY, penalty_vector);
        const SizeType dim = 2;

        const auto& r_geometry_master = GetMasterGeometry();
        const SizeType number_of_nodes_master = r_geometry_master.size();

        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType mat_size = (number_of_nodes_master+number_of_nodes_slave) * dim;
        //resizing as needed the LHS
        if(rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size )
            rLeftHandSideMatrix.resize(mat_size,mat_size,false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS


        // Shape function derivatives 
        // MASTER ------------------------------------------------
        // Initialize Jacobian
        GeometryType::JacobiansType J0_master;
        // Initialize DN_DX
        Matrix DN_DX_master(number_of_nodes_master,2);
        Matrix InvJ0_master(dim,dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De_master = r_geometry_master.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry_master.Jacobian(J0_master,this->GetIntegrationMethod());

        double DetJ0_master;

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0_master[0](0,0); Jacobian(0,1) = J0_master[0](0,1);
        Jacobian(1,0) = J0_master[0](1,0); Jacobian(1,1) = J0_master[0](1,1);

        // // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0_master,DetJ0_master);

        KRATOS_ERROR_IF(std::abs(DetJ0_master -1.0) > 1e-11) << ":::[SbmContact2DCondition]::: Error. Master: Not coincident physical and parameter spaces" 
                        << "Jacobian:   \n" << Jacobian << std::endl;

        const Matrix& r_N_master = r_geometry_master.ShapeFunctionsValues(this->GetIntegrationMethod());

        // // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX_master) = DN_De_master[0];

        // // SLAVE ------------------------------------------------
        // Initialize Jacobian
        GeometryType::JacobiansType J0_slave;
        // Initialize DN_DX
        Matrix DN_DX_slave(number_of_nodes_slave,2);
        Matrix InvJ0_slave(dim,dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De_slave = r_geometry_slave.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry_slave.Jacobian(J0_slave,this->GetIntegrationMethod());
        
        double DetJ0_slave;
        Matrix Jacobian_slave = ZeroMatrix(2,2);
        Jacobian_slave(0,0) = J0_slave[0](0,0); Jacobian_slave(0,1) = J0_slave[0](0,1);
        Jacobian_slave(1,0) = J0_slave[0](1,0); Jacobian_slave(1,1) = J0_slave[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian_slave,InvJ0_slave,DetJ0_slave);  //DetJ0 is not the true one for NURBS geometries

        KRATOS_ERROR_IF(std::abs(DetJ0_slave -1.0) > 1e-11) << ":::[SbmContact2DCondition]::: Error. Slave: Not coincident physical and parameter spaces" 
                        << "detJ" << DetJ0_slave << std::endl;

        const Matrix& r_N_slave = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        // noalias(DN_DX_slave) = prod(DN_De_slave[0],InvJ0_slave);
        noalias(DN_DX_slave) = DN_De_slave[0];

        const SizeType shift_dof = 2*number_of_nodes_master;

        // //------------------- DB MATRICES 

        // // MODIFIED
        Matrix B_master = ZeroMatrix(3,number_of_nodes_master);

        Matrix B_slave = ZeroMatrix(3,number_of_nodes_slave);

        CalculateB(B_master, DN_DX_master, number_of_nodes_master);

        CalculateB(B_slave, DN_DX_slave, number_of_nodes_slave);

        // //---------- GET CONSTITUTIVE MATRICES  
        const Matrix r_D_master = GetValue(CONSTITUTIVE_MATRIX_MASTER);

        const Matrix r_D_slave  = GetValue(CONSTITUTIVE_MATRIX_SLAVE);

        //  ----------------------------------------------------------------------------
        //  ----------------------------------------------------------------------------
        //  Shfited Boundary Method:: Retrieve infos
        //  ----------------------------------------------------------------------------
        //  ----------------------------------------------------------------------------


        KRATOS_ERROR_IF(mMasterTrueDotSurrogateNormal <= 0.0) << ":::[SbmContact2DCondition]::: Error. Master: Negative n_ntilde_master"
                    << " (condition Id " << this->Id() << ")." << std::endl;
        KRATOS_ERROR_IF(mSlaveTrueDotSurrogateNormal <= 0.0) << ":::[SbmContact2DCondition]::: Error. Slave: Negative n_ntilde_slave"
                    << " (condition Id " << this->Id() << ")." << std::endl;


        Matrix H_sum_master;
        Matrix H_sum_slave;
        ComputeTaylorSumContribution(r_geometry_master, mDistanceMaster, 2*mMasterBasisFunctionsOrder, H_sum_master);
        ComputeTaylorSumContribution(r_geometry_slave, mDistanceSlave, 2*mSlaveBasisFunctionsOrder, H_sum_slave);

        Matrix grad_H_sum_master_transposed;
        ComputeGradientTaylorExpansionContribution(r_geometry_master, mDistanceMaster, mMasterBasisFunctionsOrder, grad_H_sum_master_transposed);
        Matrix grad_H_sum_master = trans(grad_H_sum_master_transposed);

        Matrix B_sum_master = ZeroMatrix(3, number_of_nodes_master);
        CalculateB(B_sum_master, grad_H_sum_master, number_of_nodes_master);


        Matrix grad_H_sum_slave_transposed;
        ComputeGradientTaylorExpansionContribution(r_geometry_slave, mDistanceSlave, mSlaveBasisFunctionsOrder, grad_H_sum_slave_transposed);
        Matrix grad_H_sum_slave = trans(grad_H_sum_slave_transposed);

        Matrix B_sum_slave = ZeroMatrix(3, number_of_nodes_slave);
        CalculateB(B_sum_slave, grad_H_sum_slave, number_of_nodes_slave);

        const double h = std::min(mMasterCharacteristicLength, mSlaveCharacteristicLength);

        // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
        penalty_master = penalty_master / h * mMasterBasisFunctionsOrder *mMasterBasisFunctionsOrder;
        penalty_slave = penalty_slave / h * mSlaveBasisFunctionsOrder *mSlaveBasisFunctionsOrder;

        // //-----------
    
        // // Assembly
        Matrix DB_master = prod(r_D_master,B_master);
        Matrix DB_slave = prod(r_D_slave,B_slave);
        Matrix DB_sum_master = prod(r_D_master,B_sum_master);
        Matrix DB_sum_slave = prod(r_D_slave,B_sum_slave);

        // //***************
        // FLUX TERMS  */ both for active and inactive part

        // Flux Term and SBM extension term on Master 
        // -[sigma_1(u) \dot n_tilde] \dot * w_1) -- flux term
        // -[E(sigma_1(u)) \dot n_1] \dot * w_1 * (n1 \dot n_tilde)) -- extension term
        for (IndexType i = 0; i < number_of_nodes_master; i++)
        {
            for (IndexType j = 0; j < number_of_nodes_master; j++)
            {
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int iglob = 2*i+idim;
                    for (IndexType jdim = 0; jdim < 2; jdim++) 
                    {
                        const int jglob = 2*j+jdim;

                        // // FLUX STANDARD T<ERM
                        Vector sigma_u_n(2);
                        sigma_u_n[0] = (DB_master(0, jglob)* mNormalMaster[0] + DB_master(2, jglob)* mNormalMaster[1]);
                        sigma_u_n[1] = (DB_master(2, jglob)* mNormalMaster[0] + DB_master(1, jglob)* mNormalMaster[1]);
                        
                        rLeftHandSideMatrix(iglob, jglob) -= r_N_master(0,i)*sigma_u_n[idim] * mIntegrationWeightMaster;
                        
                        // FLUX EXTENSION TERM
                        Vector extension_sigma_u_n(2);
                        extension_sigma_u_n[0] = (DB_sum_master(0, jglob)* mTrueNormalMaster[0] + DB_sum_master(2, jglob)* mTrueNormalMaster[1]);
                        extension_sigma_u_n[1] = (DB_sum_master(2, jglob)* mTrueNormalMaster[0] + DB_sum_master(1, jglob)* mTrueNormalMaster[1]);

                        rLeftHandSideMatrix(iglob, jglob) += r_N_master(0,i)*extension_sigma_u_n[idim] * mMasterTrueDotSurrogateNormal 
                                                            * mIntegrationWeightMaster;
                    }
                }
            }
        }

        // Flux Term and SBM extension term on Slave 
        // -[sigma_2(u) \dot n_tilde] \dot * w_2) -- flux term
        // -[E(sigma_2(u)) \dot n_2] \dot * w_2 * (n_2 \dot n_tilde)) -- extension term
        for (IndexType i = 0; i < number_of_nodes_slave; i++)
        {
            for (IndexType j = 0; j < number_of_nodes_slave; j++)
            {
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int iglob = 2*i+idim + shift_dof;
                    for (IndexType jdim = 0; jdim < 2; jdim++) 
                    {
                        const int j_index = 2*j+jdim;
                        const int jglob = 2*j+jdim + shift_dof;

                        // // // FLUX STANDARD TERM
                        Vector sigma_u_n(2);
                        sigma_u_n[0] = (DB_slave(0, j_index)* mNormalSlave[0] + DB_slave(2, j_index)* mNormalSlave[1]);
                        sigma_u_n[1] = (DB_slave(2, j_index)* mNormalSlave[0] + DB_slave(1, j_index)* mNormalSlave[1]);
                        
                        rLeftHandSideMatrix(iglob, jglob) -= r_N_slave(0,i)*sigma_u_n[idim] * mIntegrationWeightSlave;
                        
                        // // // FLUX EXTENSION TERM
                        Vector extension_sigma_u_n(2);
                        extension_sigma_u_n[0] = (DB_sum_slave(0, j_index)* mTrueNormalSlave[0] + DB_sum_slave(2, j_index)* mTrueNormalSlave[1]);
                        extension_sigma_u_n[1] = (DB_sum_slave(2, j_index)* mTrueNormalSlave[0] + DB_sum_slave(1, j_index)* mTrueNormalSlave[1]);

                        rLeftHandSideMatrix(iglob, jglob) += r_N_slave(0,i)*extension_sigma_u_n[idim] * mSlaveTrueDotSurrogateNormal 
                                                            * mIntegrationWeightSlave;

                    }
                }
            }
        }

        // ---------------------- ACTIVE PART ------------------------------------------------ //
        // FLUX CONTINUITY
        // -E(sigma_1(u_1))*n_1*n_1 * [(n1*n_1_tilde) * (n_1 * v_1) - (n2*n_2_tilde) * (n_1 * v_2)]
        
        // weight on the master skin
        const double curvature_master = mpProjectionNodeMaster->GetValue(CURVATURE);
        const double rho_surrogate_skin_master = inner_prod(mDistanceMaster, mTrueNormalMaster)/std::abs(inner_prod(mDistanceMaster, mTrueNormalMaster)) 
                                                * norm_2(mDistanceMaster);
        double det_J_surrogate_skin_master = inner_prod(mNormalMaster, mTrueNormalMaster)/(1-curvature_master*rho_surrogate_skin_master);

        if (norm_2(mDistanceMaster) < 1e-12) {
            det_J_surrogate_skin_master = 1.0;
        }

        // weight on the slave skin
        const double curvature_slave = mpProjectionNodeSlave->GetValue(CURVATURE);
        const auto distance_slave_master = mpProjectionNodeSlave->Coordinates()-mpProjectionNodeMaster->Coordinates();

        const double rho_skin_master_skin_slave = inner_prod(distance_slave_master, mTrueNormalSlave)/std::abs(inner_prod(distance_slave_master, mTrueNormalSlave)) 
                                                * norm_2(distance_slave_master);
        double det_J_skin_master_skin_slave = -inner_prod(mTrueNormalMaster, mTrueNormalSlave)/(1-curvature_slave*rho_skin_master_skin_slave);

        if (norm_2(distance_slave_master) < 1e-12) {
            det_J_skin_master_skin_slave = 1.0;
        }

        // weight on surrogate slave
        const double rho_surrogate_skin_slave = inner_prod(mDistanceSlave, mTrueNormalSlave)/std::abs(inner_prod(mDistanceSlave, mTrueNormalSlave)) 
                                                    * norm_2(mDistanceSlave);
        double det_J_surrogate_skin_slave = inner_prod(mNormalSlave, mTrueNormalSlave)/(1-curvature_slave*rho_surrogate_skin_slave);

        if (norm_2(mDistanceSlave) < 1e-12) {
            det_J_surrogate_skin_slave = 1.0;
        }
        
        if (activation_level == 1 || activation_level == 3)
        {
            // -E(sigma_1(u_1))*n_1*n_1 * (n1*n_1_tilde) * (n_1 * v_1)
            for (IndexType i = 0; i < number_of_nodes_master; i++)
            {
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int iglob = 2*i+idim;
                    for (IndexType j = 0; j < number_of_nodes_master; j++)
                    {
                        for (IndexType jdim = 0; jdim < 2; jdim++) 
                        {
                            const int jglob = 2*j+jdim;
    
                            // // FLUX EXTENSION TERM
                            Vector extension_sigma_u_n(2);
                            extension_sigma_u_n[0] = (DB_sum_master(0, jglob)* mTrueNormalMaster[0] + DB_sum_master(2, jglob)* mTrueNormalMaster[1]);
                            extension_sigma_u_n[1] = (DB_sum_master(2, jglob)* mTrueNormalMaster[0] + DB_sum_master(1, jglob)* mTrueNormalMaster[1]);

                            const double extension_sigma_u_n_dot_n = extension_sigma_u_n[0]*mTrueNormalMaster[0] + extension_sigma_u_n[1]*mTrueNormalMaster[1];
    
                            rLeftHandSideMatrix(iglob, jglob) -= 0.5*extension_sigma_u_n_dot_n * mMasterTrueDotSurrogateNormal * mTrueNormalMaster[idim]
                                                                * r_N_master(0,i) * mIntegrationWeightMaster;
                            
                        }
                    }

                    for (IndexType j = 0; j < number_of_nodes_slave; j++)
                    {
                        for (IndexType jdim = 0; jdim < 2; jdim++) 
                        {
                            const int jglob = 2*j+jdim + shift_dof;
                            const int j_index =  2*j+jdim;
    
                            // // FLUX EXTENSION TERM
                            Vector extension_sigma_u_n(2);
                            extension_sigma_u_n[0] = (DB_sum_slave(0, j_index)* mTrueNormalSlave[0] + DB_sum_slave(2, j_index)* mTrueNormalSlave[1]);
                            extension_sigma_u_n[1] = (DB_sum_slave(2, j_index)* mTrueNormalSlave[0] + DB_sum_slave(1, j_index)* mTrueNormalSlave[1]);

                            const double extension_sigma_u_n_dot_n = extension_sigma_u_n[0]*mTrueNormalMaster[0] + extension_sigma_u_n[1]*mTrueNormalMaster[1];
    
                            rLeftHandSideMatrix(iglob, jglob) += 0.5*extension_sigma_u_n_dot_n * mMasterTrueDotSurrogateNormal * mTrueNormalMaster[idim]
                                                                * r_N_master(0,i) * mIntegrationWeightMaster;
                            
                        }
                    }
                }
            }

            // +E(sigma_1(u_1))*n_1*n_1 * (n2*n_2_tilde) * (n_1 * v_2)
            for (IndexType i = 0; i < number_of_nodes_slave; i++)
            {
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int iglob = 2*i+idim + shift_dof;
                    for (IndexType j = 0; j < number_of_nodes_master; j++)
                    {
                        for (IndexType jdim = 0; jdim < 2; jdim++) 
                        {
                            const int jglob = 2*j+jdim;
    
                            // // FLUX EXTENSION TERM
                            Vector extension_sigma_u_n(2);
                            extension_sigma_u_n[0] = (DB_sum_master(0, jglob)* mTrueNormalMaster[0] + DB_sum_master(2, jglob)* mTrueNormalMaster[1]);
                            extension_sigma_u_n[1] = (DB_sum_master(2, jglob)* mTrueNormalMaster[0] + DB_sum_master(1, jglob)* mTrueNormalMaster[1]);

                            const double extension_sigma_u_n_dot_n = extension_sigma_u_n[0]*mTrueNormalMaster[0] + extension_sigma_u_n[1]*mTrueNormalMaster[1];
    
                            rLeftHandSideMatrix(iglob, jglob) += 0.5*extension_sigma_u_n_dot_n * mSlaveTrueDotSurrogateNormal * mTrueNormalMaster[idim]
                                                                * r_N_slave(0,i) * mIntegrationWeightMaster *det_J_surrogate_skin_master/det_J_surrogate_skin_slave;
                            
                            // FIXME: just a trial
                            // rLeftHandSideMatrix(iglob, jglob) += 0.5*extension_sigma_u_n_dot_n * mSlaveTrueDotSurrogateNormal * mTrueNormalMaster[idim]
                            //                                     * r_N_slave(0,i) * mIntegrationWeightSlave;;
                            
                        }
                    }

                    for (IndexType j = 0; j < number_of_nodes_slave; j++)
                    {
                        for (IndexType jdim = 0; jdim < 2; jdim++) 
                        {
                            const int jglob = 2*j+jdim + shift_dof;
                            const int j_index =  2*j+jdim;
    
                            // // FLUX EXTENSION TERM
                            Vector extension_sigma_u_n(2);
                            extension_sigma_u_n[0] = (DB_sum_slave(0, j_index)* mTrueNormalSlave[0] + DB_sum_slave(2, j_index)* mTrueNormalSlave[1]);
                            extension_sigma_u_n[1] = (DB_sum_slave(2, j_index)* mTrueNormalSlave[0] + DB_sum_slave(1, j_index)* mTrueNormalSlave[1]);

                            const double extension_sigma_u_n_dot_n = extension_sigma_u_n[0]*mTrueNormalMaster[0] + extension_sigma_u_n[1]*mTrueNormalMaster[1];
    
                            rLeftHandSideMatrix(iglob, jglob) -= 0.5*extension_sigma_u_n_dot_n * mSlaveTrueDotSurrogateNormal * mTrueNormalMaster[idim]
                                                                * r_N_slave(0,i) * mIntegrationWeightMaster *det_J_surrogate_skin_master/det_J_surrogate_skin_slave;
                            
                            // FIXME: just a trial
                            // rLeftHandSideMatrix(iglob, jglob) -= 0.5*extension_sigma_u_n_dot_n * mSlaveTrueDotSurrogateNormal * mTrueNormalMaster[idim]
                            //                                     * r_N_slave(0,i) * mIntegrationWeightSlave;
                      
                        }
                    }
                }
            }
            // DISPLACEMENT CONTINUITY (penalty only)
            // penalty*(u_1 - u_2)(v_1 - v_2)

            for (IndexType i = 0; i < number_of_nodes_master; i++)
            {
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int iglob = 2*i+idim;
                    for (IndexType j = 0; j < number_of_nodes_master; j++)
                    {
                        for (IndexType jdim = 0; jdim < 2; jdim++) 
                        {
                            const int jglob = 2*j+jdim;
                            rLeftHandSideMatrix(iglob, jglob) += penalty_master * r_N_master(0,i) * mTrueNormalMaster[idim]
                                                                * H_sum_master(0,j) * mTrueNormalMaster[jdim] *mIntegrationWeightMaster;             
                        }
                    }

                    for (IndexType j = 0; j < number_of_nodes_slave; j++)
                    {
                        for (IndexType jdim = 0; jdim < 2; jdim++) 
                        {
                            const int jglob = 2*j+jdim + shift_dof;
                            rLeftHandSideMatrix(iglob, jglob) -= penalty_master * r_N_master(0,i) * mTrueNormalMaster[idim]
                                                                * H_sum_slave(0,j) * mTrueNormalMaster[jdim] *mIntegrationWeightMaster;             
                        }
                    }
                }
            }
            for (IndexType i = 0; i < number_of_nodes_slave; i++)
            {
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int iglob = 2*i+idim + shift_dof;
                    for (IndexType j = 0; j < number_of_nodes_master; j++)
                    {
                        for (IndexType jdim = 0; jdim < 2; jdim++) 
                        {
                            const int jglob = 2*j+jdim;
                            rLeftHandSideMatrix(iglob, jglob) -= penalty_slave * r_N_slave(0,i) * mTrueNormalMaster[idim]
                                                                * H_sum_master(0,j) * mTrueNormalMaster[jdim] *mIntegrationWeightSlave;      
                        }
                    }

                    for (IndexType j = 0; j < number_of_nodes_slave; j++)
                    {
                        for (IndexType jdim = 0; jdim < 2; jdim++) 
                        {
                            const int jglob = 2*j+jdim + shift_dof;
                            rLeftHandSideMatrix(iglob, jglob) += penalty_slave * r_N_slave(0,i) * mTrueNormalMaster[idim]
                                                                * H_sum_slave(0,j) * mTrueNormalMaster[jdim] *mIntegrationWeightSlave;             
                        }
                    }
                }
            }
        }
        

        for (unsigned int i = 0; i < number_of_nodes_slave; i++) {

            std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
            outputFile << r_geometry_slave[i].GetId() << "  " <<r_geometry_slave[i].GetDof(DISPLACEMENT_X).EquationId() <<"\n";
            
            outputFile.close();
        }

        for (unsigned int i = 0; i < number_of_nodes_master; i++) {

            std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
            outputFile << r_geometry_master[i].GetId() << "  " <<r_geometry_master[i].GetDof(DISPLACEMENT_X).EquationId() <<"\n";

            outputFile.close();
        }

    //     // /////////////////////////////////////////////////////////////////////////////////////////
     
        KRATOS_CATCH("")
    }

    void SbmContact2DCondition::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY
    
        const int activation_level = this->GetValue(ACTIVATION_LEVEL);
        const SizeType dim = 2;
    
        const auto& r_geometry_master = GetMasterGeometry();
        const SizeType number_of_nodes_master = r_geometry_master.size();
    
        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();
    
        const SizeType mat_size = (number_of_nodes_master + number_of_nodes_slave) * dim;
    
        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size, false);
        noalias(rRightHandSideVector) = ZeroVector(mat_size);


        Vector check_right_side_vector = ZeroVector(mat_size);
    
        const Matrix& r_N_master = r_geometry_master.ShapeFunctionsValues(this->GetIntegrationMethod());
        const Matrix& r_N_slave = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());

        double penalty_master = (*mpPropMaster)[PENALTY_FACTOR];
        double penalty_slave =  (*mpPropSlave)[PENALTY_FACTOR];

        const double h = std::min(mMasterCharacteristicLength, mSlaveCharacteristicLength);
        penalty_master = penalty_master / h * mMasterBasisFunctionsOrder *mMasterBasisFunctionsOrder;
        penalty_slave = penalty_slave / h * mSlaveBasisFunctionsOrder *mSlaveBasisFunctionsOrder;

        const SizeType shift_dof = dim * number_of_nodes_master;

        // compute B and DB
        Matrix DN_DX_master(number_of_nodes_master,2);

        const GeometryType::ShapeFunctionsGradientsType& DN_De_master = r_geometry_master.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

        // // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX_master) = DN_De_master[0];

        Matrix DN_DX_slave(number_of_nodes_slave,2);
        Matrix InvJ0_slave(dim,dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De_slave = r_geometry_slave.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

        noalias(DN_DX_slave) = DN_De_slave[0];

        // // MODIFIED
        Matrix B_master = ZeroMatrix(3,number_of_nodes_master);

        Matrix B_slave = ZeroMatrix(3,number_of_nodes_slave);

        CalculateB(B_master, DN_DX_master, number_of_nodes_master);

        CalculateB(B_slave, DN_DX_slave, number_of_nodes_slave);

        // //---------- GET CONSTITUTIVE MATRICES  

        // compute the residual tractions
    
        Vector traction_master = ZeroVector(3);
        Vector traction_master_true = ZeroVector(3);
        Vector traction_slave = ZeroVector(3);
        Vector traction_slave_true = ZeroVector(3);
        
        // get stress vector on surrogate master
        ConstitutiveLaw::Parameters values_surrogate_master(r_geometry_master, GetPropertiesMaster(), rCurrentProcessInfo);
        Vector old_displacement_coefficient_vector_master(2 * number_of_nodes_master);
        GetValuesVector(old_displacement_coefficient_vector_master, QuadraturePointCouplingGeometry2D<Point>::Master);
        Vector old_strain_surrogate_master = prod(B_master, old_displacement_coefficient_vector_master);

        const SizeType strain_size_surrogate_master = mpConstitutiveLawMaster->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_surrogate_master(strain_size_surrogate_master);
        ApplyConstitutiveLaw(mpConstitutiveLawMaster, old_strain_surrogate_master, values_surrogate_master, 
                            this_constitutive_variables_surrogate_master);

        const Vector& r_stress_vector_surrogate_master = values_surrogate_master.GetStressVector();

        // get stress vector on surrogate slave
        ConstitutiveLaw::Parameters values_surrogate_slave(r_geometry_slave, GetPropertiesSlave(), rCurrentProcessInfo);
        Vector old_displacement_coefficient_vector_slave(2 * number_of_nodes_slave);
        GetValuesVector(old_displacement_coefficient_vector_slave, QuadraturePointCouplingGeometry2D<Point>::Slave);
        Vector old_strain_surrogate_slave = prod(B_slave, old_displacement_coefficient_vector_slave);

        const SizeType strain_size_surrogate_slave = mpConstitutiveLawSlave->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_surrogate_slave(strain_size_surrogate_slave);
        ApplyConstitutiveLaw(mpConstitutiveLawSlave, old_strain_surrogate_slave, values_surrogate_slave, 
                            this_constitutive_variables_surrogate_slave);

        const Vector& r_stress_vector_surrogate_slave = values_surrogate_slave.GetStressVector();

        // get stress vector on skin slave ---------------------------------------------------------------------------------------
        Matrix grad_H_sum_slave_transposed;
        ComputeGradientTaylorExpansionContribution(r_geometry_slave, mDistanceSlave, mSlaveBasisFunctionsOrder, grad_H_sum_slave_transposed);
        Matrix grad_H_sum_slave = trans(grad_H_sum_slave_transposed);

        Matrix B_sum_slave = ZeroMatrix(3, number_of_nodes_slave);
        CalculateB(B_sum_slave, grad_H_sum_slave, number_of_nodes_slave);

        ConstitutiveLaw::Parameters values_skin_slave(r_geometry_slave, GetPropertiesSlave(), rCurrentProcessInfo);
        Vector old_strain_skin_slave = prod(B_sum_slave, old_displacement_coefficient_vector_slave);

        const SizeType strain_size_skin_slave = mpConstitutiveLawSlave->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_skin_slave(strain_size_skin_slave);
        ApplyConstitutiveLaw(mpConstitutiveLawSlave, old_strain_skin_slave, values_skin_slave, 
                            this_constitutive_variables_skin_slave);

        Vector& r_stress_slave = values_skin_slave.GetStressVector();
    
        traction_master[0] = r_stress_vector_surrogate_master[0] * mNormalMaster[0] + r_stress_vector_surrogate_master[2] * mNormalMaster[1];
        traction_master[1] = r_stress_vector_surrogate_master[2] * mNormalMaster[0] + r_stress_vector_surrogate_master[1] * mNormalMaster[1];

        traction_slave[0] = r_stress_vector_surrogate_slave[0] * mNormalSlave[0] + r_stress_vector_surrogate_slave[2] * mNormalSlave[1];
        traction_slave[1] = r_stress_vector_surrogate_slave[2] * mNormalSlave[0] + r_stress_vector_surrogate_slave[1] * mNormalSlave[1];

        if (this->Has(STRESS_MASTER)) {
            const Vector& r_stress_master = this->GetValue(STRESS_MASTER);
            if (r_stress_master.size() >= 3) {
                traction_master_true[0] = r_stress_master[0] * mTrueNormalMaster[0] + r_stress_master[2] * mTrueNormalMaster[1];
                traction_master_true[1] = r_stress_master[2] * mTrueNormalMaster[0] + r_stress_master[1] * mTrueNormalMaster[1];
            }
        }

        if (this->Has(STRESS_SLAVE)) {

            if (r_stress_slave.size() >= 3) {
                traction_slave_true[0] = r_stress_slave[0] * mTrueNormalSlave[0] + r_stress_slave[2] * mTrueNormalSlave[1];
                traction_slave_true[1] = r_stress_slave[2] * mTrueNormalSlave[0] + r_stress_slave[1] * mTrueNormalSlave[1];
            }
        }
    
        const double traction_master_true_normal = inner_prod(traction_master_true, mTrueNormalMaster);
        const double traction_slave_true_normal = inner_prod(traction_slave_true, mTrueNormalMaster);
    
        double gap_normal_master = 0.0;
        if (this->Has(GAP)) {
            const Vector& r_gap = -this->GetValue(GAP);
            if (r_gap.size() >= 2) {
                gap_normal_master = r_gap[0] * mTrueNormalMaster[0] + r_gap[1] * mTrueNormalMaster[1];
            }
        }

        // weight on the master skin
        const double curvature_master = mpProjectionNodeMaster->GetValue(CURVATURE);
        const double rho_surrogate_skin_master = inner_prod(mDistanceMaster, mTrueNormalMaster)/std::abs(inner_prod(mDistanceMaster, mTrueNormalMaster)) 
                                                * norm_2(mDistanceMaster);
        double det_J_surrogate_skin_master = inner_prod(mNormalMaster, mTrueNormalMaster)/(1-curvature_master*rho_surrogate_skin_master);

        if (norm_2(mDistanceMaster) < 1e-12) {
            det_J_surrogate_skin_master = 1.0;
        }

        // weight on the slave skin
        const double curvature_slave = mpProjectionNodeSlave->GetValue(CURVATURE);
        const auto distance_slave_master = mpProjectionNodeSlave->Coordinates()-mpProjectionNodeMaster->Coordinates();

        const double rho_skin_master_skin_slave = inner_prod(distance_slave_master, mTrueNormalSlave)/std::abs(inner_prod(distance_slave_master, mTrueNormalSlave)) 
                                                * norm_2(distance_slave_master);
        double det_J_skin_master_skin_slave = -inner_prod(mTrueNormalMaster, mTrueNormalSlave)/(1-curvature_slave*rho_skin_master_skin_slave);

        if (norm_2(distance_slave_master) < 1e-12) {
            det_J_skin_master_skin_slave = 1.0;
        }

        // weight on surrogate slave
        const double rho_surrogate_skin_slave = inner_prod(mDistanceSlave, mTrueNormalSlave)/std::abs(inner_prod(mDistanceSlave, mTrueNormalSlave)) 
                                                    * norm_2(mDistanceSlave);
        double det_J_surrogate_skin_slave = inner_prod(mNormalSlave, mTrueNormalSlave)/(1-curvature_slave*rho_surrogate_skin_slave);

        if (norm_2(mDistanceSlave) < 1e-12) {
            det_J_surrogate_skin_slave = 1.0;
        }


        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            for (IndexType idim = 0; idim < dim; ++idim) {
                const IndexType iglob = dim * i + idim;
                double contribution = 0.0;
    
                contribution += r_N_master(0, i) * traction_master[idim] * mIntegrationWeightMaster;
                contribution -= r_N_master(0, i) * traction_master_true[idim] * mMasterTrueDotSurrogateNormal * mIntegrationWeightMaster;

                check_right_side_vector[iglob] += contribution;
    
                if (activation_level == 1 || activation_level == 3) {
                    contribution += r_N_master(0, i) * 0.5 *
                                    (traction_master_true_normal - traction_slave_true_normal)
                                    * mTrueNormalMaster[idim] * mMasterTrueDotSurrogateNormal * mIntegrationWeightMaster;
                
                    contribution -= penalty_master * r_N_master(0, i) * mTrueNormalMaster[idim]
                        * gap_normal_master * mIntegrationWeightMaster;

                }
    
                rRightHandSideVector[iglob] += contribution;
            }
        }
        
        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            for (IndexType idim = 0; idim < dim; ++idim) {
                const IndexType iglob = shift_dof + dim * i + idim;
                double contribution = 0.0;
                
                contribution += r_N_slave(0, i) * traction_slave[idim] * mIntegrationWeightSlave;
                contribution -= r_N_slave(0, i) * traction_slave_true[idim] * mSlaveTrueDotSurrogateNormal * mIntegrationWeightSlave;

                check_right_side_vector[iglob] += contribution;
    
                if (activation_level == 1 || activation_level == 3) {
                    contribution -= r_N_slave(0, i) * 0.5 *
                                    (traction_master_true_normal - traction_slave_true_normal) * det_J_surrogate_skin_master/det_J_surrogate_skin_slave
                                    * mTrueNormalMaster[idim] * mSlaveTrueDotSurrogateNormal * mIntegrationWeightMaster;

                    // FIXME: just a trial
                    // contribution -= r_N_slave(0, i) * 0.5 *
                    //                 (traction_master_true_normal - traction_slave_true_normal) 
                    //                 * mTrueNormalMaster[idim] * mSlaveTrueDotSurrogateNormal * mIntegrationWeightSlave;
                                    
                    contribution += penalty_slave * r_N_slave(0, i) * mTrueNormalMaster[idim]
                        * gap_normal_master * mIntegrationWeightSlave;

                }
    
                rRightHandSideVector[iglob] += contribution;
            }
        }

        // Matrix rLeftHandSideMatrix;
        // noalias(rRightHandSideVector) = ZeroVector(mat_size);
        // CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
        // Vector solution_coefficient_vector; //(mat_size);
        // GetValuesVector(solution_coefficient_vector, 2);
        // rRightHandSideVector -= prod(rLeftHandSideMatrix, solution_coefficient_vector);
        

        KRATOS_CATCH("")
    }

    int SbmContact2DCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT((*mpPropMaster).Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void SbmContact2DCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry_master = GetMasterGeometry();
        const SizeType number_of_nodes_master = r_geometry_master.size();

        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType number_of_nodes = (number_of_nodes_master + number_of_nodes_slave);

        if (rResult.size() != 2 * number_of_nodes)
            rResult.resize(2 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry_master[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }

        const int shift = number_of_nodes_master*2;
        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry_slave[i];
            rResult[index+shift] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1+shift] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

    void SbmContact2DCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry_master = GetMasterGeometry();
        const SizeType number_of_nodes_master = r_geometry_master.size();

        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * (number_of_nodes_master+number_of_nodes_slave));

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }

    };


    void SbmContact2DCondition::GetValuesVector(
        Vector& rValues, IndexType index) const
    {
        if (index == QuadraturePointCouplingGeometry2D<Point>::Master) {
            
            SizeType number_of_control_points = GetMasterGeometry().size();
            const SizeType mat_size = number_of_control_points * 2;

            if (rValues.size() != mat_size)
                rValues.resize(mat_size, false);

            for (IndexType i = 0; i < number_of_control_points; ++i)
            {
                const array_1d<double, 3>& displacement = GetMasterGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                IndexType index = i * 2;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }
        } else if (index == QuadraturePointCouplingGeometry2D<Point>::Slave) {
            SizeType number_of_control_points = GetSlaveGeometry().size();
            const SizeType mat_size = number_of_control_points * 2;

            if (rValues.size() != mat_size)
                rValues.resize(mat_size, false);

            for (IndexType i = 0; i < number_of_control_points; ++i)
            {
                const array_1d<double, 3>& displacement = GetSlaveGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                IndexType index = i * 2;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }
        } else if (index == 2) {
            
            SizeType number_of_control_points_master = GetMasterGeometry().size();
            const SizeType mat_size_master = number_of_control_points_master * 2;
            SizeType number_of_control_points_slave = GetSlaveGeometry().size();
            const SizeType mat_size_slave = number_of_control_points_slave * 2;
            const SizeType mat_size = mat_size_master + mat_size_slave;

            if (rValues.size() != mat_size)
                rValues.resize(mat_size, false);

            for (IndexType i = 0; i < number_of_control_points_master; ++i)
            {
                const array_1d<double, 3>& displacement = GetMasterGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                IndexType index = i * 2;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }

            for (IndexType i = 0; i < number_of_control_points_slave; ++i)
            {
                const array_1d<double, 3>& displacement = GetSlaveGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                IndexType index = i * 2 + mat_size_master;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }
        }
    }

    void SbmContact2DCondition::GetStrainVector(
        Vector& rStrainVector, IndexType index) const
    {
        const auto& r_geometry = GetGeometry().GetGeometryPart(index);

        const SizeType dim = 2;//r_geometry.WorkingSpaceDimension();

        const SizeType number_of_nodes = r_geometry.size();


        // Shape function derivatives 
        // Initialize Jacobian
        GeometryType::JacobiansType J0;

        // Initialize DN_DX
        Matrix DN_DX(number_of_nodes,2);
        Matrix InvJ0(dim,dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        double DetJ0;
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);


        // retrieve shape function values
        Vector displacements(2*number_of_nodes);
        GetValuesVector(displacements, index);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = DN_De[0];

        // MODIFIED
        Matrix B = ZeroMatrix(3,number_of_nodes);

        CalculateB(B, DN_DX, number_of_nodes);

        if (rStrainVector.size() != dim) rStrainVector.resize(dim);

        rStrainVector = prod(B, displacements);
    }

    /**
     * @brief 
     * 
     * @param rStrainVector 
     * @param index 
     * @param rCurrentProcessInfo 
     */
    void SbmContact2DCondition::SetConstitutiveVariables(
        Vector& rStrainVector, IndexType index, const Kratos::ProcessInfo& rCurrentProcessInfo) 
    {
        const auto& r_geometry = GetGeometry().GetGeometryPart(index);
        const SizeType dim = 2;//r_geometry.WorkingSpaceDimension();

        ConstitutiveLaw::Pointer rpConstitutiveLaw = GetConstitutiveLaw(index);

        PropertiesType r_prop = GetProperty(index);

        ConstitutiveLaw::Parameters Values(r_geometry, r_prop, rCurrentProcessInfo);

        const SizeType strain_size = rpConstitutiveLaw->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);
    
        Values.SetStrainVector(rStrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);

        Values.SetConstitutiveMatrix(this_constitutive_variables.D);
        rpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

        if (index == 0) {
            this->SetValue(CONSTITUTIVE_MATRIX_MASTER, Values.GetConstitutiveMatrix());
            this->SetValue(STRESS_MASTER, Values.GetStressVector());
        } 
        else if (index == 1)
        {
            this->SetValue(CONSTITUTIVE_MATRIX_SLAVE, Values.GetConstitutiveMatrix());
            this->SetValue(STRESS_SLAVE, Values.GetStressVector());
        }
        
    }

    void SbmContact2DCondition::SetGap() 
    {
        NodeType projection_node_master;         NodeType projection_node_slave;
        //
        projection_node_master = GetMasterGeometry().GetValue(NEIGHBOUR_NODES)[0];
        projection_node_slave = projection_node_master.GetValue(NEIGHBOUR_NODES)[0];

        Vector skin_master_coord_deformed = projection_node_master + GetValue(DISPLACEMENT_MASTER);
        Vector skin_slave_coord_deformed  = projection_node_slave + GetValue(DISPLACEMENT_SLAVE);
        
        SetValue(SKIN_MASTER_COORDINATES, projection_node_master);
        SetValue(SKIN_SLAVE_COORDINATES, projection_node_slave);

        Vector gap_deformed = skin_slave_coord_deformed - skin_master_coord_deformed;

        SetValue(GAP, gap_deformed);
    }

    void SbmContact2DCondition::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX,
        const SizeType number_of_control_points) const
    {
        
        const SizeType mat_size = number_of_control_points * 2;

        if (rB.size1() != 3 || rB.size2() != mat_size)
            rB.resize(3, mat_size);
        noalias(rB) = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 2;
            IndexType dirr = r % 2;

            rB(0, r) = r_DN_DX(kr,0) * (1-dirr);
            rB(1, r) = r_DN_DX(kr,1) * dirr;
            rB(2, r) = r_DN_DX(kr,0) * (dirr) + r_DN_DX(kr,1) * (1-dirr);
        }
    }


    void SbmContact2DCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters_master(
            GetMasterGeometry(), *mpPropMaster, rCurrentProcessInfo);
        
        // Master
        const auto& r_geometry_master = GetMasterGeometry();
        const SizeType number_of_nodes_master = r_geometry_master.size();

        // Slave
        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();


        NodeType projection_node_master;         NodeType projection_node_slave;
        //
        projection_node_master = r_geometry_master.GetValue(NEIGHBOUR_NODES)[0];
        projection_node_slave = projection_node_master.GetValue(NEIGHBOUR_NODES)[0];
        const auto& slave_neighbour_geometries = projection_node_slave.GetValue(NEIGHBOUR_GEOMETRIES);
        KRATOS_ERROR_IF(slave_neighbour_geometries.empty()) << ":::[SbmContact2DCondition]::: Missing neighbour geometries on projection node "
                                                            << projection_node_slave.Id() << " (condition Id " << this->Id() << ")." << std::endl;
        slave_neighbour_geometries[0]->SetValue(IS_ASSEMBLED, false);

        mpConstitutiveLawMaster->FinalizeMaterialResponse(constitutive_law_parameters_master, ConstitutiveLaw::StressMeasure_Cauchy);

        ConstitutiveLaw::Parameters constitutive_law_parameters_slave(
            GetSlaveGeometry(), *mpPropSlave, rCurrentProcessInfo);

        mpConstitutiveLawSlave->FinalizeMaterialResponse(constitutive_law_parameters_slave, ConstitutiveLaw::StressMeasure_Cauchy);

        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        Vector normal_true_master = GetValue(NORMAL_SKIN_MASTER);
        Vector stress_vector_master_on_true = GetValue(STRESS_MASTER);

        Vector sigma_n(2);

        sigma_n[0] = stress_vector_master_on_true[0]*normal_true_master[0] + stress_vector_master_on_true[2]*normal_true_master[1];
        sigma_n[1] = stress_vector_master_on_true[2]*normal_true_master[0] + stress_vector_master_on_true[1]*normal_true_master[1];

        double sigma_n_n = sigma_n[0]*normal_true_master[0] + sigma_n[1]*normal_true_master[1];


        Vector normal_true_slave = GetValue(NORMAL_SKIN_SLAVE);
        Vector stress_vector_slave_on_true = GetValue(STRESS_SLAVE);

        Vector sigma_n_slave(2);

        sigma_n_slave[0] = stress_vector_slave_on_true[0]*normal_true_slave[0] + stress_vector_slave_on_true[2]*normal_true_slave[1];
        sigma_n_slave[1] = stress_vector_slave_on_true[2]*normal_true_slave[0] + stress_vector_slave_on_true[1]*normal_true_slave[1];

        double sigma_n_n_slave = sigma_n_slave[0]*normal_true_slave[0] + sigma_n_slave[1]*normal_true_slave[1];
        
        // if (this->Has(ACTIVATION_LEVEL)) {
        //     int activation_level = this->GetValue(ACTIVATION_LEVEL);
        //     if (activation_level == 1 || activation_level == 3) {
        //         KRATOS_WATCH(stress_vector_slave_on_true)
        //         KRATOS_WATCH(stress_vector_master_on_true)
        //         KRATOS_WATCH(sigma_n_n)
        //         KRATOS_WATCH(sigma_n_n_slave)
        //         KRATOS_WATCH(sigma_n_n+sigma_n_n_slave)

        //         KRATOS_WATCH(r_geometry_master.Center())

        //         KRATOS_WATCH("----------------")
        //     }
        // }
        
        // KRATOS_WATCH("----------")

        SetValue(NORMAL_STRESS, sigma_n);



        // //---------------------
        // // Set the stress vector on the true boundary
        // //---------------------
        const auto& master_center = r_geometry_master.Center();
        const auto& slave_center = r_geometry_slave.Center();

        const std::vector<double> integration_weight_list_master = GetValue(INTEGRATION_WEIGHTS_MASTER);
        const std::vector<Vector> integration_point_list_master = GetValue(INTEGRATION_POINTS_MASTER);
        const SizeType number_of_integration_points_on_true_master = integration_weight_list_master.size();

        const std::vector<double> integration_weight_list_slave = GetValue(INTEGRATION_WEIGHTS_SLAVE);
        const std::vector<Vector> integration_point_list_slave = GetValue(INTEGRATION_POINTS_SLAVE);
        const SizeType number_of_integration_points_on_true_slave = integration_weight_list_slave.size();

        Matrix values_on_true_boundary_master;
        if (number_of_integration_points_on_true_master > 0) {
            values_on_true_boundary_master = ZeroMatrix(number_of_integration_points_on_true_master, 6);

            Vector coefficient_master(mMasterDim * number_of_nodes_master);
            GetValuesVector(coefficient_master, QuadraturePointCouplingGeometry2D<Point>::Master);

            for (IndexType i = 0; i < number_of_integration_points_on_true_master; ++i) {
                const Vector& r_integration_point = integration_point_list_master[i];
                const double integration_weight = integration_weight_list_master[i];

                Vector distance_vector_master = ZeroVector(mMasterDim);
                distance_vector_master[0] = r_integration_point[0] - master_center.X();
                distance_vector_master[1] = r_integration_point[1] - master_center.Y();

                Matrix grad_H_sum_master_transposed;
                ComputeGradientTaylorExpansionContribution(r_geometry_master, distance_vector_master, mMasterBasisFunctionsOrder, grad_H_sum_master_transposed);
                Matrix grad_H_sum_master = trans(grad_H_sum_master_transposed);

                Matrix B_master;
                CalculateB(B_master, grad_H_sum_master, number_of_nodes_master);

                Vector strain_master = prod(B_master, coefficient_master);

                ConstitutiveLaw::Parameters values_master_true(r_geometry_master, GetPropertiesMaster(), rCurrentProcessInfo);
                ConstitutiveVariables constitutive_variables_master(mpConstitutiveLawMaster->GetStrainSize());
                ApplyConstitutiveLaw(mpConstitutiveLawMaster, strain_master, values_master_true, constitutive_variables_master);

                const Vector& stress_vector_master_on_true = values_master_true.GetStressVector();

                Vector normal_stress_master = ZeroVector(3);
                normal_stress_master[0] = stress_vector_master_on_true[0] * normal_true_master[0] + stress_vector_master_on_true[2] * normal_true_master[1];
                normal_stress_master[1] = stress_vector_master_on_true[2] * normal_true_master[0] + stress_vector_master_on_true[1] * normal_true_master[1];

                double normal_stress = normal_stress_master[0] * normal_true_master[0] + normal_stress_master[1] * normal_true_master[1];
                const double shear_stress = -normal_stress_master[0] * normal_true_master[1] + normal_stress_master[1] * normal_true_master[0];
                
                // FIXME: 
                auto gap = GetValue(GAP);
                double normal_gap = gap[0] * normal_true_master[0] + gap[1] * normal_true_master[1];
                auto penalty_master = (*mpPropMaster)[PENALTY_FACTOR];
                const double h = std::min(mMasterCharacteristicLength, mSlaveCharacteristicLength);

                // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
                // penalty_master = penalty_master / h * mMasterBasisFunctionsOrder *mMasterBasisFunctionsOrder;
                // if (normal_gap < 0.0) {
                //     normal_stress = penalty_master * normal_gap;
                // } else {
                //     normal_stress = 0;
                // }

                values_on_true_boundary_master(i, 0) = integration_weight;
                values_on_true_boundary_master(i, 1) = r_integration_point.size() > 0 ? r_integration_point[0] : 0.0;
                values_on_true_boundary_master(i, 2) = r_integration_point.size() > 1 ? r_integration_point[1] : 0.0;
                values_on_true_boundary_master(i, 3) = r_integration_point.size() > 2 ? r_integration_point[2] : 0.0;
                values_on_true_boundary_master(i, 4) = normal_stress;
                values_on_true_boundary_master(i, 5) = shear_stress;
            }
        }

        Matrix values_on_true_boundary_slave;
        if (number_of_integration_points_on_true_slave > 0) {
            values_on_true_boundary_slave = ZeroMatrix(number_of_integration_points_on_true_slave, 6);

            Vector coefficient_slave(mSlaveDim * number_of_nodes_slave);
            GetValuesVector(coefficient_slave, QuadraturePointCouplingGeometry2D<Point>::Slave);

            for (IndexType i = 0; i < number_of_integration_points_on_true_slave; ++i) {
                const Vector& r_integration_point = integration_point_list_slave[i];
                const double integration_weight = integration_weight_list_slave[i];

                Vector distance_vector_slave = ZeroVector(mSlaveDim);
                distance_vector_slave[0] = r_integration_point[0] - slave_center.X();
                distance_vector_slave[1] = r_integration_point[1] - slave_center.Y();

                Matrix grad_H_sum_slave_transposed;
                ComputeGradientTaylorExpansionContribution(r_geometry_slave, distance_vector_slave, mSlaveBasisFunctionsOrder, grad_H_sum_slave_transposed);
                Matrix grad_H_sum_slave = trans(grad_H_sum_slave_transposed);

                Matrix B_slave;
                CalculateB(B_slave, grad_H_sum_slave, number_of_nodes_slave);

                Vector strain_slave = prod(B_slave, coefficient_slave);

                ConstitutiveLaw::Parameters values_slave_true(r_geometry_slave, GetPropertiesSlave(), rCurrentProcessInfo);
                ConstitutiveVariables constitutive_variables_slave(mpConstitutiveLawSlave->GetStrainSize());
                ApplyConstitutiveLaw(mpConstitutiveLawSlave, strain_slave, values_slave_true, constitutive_variables_slave);

                const Vector& stress_vector_slave_on_true = values_slave_true.GetStressVector();

                Vector normal_stress_slave_vector = ZeroVector(3);
                normal_stress_slave_vector[0] = stress_vector_slave_on_true[0] * normal_true_slave[0] + stress_vector_slave_on_true[2] * normal_true_slave[1];
                normal_stress_slave_vector[1] = stress_vector_slave_on_true[2] * normal_true_slave[0] + stress_vector_slave_on_true[1] * normal_true_slave[1];

                double normal_stress = normal_stress_slave_vector[0] * normal_true_slave[0] + normal_stress_slave_vector[1] * normal_true_slave[1];
                const double shear_stress = -normal_stress_slave_vector[0] * normal_true_slave[1] + normal_stress_slave_vector[1] * normal_true_slave[0];
                

                // FIXME: 
                auto gap = GetValue(GAP);
                double normal_gap = gap[0] * normal_true_master[0] + gap[1] * normal_true_master[1];
                auto penalty_master = (*mpPropMaster)[PENALTY_FACTOR];
                const double h = std::min(mMasterCharacteristicLength, mSlaveCharacteristicLength);

                // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
                // penalty_master = penalty_master / h * mMasterBasisFunctionsOrder *mMasterBasisFunctionsOrder;
                // if (normal_gap < 0.0) {
                //     normal_stress = penalty_master * normal_gap;
                // } else {
                //     normal_stress = 0;
                // }
                values_on_true_boundary_slave(i, 0) = integration_weight;
                values_on_true_boundary_slave(i, 1) = r_integration_point.size() > 0 ? r_integration_point[0] : 0.0;
                values_on_true_boundary_slave(i, 2) = r_integration_point.size() > 1 ? r_integration_point[1] : 0.0;
                values_on_true_boundary_slave(i, 3) = r_integration_point.size() > 2 ? r_integration_point[2] : 0.0;
                values_on_true_boundary_slave(i, 4) = normal_stress;
                values_on_true_boundary_slave(i, 5) = shear_stress;
            }
        }

        this->SetValue(RESULTS_ON_TRUE_BOUNDARY_MASTER, values_on_true_boundary_master);
        this->SetValue(RESULTS_ON_TRUE_BOUNDARY_SLAVE, values_on_true_boundary_slave);
        this->SetValue(RESULTS_ON_TRUE_BOUNDARY, Matrix());
    }

    void SbmContact2DCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){

        // InitializeMaterial();
        ConstitutiveLaw::Parameters constitutive_law_parameters_master(
            GetMasterGeometry(), (*mpPropMaster), rCurrentProcessInfo);

        mpConstitutiveLawMaster->InitializeMaterialResponse(constitutive_law_parameters_master, ConstitutiveLaw::StressMeasure_Cauchy);

        // InitializeMaterial(); //slave
        ConstitutiveLaw::Parameters constitutive_law_parameters_slave(
            GetSlaveGeometry(), (*mpPropSlave), rCurrentProcessInfo);

        mpConstitutiveLawSlave->InitializeMaterialResponse(constitutive_law_parameters_slave, ConstitutiveLaw::StressMeasure_Cauchy);
        
        this->InitializeNonLinearIteration(rCurrentProcessInfo);

        // #############################################################
        SetValue(ACTIVATION_LEVEL, 0);  
        SetValue(YOUNG_MODULUS_MASTER, (*mpPropMaster)[YOUNG_MODULUS]);
        SetValue(YOUNG_MODULUS_SLAVE, (*mpPropSlave)[YOUNG_MODULUS]);
    }

    /**
     * @brief 
     * 
     * @param rCurrentProcessInfo 
     */
    void SbmContact2DCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo){

        const auto& r_geometry_master = GetMasterGeometry();
        const auto& r_geometry_slave = GetSlaveGeometry();

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();
        const SizeType dim = mMasterDim;

        Vector coefficient_master(dim * number_of_nodes_master);
        GetValuesVector(coefficient_master, QuadraturePointCouplingGeometry2D<Point>::Master);

        Vector coefficient_slave(dim * number_of_nodes_slave);
        GetValuesVector(coefficient_slave, QuadraturePointCouplingGeometry2D<Point>::Slave);

        Matrix H_sum_master;
        Matrix H_sum_slave;
        ComputeTaylorSumContribution(r_geometry_master, mDistanceMaster, 2*mMasterBasisFunctionsOrder, H_sum_master);
        ComputeTaylorSumContribution(r_geometry_slave, mDistanceSlave, 2*mSlaveBasisFunctionsOrder, H_sum_slave);

        Matrix H_master_true = ZeroMatrix(dim, dim * number_of_nodes_master);
        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            for (IndexType idim = 0; idim < dim; ++idim) {
                H_master_true(idim, dim * i + idim) = H_sum_master(0, i);
            }
        }

        Vector displacement_master_true_sub = prod(H_master_true, coefficient_master);
        Vector displacement_master_true = ZeroVector(3);
        displacement_master_true[0] = displacement_master_true_sub[0];
        displacement_master_true[1] = displacement_master_true_sub[1];
        this->SetValue(DISPLACEMENT_MASTER, displacement_master_true);

        Matrix H_slave_true = ZeroMatrix(dim, dim * number_of_nodes_slave);
        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            for (IndexType idim = 0; idim < dim; ++idim) {
                H_slave_true(idim, dim * i + idim) = H_sum_slave(0, i);
            }
        }

        Vector displacement_slave_true_sub = prod(H_slave_true, coefficient_slave);
        Vector displacement_slave_true = ZeroVector(3);
        displacement_slave_true[0] = displacement_slave_true_sub[0];
        displacement_slave_true[1] = displacement_slave_true_sub[1];
        this->SetValue(DISPLACEMENT_SLAVE, displacement_slave_true);

        Matrix grad_H_sum_master_transposed;
        ComputeGradientTaylorExpansionContribution(r_geometry_master, mDistanceMaster, mMasterBasisFunctionsOrder, grad_H_sum_master_transposed);
        Matrix grad_H_sum_master = trans(grad_H_sum_master_transposed);

        Matrix B_sum_master = ZeroMatrix(3, number_of_nodes_master);
        CalculateB(B_sum_master, grad_H_sum_master, number_of_nodes_master);

        Vector strain_vector_master = prod(B_sum_master, coefficient_master);
        this->SetValue(STRAIN_MASTER, strain_vector_master);
        this->SetConstitutiveVariables(strain_vector_master, 0, rCurrentProcessInfo);

        Matrix grad_H_sum_slave_transposed;
        ComputeGradientTaylorExpansionContribution(r_geometry_slave, mDistanceSlave, mSlaveBasisFunctionsOrder, grad_H_sum_slave_transposed);
        Matrix grad_H_sum_slave = trans(grad_H_sum_slave_transposed);

        Matrix B_sum_slave = ZeroMatrix(3, number_of_nodes_slave);
        CalculateB(B_sum_slave, grad_H_sum_slave, number_of_nodes_slave);
        Vector strain_vector_slave = prod(B_sum_slave, coefficient_slave);
        this->SetValue(STRAIN_SLAVE, strain_vector_slave);
        this->SetConstitutiveVariables(strain_vector_slave, 1, rCurrentProcessInfo);

        SetGap();
    }

    /**
     * @brief 
     * 
     * @param rCurrentProcessInfo 
     */
    void SbmContact2DCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo){

        this->InitializeNonLinearIteration(rCurrentProcessInfo);
    }


    


    //----------------------------------------------------------------------------------
    const Matrix SbmContact2DCondition::GetConstitutiveMatrix(IndexType index, Matrix& r_B, GeometryType r_geometry,
                                                                   Vector& old_displacement, const Kratos::ProcessInfo& rCurrentProcessInfo,
                                                                   Vector& stress_vector) {

        ConstitutiveLaw::Pointer rpConstitutiveLaw = GetConstitutiveLaw(index);

        PropertiesType r_prop = GetProperty(index);

        ConstitutiveLaw::Parameters Values(r_geometry, r_prop, rCurrentProcessInfo);

        const SizeType strain_size = rpConstitutiveLaw->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector old_strain = prod(r_B,old_displacement);
    
        // Values.SetStrainVector(this_constitutive_variables.rStrainVector);
        Values.SetStrainVector(old_strain);
        Values.SetStressVector(this_constitutive_variables.StressVector);

        Values.SetConstitutiveMatrix(this_constitutive_variables.D);
        rpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

        
        stress_vector = Values.GetStressVector();

        return Values.GetConstitutiveMatrix();
        
    }



    void SbmContact2DCondition::GetDeformed(const Matrix& N, Vector& reference_position, Vector& displacement, Vector& deformed_position) {
        
        // compute the old displacement for the two points in the pair
        double displacement_x = 0.0; double displacement_y = 0.0;
        
        for (IndexType i = 0; i < GetMasterGeometry().size(); ++i)
        {
            // KRATOS_WATCH(r_geometry[i])
            double output_solution_step_value_x = displacement[2*i];
            double output_solution_step_value_y = displacement[2*i+1];
            displacement_x += N(0, i) * output_solution_step_value_x;
            displacement_y += N(0, i) * output_solution_step_value_y;
        } 
        deformed_position[0] = reference_position[0] + displacement_x; 
        deformed_position[1] = reference_position[1] + displacement_y;
    }


    void SbmContact2DCondition::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        if (rValues.size() != 1)
            rValues.resize(1);
        if (rVariable == INTEGRATION_WEIGHT)
            rValues[0] = this->GetValue(INTEGRATION_WEIGHT);
        else if (rVariable == NORMAL_GAP)
            rValues[0] = this->GetValue(NORMAL_GAP);
    }

    void SbmContact2DCondition::CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints not ready for the SBM!!";
        const SizeType dimension = 2;//GetMasterGeometry().WorkingSpaceDimension();

        if (rValues.size() != dimension)
            rValues.resize(dimension);
        if (rVariable == NORMAL_STRESS)
        {
            Vector normal_physical_space = this->GetValue(NORMAL);
            Vector sigma_voigt_master = GetValue(STRESS_MASTER);

            Vector normal_stress_master(2);

            normal_stress_master[0] = sigma_voigt_master[0]*normal_physical_space[0] + sigma_voigt_master[2]*normal_physical_space[1];
            normal_stress_master[1] = sigma_voigt_master[2]*normal_physical_space[0] + sigma_voigt_master[1]*normal_physical_space[1];
            
            this->SetValue(NORMAL_STRESS, normal_stress_master);

            rValues[0] = normal_stress_master;
        }
        else if (rVariable == NORMAL_MASTER) 
            rValues[0] = this->GetValue(NORMAL_MASTER);
    }




     // Function to compute a single term in the Taylor expansion
    double SbmContact2DCondition::ComputeTaylorTerm(const double derivative, const double dx, const int n_k, const double dy, const int k) const
    {
        return derivative * std::pow(dx, n_k) * std::pow(dy, k)
            / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
    }


    void SbmContact2DCondition::ApplyConstitutiveLaw(
        ConstitutiveLaw::Pointer pConstitutiveLaw,
        Vector& rStrain,
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiVariables) const
    {
        KRATOS_DEBUG_ERROR_IF(pConstitutiveLaw == nullptr)
            << "Null constitutive law pointer in SbmContact2DCondition::ApplyConstitutiveLaw" << std::endl;

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions = rValues.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        rValues.SetStrainVector(rStrain);
        rValues.SetStressVector(rConstitutiVariables.StressVector);
        rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

        pConstitutiveLaw->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy);
    }

} // Namespace Kratos
