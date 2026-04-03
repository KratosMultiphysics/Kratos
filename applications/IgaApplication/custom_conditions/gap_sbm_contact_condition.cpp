
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
#include "custom_conditions/gap_sbm_contact_condition.h"

namespace Kratos
{

void GapSbmContactCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    mSideIdentifier = GetValue(IDENTIFIER);
    KRATOS_ERROR_IF(mSideIdentifier != "MASTER" && mSideIdentifier != "SLAVE" && mSideIdentifier != "INACTIVE")
        << "\"GapSbmContactCondition\" #" << Id()
        << " has IDENTIFIER=\"" << mSideIdentifier
        << "\". Expected \"MASTER\", \"SLAVE\" or \"INACTIVE\"." << std::endl;

    InitializeMemberVariables();
    InitializeSbmMemberVariables();
    InitializeMaterial();
}


void GapSbmContactCondition::InitializeMaterial()
{
    KRATOS_TRY
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
            return;
        }
        const GeometryType& r_geometry = GetSurrogateGeometry();
        const Properties& r_properties = GetProperties();
        const SizeType number_of_control_points = r_geometry.size();
        Vector N_sum_vec = ZeroVector(number_of_control_points);
        ComputeTaylorExpansionContribution(
            r_geometry,
            mDistanceVectorSkinReferenceMaster,
            mBasisFunctionsOrderMaster,
            N_sum_vec);
        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, N_sum_vec);
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );

}

void GapSbmContactCondition::InitializeMemberVariables()
{
    // // Compute class memeber variables
    const auto& r_geometry = GetGeometry();

    const auto& r_projected_geometry = GetSurrogateGeometry();
    const auto& r_DN_De = r_projected_geometry.ShapeFunctionsLocalGradients(r_projected_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2) << "GapSbmContactCondition momentarily only supports 2D conditions, but the current dimension is" << mDim << std::endl;
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrderMaster = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrderMaster = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    mBasisFunctionsOrderMaster *= 2; 

    // Compute the normals
    mNormalParameterSpace = r_geometry.Normal(0, GetIntegrationMethod());
    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);
    mNormalPhysicalSpaceMaster = mNormalParameterSpace;

    SetValue(NORMAL, mNormalPhysicalSpaceMaster);

    // calculate the integration weight
    // reading integration point
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

    // const double integration_weight = r_integration_points[0].Weight() * std::abs(detJ0) * thickness;

    const double integration_weight = r_integration_points[0].Weight()*thickness;

    SetValue(INTEGRATION_WEIGHT, integration_weight);
}

void GapSbmContactCondition::InitializeSbmMemberVariables()
{
    const auto& r_geometry = this->GetGeometry();
    const auto& r_surrogate_geometry = GetSurrogateGeometry();

    // FIXME: simply remove
    if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
        return;
    }

    mDistanceVectorSkinReferenceMaster.resize(3);
    noalias(mDistanceVectorSkinReferenceMaster) = r_geometry.Center().Coordinates() - r_surrogate_geometry.Center().Coordinates();


    const auto& p_skin_node_slave = pGetProjectionNode();
    const auto& p_surrogate_geometry_slave = p_skin_node_slave->GetValue(NEIGHBOUR_GEOMETRIES)[0];


    mDistanceVectorSkinReferenceSlave.resize(3);
    noalias(mDistanceVectorSkinReferenceSlave) = p_skin_node_slave->Coordinates() - p_surrogate_geometry_slave->Center().Coordinates();
    mNormalPhysicalSpaceSlave = p_surrogate_geometry_slave->Normal(0, GetIntegrationMethod());

    mNormalPhysicalSpaceSlave = mNormalPhysicalSpaceSlave / MathUtils<double>::Norm(mNormalPhysicalSpaceSlave);

    mNormalPhysicalSpaceSlave = -mNormalPhysicalSpaceMaster; //TODO:


    const auto& r_DN_De_slave = p_surrogate_geometry_slave->ShapeFunctionsLocalGradients(
        p_surrogate_geometry_slave->GetDefaultIntegrationMethod());

    if (mDim == 3) {
        mBasisFunctionsOrderSlave = std::cbrt(r_DN_De_slave[0].size1()) - 1;
    } else {
        mBasisFunctionsOrderSlave = std::sqrt(r_DN_De_slave[0].size1()) - 1;
    }

    mBasisFunctionsOrderSlave *= 2;

    SetValue(SKIN_MASTER_COORDINATES, GetGeometry().Center().Coordinates());
    SetValue(SKIN_SLAVE_COORDINATES, pGetProjectionNode()->Coordinates());
}

void GapSbmContactCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    // FIXME: simply remove
    if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
        return;
    }

    const SizeType mat_size = GetSurrogateGeometry().size() * 2;

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    SetValue(SKIN_MASTER_COORDINATES, GetGeometry().Center().Coordinates());
    SetValue(SKIN_SLAVE_COORDINATES, pGetProjectionNode()->Coordinates());

    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void GapSbmContactCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    // FIXME: simply remove
    if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
        return;
    }

    const auto& r_skin_geometry_master = this->GetGeometry();
    const auto& r_surrogate_geometry_master = GetSurrogateGeometry();
    const unsigned int number_of_control_points_master = r_surrogate_geometry_master.size();

    const auto& r_skin_node_slave = *pGetProjectionNode();
    const auto& r_surrogate_geometry_slave = *r_skin_node_slave.GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const unsigned int number_of_control_points_slave = r_surrogate_geometry_slave.size();


    // reading integration points and local gradients
    const std::size_t mat_size_master = number_of_control_points_master * mDim;
    const std::size_t mat_size_slave = number_of_control_points_slave * mDim;
    const std::size_t mat_size = (number_of_control_points_master + number_of_control_points_slave) * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    if (this->GetValue(ACTIVATION_LEVEL) == 0)
    {
        return;
    }

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec_master = ZeroVector(number_of_control_points_master);
    ComputeTaylorExpansionContribution(
        r_surrogate_geometry_master,
        mDistanceVectorSkinReferenceMaster,
        mBasisFunctionsOrderMaster,
        N_sum_vec_master);

    // slave side
    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec_slave = ZeroVector(number_of_control_points_slave);
    ComputeTaylorExpansionContribution(
        r_surrogate_geometry_slave,
        mDistanceVectorSkinReferenceSlave,
        mBasisFunctionsOrderSlave,
        N_sum_vec_slave);

    // master true gradient by shifted basis functions
    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed_master = ZeroMatrix(3, number_of_control_points_master);
    ComputeGradientTaylorExpansionContribution(
        r_surrogate_geometry_master,
        mDistanceVectorSkinReferenceMaster,
        mBasisFunctionsOrderMaster,
        grad_N_sum_transposed_master);
    Matrix grad_N_sum_master = trans(grad_N_sum_transposed_master);

    Matrix B_sum_master = ZeroMatrix(mDim,mat_size_master);
    CalculateB(B_sum_master, grad_N_sum_master, number_of_control_points_master);

    // obtain the tangent constitutive matrix at the true position
    ConstitutiveLaw::Parameters values_true_master(r_skin_geometry_master, GetProperties(), rCurrentProcessInfo);
    Vector displacement_coefficient_vector_master(mat_size_master);
    GetSolutionCoefficientVector(displacement_coefficient_vector_master, 0);
    Vector strain_true_master = prod(B_sum_master, displacement_coefficient_vector_master);
    const std::size_t strain_size_true_master = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true_master(strain_size_true_master);
    ApplyConstitutiveLaw(mat_size_master, strain_true_master, values_true_master, this_constitutive_variables_true_master);

    const Matrix& r_D_on_true_master = values_true_master.GetConstitutiveMatrix();

    Matrix DB_sum_master = prod(r_D_on_true_master, B_sum_master);

    // Differential area
    mNitschePenalty = -1; //FIXME:
    // double penalty_integration = mPenalty * integration_weight;

    // ASSEMBLE
    const std::size_t shift_dof = mat_size_master;
    //-----------------------------------------------------
    // MASTER
    for (IndexType i = 0; i < number_of_control_points_master; i++) {
        for (IndexType idim = 0; idim < 2; idim++) {
            const int id1 = 2*idim;
            const int iglob = 2*i+idim;

            // PENALTY TERM
            // rLeftHandSideMatrix(iglob, 2*j+idim) += N_sum_vec_master(i)*N_sum_vec_master(j)* penalty_integration;

            Vector extended_sigma_w_n = ZeroVector(3);
            extended_sigma_w_n[0] = (DB_sum_master(0, iglob)* mNormalPhysicalSpaceMaster[0] + DB_sum_master(2, iglob)* mNormalPhysicalSpaceMaster[1]);
            extended_sigma_w_n[1] = (DB_sum_master(2, iglob)* mNormalPhysicalSpaceMaster[0] + DB_sum_master(1, iglob)* mNormalPhysicalSpaceMaster[1]);

            const double extended_sigma_w_n_dot_n = extended_sigma_w_n[0]*mNormalPhysicalSpaceMaster[0] + extended_sigma_w_n[1]*mNormalPhysicalSpaceMaster[1];
            for (IndexType j = 0; j < number_of_control_points_master; j++) {
                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob = 2*j+jdim;

                    Vector extension_sigma_u_n(2);
                    extension_sigma_u_n[0] = (DB_sum_master(0, jglob)* mNormalPhysicalSpaceMaster[0] + DB_sum_master(2, jglob)* mNormalPhysicalSpaceMaster[1]);
                    extension_sigma_u_n[1] = (DB_sum_master(2, jglob)* mNormalPhysicalSpaceMaster[0] + DB_sum_master(1, jglob)* mNormalPhysicalSpaceMaster[1]);

                    const double extension_sigma_u_n_dot_n = extension_sigma_u_n[0]*mNormalPhysicalSpaceMaster[0] + extension_sigma_u_n[1]*mNormalPhysicalSpaceMaster[1];

                    // FLUX 
                    // -[sigma(u) \dot n] \dot n * (w1 \dot n1 ....           + w2 \dot n2)
                    //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= extension_sigma_u_n_dot_n * N_sum_vec_master(i) * mNormalPhysicalSpaceMaster[idim] * integration_weight;

                    // // PENALTY FREE FOR NON PENETRABILITY g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (u_1 \dot n1 + u_2 \dot n2)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= mNitschePenalty *extended_sigma_w_n_dot_n
                                                         * (N_sum_vec_master(j) * mNormalPhysicalSpaceMaster[jdim]) * integration_weight;
                }
            }

            // SLAVE NON PENETRABILITY PART
            for (IndexType j = 0; j < number_of_control_points_master; j++) {
                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob = 2*j+jdim + shift_dof;

                    // // PENALTY FREE FOR NON PENETRABILITY g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (u_1 \dot n1 + u_2 \dot n2)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= mNitschePenalty *extended_sigma_w_n_dot_n
                                                         * (N_sum_vec_slave(j) * mNormalPhysicalSpaceSlave[jdim]) * integration_weight;
                }

            }
        }
    }

    //SLAVE
    for (IndexType i = 0; i < number_of_control_points_slave; i++) {
        for (IndexType idim = 0; idim < 2; idim++) {
            const int id1 = 2*idim;
            const int iloc = 2*i+idim;
            const int iglob = 2*i+idim + shift_dof;

            // PENALTY TERM
            // rLeftHandSideMatrix(iglob, 2*j+idim) += N_sum_vec_slave(i)*N_sum_vec_slave(j)* penalty_integration;

            for (IndexType j = 0; j < number_of_control_points_master; j++) {
                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob = 2*j+jdim;

                    Vector extension_sigma_u_n(2);
                    extension_sigma_u_n[0] = (DB_sum_master(0, jglob)* mNormalPhysicalSpaceMaster[0] + DB_sum_master(2, jglob)* mNormalPhysicalSpaceMaster[1]);
                    extension_sigma_u_n[1] = (DB_sum_master(2, jglob)* mNormalPhysicalSpaceMaster[0] + DB_sum_master(1, jglob)* mNormalPhysicalSpaceMaster[1]);

                    const double extension_sigma_u_n_dot_n = extension_sigma_u_n[0]*mNormalPhysicalSpaceMaster[0] + extension_sigma_u_n[1]*mNormalPhysicalSpaceMaster[1];

                    // FLUX 
                    // -[sigma(u) \dot n] \dot n * (w1 \dot n1 ....           + w2 \dot n2)
                    //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= extension_sigma_u_n_dot_n * N_sum_vec_slave(i) * mNormalPhysicalSpaceSlave[idim] * integration_weight;
                }

            }
        }
    }

    KRATOS_CATCH("")
}

void GapSbmContactCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    // FIXME: simply remove
    if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
        return;
    }

    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const auto& r_true_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();

    // reading integration points and local gradients
    const SizeType mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution(
        r_surrogate_geometry,
        mDistanceVectorSkinReferenceMaster,
        mBasisFunctionsOrderMaster,
        N_sum_vec);



    // FIXME:
    MatrixType lhs;
    CalculateLeftHandSide(lhs, rCurrentProcessInfo);

    const std::size_t size = lhs.size1();
    KRATOS_ERROR_IF(size == 0) << "SolidCouplingCondition::CalculateRightHandSide: The left hand side matrix has zero size." << std::endl;

    Vector solution_coefficient_vector;
    GetSolutionCoefficientVector(solution_coefficient_vector, 2);

    if (rRightHandSideVector.size() != size) {
        rRightHandSideVector.resize(size, false);
    }
    noalias(rRightHandSideVector) = -prod(lhs, solution_coefficient_vector);

    KRATOS_CATCH("")
}


    void GapSbmContactCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
            return;
        }
        const auto& r_surrogate_geometry_master = GetSurrogateGeometry();
        const SizeType number_of_control_points_master = r_surrogate_geometry_master.size();

        const auto& r_skin_node_slave = *pGetProjectionNode();
        const auto& r_surrogate_geometry_slave = *r_skin_node_slave.GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points_slave = r_surrogate_geometry_slave.size();

        const SizeType number_of_control_points = number_of_control_points_master + number_of_control_points_slave;

        if (rResult.size() != 2 * number_of_control_points)
            rResult.resize(2 * number_of_control_points, false);

        for (IndexType i = 0; i < number_of_control_points_master; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_surrogate_geometry_master[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }

        const IndexType shift = number_of_control_points_master * 2;
        for (IndexType i = 0; i < number_of_control_points_slave; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_surrogate_geometry_slave[i];
            rResult[index + shift] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1 + shift] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

    void GapSbmContactCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
            return;
        }
        const auto& r_surrogate_geometry_master = GetSurrogateGeometry();
        const SizeType number_of_control_points_master = r_surrogate_geometry_master.size();

        const auto& r_skin_node_slave = *pGetProjectionNode();
        const auto& r_surrogate_geometry_slave = *r_skin_node_slave.GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points_slave = r_surrogate_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * (number_of_control_points_master + number_of_control_points_slave));

        for (IndexType i = 0; i < number_of_control_points_master; ++i) {
            const auto& r_node = r_surrogate_geometry_master[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }

        for (IndexType i = 0; i < number_of_control_points_slave; ++i) {
            const auto& r_node = r_surrogate_geometry_slave[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }
    };


    void GapSbmContactCondition::GetSolutionCoefficientVector(
        Vector& rValues,
        IndexType index) const
    {
        if (index == 0) {
            const auto& r_surrogate_geometry_master = GetSurrogateGeometry();
            const SizeType number_of_control_points_master = r_surrogate_geometry_master.size();
            const SizeType mat_size_master = number_of_control_points_master * 2;

            if (rValues.size() != mat_size_master)
                rValues.resize(mat_size_master, false);

            for (IndexType i = 0; i < number_of_control_points_master; ++i)
            {
                const array_1d<double, 3>& displacement = r_surrogate_geometry_master[i].GetSolutionStepValue(DISPLACEMENT);
                const IndexType index = i * 2;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }
        } else if (index == 1) {
            const auto& r_skin_node_slave = *pGetProjectionNode();
            const auto& r_surrogate_geometry_slave = *r_skin_node_slave.GetValue(NEIGHBOUR_GEOMETRIES)[0];
            const SizeType number_of_control_points_slave = r_surrogate_geometry_slave.size();
            const SizeType mat_size_slave = number_of_control_points_slave * 2;

            if (rValues.size() != mat_size_slave)
                rValues.resize(mat_size_slave, false);

            for (IndexType i = 0; i < number_of_control_points_slave; ++i)
            {
                const array_1d<double, 3>& displacement = r_surrogate_geometry_slave[i].GetSolutionStepValue(DISPLACEMENT);
                const IndexType index = i * 2;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }
        } else if (index == 2) {
            const auto& r_surrogate_geometry_master = GetSurrogateGeometry();
            const SizeType number_of_control_points_master = r_surrogate_geometry_master.size();
            const SizeType mat_size_master = number_of_control_points_master * 2;

            const auto& r_skin_node_slave = *pGetProjectionNode();
            const auto& r_surrogate_geometry_slave = *r_skin_node_slave.GetValue(NEIGHBOUR_GEOMETRIES)[0];
            const SizeType number_of_control_points_slave = r_surrogate_geometry_slave.size();
            const SizeType mat_size_slave = number_of_control_points_slave * 2;

            const SizeType mat_size = mat_size_master + mat_size_slave;
            if (rValues.size() != mat_size)
                rValues.resize(mat_size, false);

            for (IndexType i = 0; i < number_of_control_points_master; ++i)
            {
                const array_1d<double, 3>& displacement = r_surrogate_geometry_master[i].GetSolutionStepValue(DISPLACEMENT);
                const IndexType index = i * 2;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }

            for (IndexType i = 0; i < number_of_control_points_slave; ++i)
            {
                const array_1d<double, 3>& displacement = r_surrogate_geometry_slave[i].GetSolutionStepValue(DISPLACEMENT);
                const IndexType index = i * 2 + mat_size_master;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }
        }
    }

    void GapSbmContactCondition::SetGap()
    {
        const auto& r_master_geometry = GetGeometry();
        const auto p_slave_node = pGetProjectionNode();

        Vector skin_master_coord(3);
        noalias(skin_master_coord) = r_master_geometry.Center().Coordinates();

        Vector skin_slave_coord(3);
        noalias(skin_slave_coord) = p_slave_node->Coordinates();

        Vector displacement_master = ZeroVector(3);
        Vector displacement_slave = ZeroVector(3);

        if (this->Has(DISPLACEMENT_MASTER)) {
            displacement_master = GetValue(DISPLACEMENT_MASTER);
        }
        if (this->Has(DISPLACEMENT_SLAVE)) {
            displacement_slave = GetValue(DISPLACEMENT_SLAVE);
        }

        const Vector skin_master_coord_deformed = skin_master_coord + displacement_master;
        const Vector skin_slave_coord_deformed = skin_slave_coord + displacement_slave;

        SetValue(SKIN_MASTER_COORDINATES, skin_master_coord);
        SetValue(SKIN_SLAVE_COORDINATES, skin_slave_coord);

        const Vector gap_deformed = skin_slave_coord_deformed - skin_master_coord_deformed;
        SetValue(GAP, gap_deformed);
    }

    void GapSbmContactCondition::CalculateB(
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

    void GapSbmContactCondition::ApplyConstitutiveLaw(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
                                        ConstitutiveVariables& rConstitutiVariables)
    {
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=rValues.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
        rValues.SetStrainVector(rStrain);
        rValues.SetStressVector(rConstitutiVariables.StressVector);
        rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

        mpConstitutiveLaw->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy); 
    }


    void GapSbmContactCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
            return;
        }
        const auto& r_surrogate_geometry_master = GetSurrogateGeometry();
        const SizeType number_of_control_points_master = r_surrogate_geometry_master.size();
        const SizeType mat_size_master = number_of_control_points_master * mDim;

        const auto& r_skin_node_slave = *pGetProjectionNode();
        const auto& r_surrogate_geometry_slave = *r_skin_node_slave.GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points_slave = r_surrogate_geometry_slave.size();
        const SizeType mat_size_slave = number_of_control_points_slave * mDim;

        ConstitutiveLaw::Parameters constitutive_law_parameters(
            r_surrogate_geometry_master, GetProperties(), rCurrentProcessInfo);

        mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
        Vector old_displacement(mat_size_master);
        GetSolutionCoefficientVector(old_displacement, 0);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points_master);
        ComputeGradientTaylorExpansionContribution(
            r_surrogate_geometry_master,
            mDistanceVectorSkinReferenceMaster,
            mBasisFunctionsOrderMaster,
            grad_N_sum_transposed);
        Matrix grad_N_sum = trans(grad_N_sum_transposed);

        Matrix B_sum = ZeroMatrix(mDim,mat_size_master);
        CalculateB(B_sum, grad_N_sum, number_of_control_points_master);

        // obtain the tangent constitutive matrix at the true position
        ConstitutiveLaw::Parameters values_true(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        Vector old_displacement_coefficient_vector(mat_size_master);
        GetSolutionCoefficientVector(old_displacement_coefficient_vector, 0);
        Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);

        const SizeType strain_size_true = mpConstitutiveLaw->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
        ApplyConstitutiveLaw(mat_size_master, old_strain_on_true, values_true, this_constitutive_variables_true);

        const Vector sigma = values_true.GetStressVector();
        Vector sigma_n = ZeroVector(3);

        sigma_n[0] = sigma[0]*mNormalPhysicalSpaceMaster[0] + sigma[2]*mNormalPhysicalSpaceMaster[1];
        sigma_n[1] = sigma[2]*mNormalPhysicalSpaceMaster[0] + sigma[1]*mNormalPhysicalSpaceMaster[1];

        SetValue(NORMAL_STRESS, sigma_n);

        SetValue(CAUCHY_STRESS_XX, sigma[0]);
        SetValue(CAUCHY_STRESS_YY, sigma[1]);
        SetValue(CAUCHY_STRESS_XY, sigma[2]);
        // //---------------------
    }

void GapSbmContactCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
        return;
    }
    //--------------------------------------------------------------------------------------------
    // calculate the constitutive law response
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
    SetValue(ACTIVATION_LEVEL, 0);  
}

void GapSbmContactCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    // FIXME: simply remove
    if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
        return;
    }

    const auto& r_surrogate_geometry_master = GetSurrogateGeometry();
    const SizeType number_of_nodes_master = r_surrogate_geometry_master.size();

    const auto& r_skin_node_slave = *pGetProjectionNode();
    const auto& r_surrogate_geometry_slave = *r_skin_node_slave.GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const SizeType number_of_nodes_slave = r_surrogate_geometry_slave.size();

    const SizeType dim = mDim;

    Vector coefficient_master(dim * number_of_nodes_master);
    GetSolutionCoefficientVector(coefficient_master, 0);

    Vector coefficient_slave(dim * number_of_nodes_slave);
    GetSolutionCoefficientVector(coefficient_slave, 1);

    Vector H_sum_master;
    ComputeTaylorExpansionContribution(
        r_surrogate_geometry_master,
        mDistanceVectorSkinReferenceMaster,
        mBasisFunctionsOrderMaster,
        H_sum_master);
        
    Vector H_sum_slave;
    ComputeTaylorExpansionContribution(
        r_surrogate_geometry_slave,
        mDistanceVectorSkinReferenceSlave,
        mBasisFunctionsOrderSlave,
        H_sum_slave);

    Matrix H_master_true = ZeroMatrix(dim, dim * number_of_nodes_master);
    for (IndexType i = 0; i < number_of_nodes_master; ++i) {
        for (IndexType idim = 0; idim < dim; ++idim) {
            H_master_true(idim, dim * i + idim) = H_sum_master(i);
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
            H_slave_true(idim, dim * i + idim) = H_sum_slave(i);
        }
    }

    Vector displacement_slave_true_sub = prod(H_slave_true, coefficient_slave);
    Vector displacement_slave_true = ZeroVector(3);
    displacement_slave_true[0] = displacement_slave_true_sub[0];
    displacement_slave_true[1] = displacement_slave_true_sub[1];
    this->SetValue(DISPLACEMENT_SLAVE, displacement_slave_true);

    Matrix grad_H_sum_master_transposed;
    ComputeGradientTaylorExpansionContribution(
        r_surrogate_geometry_master,
        mDistanceVectorSkinReferenceMaster,
        mBasisFunctionsOrderMaster,
        grad_H_sum_master_transposed);
    Matrix grad_H_sum_master = trans(grad_H_sum_master_transposed);

    Matrix B_sum_master;
    CalculateB(B_sum_master, grad_H_sum_master, number_of_nodes_master);

    Vector strain_vector_master = prod(B_sum_master, coefficient_master);
    this->SetValue(STRAIN_MASTER, strain_vector_master);

    ConstitutiveLaw::Parameters values_master(r_surrogate_geometry_master, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables_master(mpConstitutiveLaw->GetStrainSize());
    ApplyConstitutiveLaw(dim * number_of_nodes_master, strain_vector_master, values_master, constitutive_variables_master);
    this->SetValue(STRESS_MASTER, values_master.GetStressVector());

    Matrix grad_H_sum_slave_transposed;
    ComputeGradientTaylorExpansionContribution(
        r_surrogate_geometry_slave,
        mDistanceVectorSkinReferenceSlave,
        mBasisFunctionsOrderSlave,
        grad_H_sum_slave_transposed);
    Matrix grad_H_sum_slave = trans(grad_H_sum_slave_transposed);

    Matrix B_sum_slave;
    CalculateB(B_sum_slave, grad_H_sum_slave, number_of_nodes_slave);

    Vector strain_vector_slave = prod(B_sum_slave, coefficient_slave);
    this->SetValue(STRAIN_SLAVE, strain_vector_slave);

    ConstitutiveLaw::Parameters values_slave(r_surrogate_geometry_slave, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables_slave(mpConstitutiveLaw->GetStrainSize());
    ApplyConstitutiveLaw(dim * number_of_nodes_slave, strain_vector_slave, values_slave, constitutive_variables_slave);
    this->SetValue(STRESS_SLAVE, values_slave.GetStressVector());

    this->SetValue(NORMAL_SKIN_MASTER, mNormalPhysicalSpaceMaster);
    this->SetValue(NORMAL_SKIN_SLAVE, mNormalPhysicalSpaceSlave);

    SetGap();
}

void GapSbmContactCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    this->InitializeNonLinearIteration(rCurrentProcessInfo);
}

void GapSbmContactCondition::ComputeTaylorExpansionContribution(
    const GeometryType& rGeometry,
    const Vector& rDistanceVector,
    const SizeType BasisFunctionsOrder,
    Vector& H_sum_vec)
{
    const SizeType number_of_control_points = rGeometry.PointsNumber();
    const Matrix& r_N = rGeometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    
    if (H_sum_vec.size() != number_of_control_points) {
        H_sum_vec.resize(number_of_control_points, false);
    }
    noalias(H_sum_vec) = ZeroVector(number_of_control_points);

    if (BasisFunctionsOrder == 0) {
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            H_sum_vec(i) = r_N(0, i);
        }
        return;
    }

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(BasisFunctionsOrder);
    for (IndexType n = 1; n <= BasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        double H_taylor_term = 0.0;

        if (mDim == 2) {
            for (IndexType n = 1; n <= BasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n; k++) {
                    const IndexType n_k = n - k;
                    const double derivative = r_shape_function_derivatives(i, k);
                    H_taylor_term += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType n = 1; n <= BasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                IndexType countDerivativeId = 0;
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        const IndexType k_z = n - k_x - k_y;
                        const double derivative = r_shape_function_derivatives(i, countDerivativeId);
                        H_taylor_term += ComputeTaylorTerm3D(
                            derivative,
                            rDistanceVector[0],
                            k_x,
                            rDistanceVector[1],
                            k_y,
                            rDistanceVector[2],
                            k_z);
                        countDerivativeId++;
                    }
                }
            }
        }

        H_sum_vec(i) = H_taylor_term + r_N(0, i);
    }
}

void GapSbmContactCondition::ComputeGradientTaylorExpansionContribution(
    const GeometryType& rGeometry,
    const Vector& rDistanceVector,
    const SizeType BasisFunctionsOrder,
    Matrix& grad_H_sum)
{
    const SizeType number_of_control_points = rGeometry.PointsNumber();
    const auto& r_DN_De = rGeometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(BasisFunctionsOrder);
    for (IndexType n = 1; n <= BasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    if (grad_H_sum.size1() != 3 || grad_H_sum.size2() != number_of_control_points) {
        grad_H_sum.resize(3, number_of_control_points);
    }

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        double H_taylor_term_X = 0.0;
        double H_taylor_term_Y = 0.0;
        double H_taylor_term_Z = 0.0;

        if (mDim == 2) {
            for (IndexType n = 2; n <= BasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n-1; k++) {
                    const IndexType n_k = n - 1 - k;
                    const double derivative = shapeFunctionDerivatives(i, k);
                    H_taylor_term_X += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n-1; k++) {
                    const IndexType n_k = n - 1 - k;
                    const double derivative = shapeFunctionDerivatives(i, k + 1);
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType n = 2; n <= BasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
                IndexType countDerivativeId = 0;
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        const IndexType k_z = n - k_x - k_y;
                        const double derivative = shapeFunctionDerivatives(i, countDerivativeId);

                        if (k_x >= 1) {
                            H_taylor_term_X += ComputeTaylorTerm3D(
                                derivative,
                                rDistanceVector[0],
                                k_x - 1,
                                rDistanceVector[1],
                                k_y,
                                rDistanceVector[2],
                                k_z);
                        }
                        if (k_y >= 1) {
                            H_taylor_term_Y += ComputeTaylorTerm3D(
                                derivative,
                                rDistanceVector[0],
                                k_x,
                                rDistanceVector[1],
                                k_y - 1,
                                rDistanceVector[2],
                                k_z);
                        }
                        if (k_z >= 1) {
                            H_taylor_term_Z += ComputeTaylorTerm3D(
                                derivative,
                                rDistanceVector[0],
                                k_x,
                                rDistanceVector[1],
                                k_y,
                                rDistanceVector[2],
                                k_z - 1);
                        }
                        countDerivativeId++;
                    }
                }
            }
        }

        grad_H_sum(0, i) = H_taylor_term_X + r_DN_De[0](i, 0);
        grad_H_sum(1, i) = H_taylor_term_Y + r_DN_De[0](i, 1);
        if (mDim == 3)
            grad_H_sum(2, i) = H_taylor_term_Z + r_DN_De[0](i, 2);
        else
            grad_H_sum(2, i) = 0;
    }
}

// Function to compute a single term in the Taylor expansion
double GapSbmContactCondition::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    const double result = derivative * std::pow(dx, n_k) * std::pow(dy, k)
        / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));
    return result;
}

double GapSbmContactCondition::ComputeTaylorTerm3D(
    const double derivative, 
    const double dx, 
    const IndexType k_x, 
    const double dy, 
    const IndexType k_y, 
    const double dz, 
    const IndexType k_z)
{   
    const double result = derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z)
        / (MathUtils<double>::Factorial(k_x) * MathUtils<double>::Factorial(k_y) * MathUtils<double>::Factorial(k_z));
    return result;
}

const GapSbmContactCondition::NodeType::Pointer GapSbmContactCondition::pGetProjectionNode() const
{
    KRATOS_ERROR_IF_NOT(GetGeometry().Has(PROJECTION_NODE))
        << "\"GapSbmContactCondition\" #" << Id() << " has no PROJECTION_NODE set. "
        << "Geometry center: " << GetGeometry().Center() << std::endl;
    return GetGeometry().GetValue(PROJECTION_NODE);
}

const GapSbmContactCondition::GeometryType& GapSbmContactCondition::GetProjectionReferenceGeometry() const
{
    auto p_projection_node = pGetProjectionNode();
    KRATOS_ERROR_IF_NOT(p_projection_node->Has(NEIGHBOUR_GEOMETRIES))
        << "\"GapSbmContactCondition\" #" << Id() << " projection node has no NEIGHBOUR_GEOMETRIES." << std::endl;
    const auto& r_neighbour_geometries = p_projection_node->GetValue(NEIGHBOUR_GEOMETRIES);
    KRATOS_ERROR_IF(r_neighbour_geometries.empty())
        << "\"GapSbmContactCondition\" #" << Id() << " projection node has empty NEIGHBOUR_GEOMETRIES." << std::endl;
    return *r_neighbour_geometries[0];
}

} // Namespace Kratos
