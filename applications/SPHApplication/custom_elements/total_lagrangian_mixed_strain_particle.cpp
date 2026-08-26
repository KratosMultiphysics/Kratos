#include "custom_elements/total_lagrangian_mixed_strain_particle.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{

template<class TKernelType, std::size_t TDim>
Element::Pointer TotalLagrangianMixedStrainParticle<TKernelType, TDim>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
) const
{
    KRATOS_TRY
    
    TotalLagrangianMixedStrainParticle<TKernelType, TDim>::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangianMixedStrainParticle<TKernelType, TDim>>
        (NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    p_new_elem->SetConstitutiveLaw(this->mThisConstitutiveLaw);

    return p_new_elem;

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType dofs_per_node = TDim + (TDim * TDim);

    if (rResult.size() != dofs_per_node * number_of_neighbours)
        rResult.resize(dofs_per_node * number_of_neighbours, false);

    if constexpr (TDim == 2){
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * TDim;
            const SizeType F_block = i * TDim * TDim + TDim * number_of_neighbours;

            const IndexType vpos = r_geom[0].GetDofPosition(VELOCITY_X);
            const IndexType Fpos = r_geom[0].GetDofPosition(DEFORMATION_GRADIENT_XX);

            rResult[v_block    ] = r_geom[0].GetDof(VELOCITY_X, vpos    ).EquationId();
            rResult[v_block + 1] = r_geom[0].GetDof(VELOCITY_Y, vpos + 1).EquationId();
            rResult[F_block    ] = r_geom[0].GetDof(DEFORMATION_GRADIENT_XX, Fpos    ).EquationId();
            rResult[F_block + 1] = r_geom[0].GetDof(DEFORMATION_GRADIENT_YY, Fpos + 1).EquationId();
            rResult[F_block + 2] = r_geom[0].GetDof(DEFORMATION_GRADIENT_XY, Fpos + 2).EquationId();
            rResult[F_block + 3] = r_geom[0].GetDof(DEFORMATION_GRADIENT_YX, Fpos + 3).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * TDim;
            const SizeType F_block = i * TDim * TDim + TDim * number_of_neighbours;

            const IndexType vpos = r_geom[0].GetDofPosition(VELOCITY_X);
            const IndexType Fpos = r_geom[0].GetDofPosition(DEFORMATION_GRADIENT_XX);

            rResult[v_block    ] = r_geom[0].GetDof(VELOCITY_X, vpos    ).EquationId();
            rResult[v_block + 1] = r_geom[0].GetDof(VELOCITY_Y, vpos + 1).EquationId();
            rResult[v_block + 2] = r_geom[0].GetDof(VELOCITY_Z, vpos + 2).EquationId();
            rResult[F_block    ] = r_geom[0].GetDof(DEFORMATION_GRADIENT_XX, Fpos    ).EquationId();
            rResult[F_block + 1] = r_geom[0].GetDof(DEFORMATION_GRADIENT_YY, Fpos + 1).EquationId();
            rResult[F_block + 2] = r_geom[0].GetDof(DEFORMATION_GRADIENT_ZZ, Fpos + 2).EquationId(); 
            rResult[F_block + 3] = r_geom[0].GetDof(DEFORMATION_GRADIENT_XY, Fpos + 3).EquationId();
            rResult[F_block + 4] = r_geom[0].GetDof(DEFORMATION_GRADIENT_XZ, Fpos + 4).EquationId();
            rResult[F_block + 5] = r_geom[0].GetDof(DEFORMATION_GRADIENT_YX, Fpos + 5).EquationId();
            rResult[F_block + 6] = r_geom[0].GetDof(DEFORMATION_GRADIENT_YZ, Fpos + 6).EquationId();
            rResult[F_block + 7] = r_geom[0].GetDof(DEFORMATION_GRADIENT_ZX, Fpos + 7).EquationId();
            rResult[F_block + 8] = r_geom[0].GetDof(DEFORMATION_GRADIENT_ZY, Fpos + 8).EquationId(); 
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType dofs_per_node = TDim + (TDim * TDim);

    if(rElementalDofList.size() != dofs_per_node * number_of_neighbours)
        rElementalDofList.resize(dofs_per_node * number_of_neighbours);
    
    if constexpr (TDim == 2){
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * TDim;
            const SizeType F_block = i * TDim * TDim + TDim * number_of_neighbours;

            rElementalDofList[v_block    ] = r_geom[0].pGetDof(VELOCITY_X);
            rElementalDofList[v_block + 1] = r_geom[0].pGetDof(VELOCITY_Y);
            rElementalDofList[F_block    ] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_XX);
            rElementalDofList[F_block + 1] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_YY);
            rElementalDofList[F_block + 2] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_XY);
            rElementalDofList[F_block + 3] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_YX);
        }
    } else {
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * TDim;
            const SizeType F_block = i * TDim * TDim + TDim * number_of_neighbours;

            rElementalDofList[v_block    ] = r_geom[0].pGetDof(VELOCITY_X);
            rElementalDofList[v_block + 1] = r_geom[0].pGetDof(VELOCITY_Y);
            rElementalDofList[v_block + 2] = r_geom[0].pGetDof(VELOCITY_Z);
            rElementalDofList[F_block    ] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_XX);
            rElementalDofList[F_block + 1] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_YY);
            rElementalDofList[F_block + 2] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_ZZ); 
            rElementalDofList[F_block + 3] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_XY);
            rElementalDofList[F_block + 4] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_XZ);
            rElementalDofList[F_block + 5] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_YX);
            rElementalDofList[F_block + 6] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_YZ);
            rElementalDofList[F_block + 7] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_ZX);
            rElementalDofList[F_block + 8] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_ZY); 
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_TRY
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType dofs_per_node = TDim + (TDim * TDim);
    const SizeType mat_size = dofs_per_node * number_of_neighbours;

    if (rValues.size() != mat_size) 
        rValues.resize(mat_size, false);

    if constexpr (TDim == 2){
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * TDim;
            const SizeType F_block = i * TDim * TDim + TDim * number_of_neighbours;

            const array_1d<double, 3>& velocity = r_geom[0].FastGetSolutionStepValue(VELOCITY, Step);
            
            rValues[v_block    ] = velocity[0];
            rValues[v_block + 1] = velocity[1];
            rValues[F_block    ] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XX, Step);
            rValues[F_block + 1] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YY, Step);
            rValues[F_block + 2] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XY, Step);
            rValues[F_block + 3] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YX, Step);
        }
    } else {
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * TDim;
            const SizeType F_block = i * TDim * TDim + TDim * number_of_neighbours;

            const array_1d<double, 3>& velocity = r_geom[0].FastGetSolutionStepValue(VELOCITY, Step);

            rValues[v_block    ] = velocity[0];
            rValues[v_block + 1] = velocity[1];
            rValues[v_block + 2] = velocity[2];
            rValues[F_block    ] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XX, Step);
            rValues[F_block + 1] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YY, Step);
            rValues[F_block + 2] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_ZZ, Step); 
            rValues[F_block + 3] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XY, Step);
            rValues[F_block + 4] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XZ, Step);
            rValues[F_block + 5] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YX, Step);
            rValues[F_block + 6] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YZ, Step);
            rValues[F_block + 7] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_ZX, Step);
            rValues[F_block + 8] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_ZY, Step);
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::GetFirstDerivativesVector(VectorType& rValues, int step) const
{
    KRATOS_TRY
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType dofs_per_node = TDim + (TDim * TDim);
    const SizeType mat_size = dofs_per_node * number_of_neighbours;

    if (rValues.size() != mat_size) rValues.resize(mat_size, false);

    if constexpr (TDim == 2){
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * TDim;
            const SizeType F_block = i * TDim * TDim + TDim * number_of_neighbours;

            const array_1d<double, 3>& acceleration = r_geom[0].FastGetSolutionStepValue(ACCELERATION, step);
            
            rValues[v_block    ] = acceleration[0];
            rValues[v_block + 1] = acceleration[1];
            rValues[F_block    ] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_XX, step);
            rValues[F_block + 1] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_YY, step);
            rValues[F_block + 2] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_XY, step);
            rValues[F_block + 3] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_YX, step);
        }
    } else {
        for(IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * TDim;
            const SizeType F_block = i * TDim * TDim + TDim * number_of_neighbours;

            const array_1d<double, 3>& acceleration = r_geom[0].FastGetSolutionStepValue(ACCELERATION, step);

            rValues[v_block    ] = acceleration[0];
            rValues[v_block + 1] = acceleration[1];
            rValues[v_block + 2] = acceleration[2];
            rValues[F_block    ] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_XX, step);
            rValues[F_block + 1] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_YY, step);
            rValues[F_block + 2] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_ZZ, step);
            rValues[F_block + 3] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_XY, step);
            rValues[F_block + 4] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_XZ, step);
            rValues[F_block + 5] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_YX, step);
            rValues[F_block + 6] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_YZ, step);
            rValues[F_block + 7] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_ZX, step);
            rValues[F_block + 8] = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_ZY, step);
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateAll(
    MatrixType& rLHS, 
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const auto& r_geom = this->GetGeometry(); 
    const auto& r_prop = this->GetProperties();
    const SizeType strain_size = this->mThisConstitutiveLaw->GetStrainSize();
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType block_size = TDim + TDim * TDim; 
    const SizeType mat_size = block_size * number_of_neighbours;

    double theta = 1.0;
    double inverse_delta_time = 0.0;
    MatrixType mass_matrix;

    theta = rProcessInfo[TIME_INTEGRATION_THETA];
    const double delta_time = rProcessInfo[DELTA_TIME];
    inverse_delta_time = 1.0 / delta_time;
    CalculateMassMatrix(mass_matrix, rProcessInfo);

    if (CalculateStiffnessMatrixFlag){
        if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) rLHS.resize(mat_size, mat_size, false);
        noalias(rLHS) = ZeroMatrix(mat_size, mat_size);
    }

    if (CalculateResidualVectorFlag){
        if (rRHS.size() != mat_size) rRHS.resize(mat_size, false);
        noalias(rRHS) = ZeroVector(mat_size);
    }

    KinematicVariables this_kinematic_variables(strain_size, TDim, number_of_neighbours);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    ConstitutiveLaw::Pointer p_historical_constitutive_law;
    if (CalculateResidualVectorFlag && theta < 1.0) {
        p_historical_constitutive_law = this->mThisConstitutiveLaw->Clone();
    }

    ConstitutiveLaw::Parameters cl_values(r_geom, r_prop, rProcessInfo); 
    Flags& r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if (CalculateStiffnessMatrixFlag){
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }
    cl_values.SetStrainVector(this_constitutive_variables.StrainVector);

    const double gauss_weight = r_geom[0].GetValue(VOLUME);

    CalculateKinematicVariables(this_kinematic_variables, rProcessInfo);
    this->CalculateConstitutiveVariables(this_constitutive_variables, this_kinematic_variables, cl_values, ConstitutiveLaw::StressMeasure_PK2);

    if (CalculateStiffnessMatrixFlag){
        // Initializing the blocks
        MatrixType K11(TDim * number_of_neighbours, TDim * number_of_neighbours); K11.clear(); 
        MatrixType K12(TDim * number_of_neighbours, TDim * TDim * number_of_neighbours); K12.clear();
        MatrixType K21(TDim * TDim * number_of_neighbours, TDim * number_of_neighbours); K21.clear();
        MatrixType K22(TDim * TDim * number_of_neighbours, TDim * TDim * number_of_neighbours); K22.clear();

        
        CalculateAndAddKg(K12, this_kinematic_variables.DW_DX, this_constitutive_variables.StressVector, gauss_weight);
        CalculateAndAddKm(K12, this_kinematic_variables, this_constitutive_variables.C, gauss_weight);

        CalculateGeometricalTangentMatrix(K21, this_kinematic_variables, gauss_weight);
        
        if (rProcessInfo[DISSIPATION_COEFFICIENT] != 0.0)
            CalculateAndAddUpwindStabilizationTangent(K11, this_kinematic_variables, rProcessInfo);

        K11 *= theta; K12 *= theta; K21 *= -theta;

        AssembleLHS(rLHS, K11, K12, K21, K22);

        noalias(rLHS) += inverse_delta_time * mass_matrix;
    }

    if (CalculateResidualVectorFlag){
        VectorType RHSv(TDim * number_of_neighbours); RHSv.clear();
        VectorType RHSF(TDim * TDim * number_of_neighbours); RHSF.clear();

        CalculateLinearMomentumResidualVector(RHSv, this_kinematic_variables, rProcessInfo, this_constitutive_variables.StressVector, gauss_weight);
        CalculateGeometricalResidualVector(RHSF, this_kinematic_variables, gauss_weight);
        if (rProcessInfo[DISSIPATION_COEFFICIENT] != 0.0)
            CalculateAndAddUpwindStabilizationResidual(RHSv, this_kinematic_variables, rProcessInfo);

        AssembleRHS(rRHS, RHSv, RHSF);
        rRHS *= theta;

        if (theta < 1.0) {
            KinematicVariables historical_kinematic_variables(strain_size, TDim, number_of_neighbours);
            ConstitutiveVariables historical_constitutive_variables(strain_size);
            ConstitutiveLaw::Parameters historical_cl_values(r_geom, r_prop, rProcessInfo);
            Flags& r_historical_cl_options = historical_cl_values.GetOptions();
            r_historical_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
            r_historical_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            r_historical_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
            historical_cl_values.SetStrainVector(historical_constitutive_variables.StrainVector);

            CalculateKinematicVariables(historical_kinematic_variables, rProcessInfo, 1);
            this->SetConstitutiveLawVariables(historical_constitutive_variables, historical_kinematic_variables, historical_cl_values);
            p_historical_constitutive_law->CalculateMaterialResponse(historical_cl_values, ConstitutiveLaw::StressMeasure_PK2);

            VectorType historical_RHSv(TDim * number_of_neighbours); historical_RHSv.clear();
            VectorType historical_RHSF(TDim * TDim * number_of_neighbours); historical_RHSF.clear();
            VectorType historical_spatial_rhs(mat_size); historical_spatial_rhs.clear();

            CalculateLinearMomentumResidualVector(historical_RHSv, historical_kinematic_variables, rProcessInfo, historical_constitutive_variables.StressVector, gauss_weight, 1);
            CalculateGeometricalResidualVector(historical_RHSF, historical_kinematic_variables, gauss_weight, 1);
            if (rProcessInfo[DISSIPATION_COEFFICIENT] != 0.0)
                CalculateAndAddUpwindStabilizationResidual(historical_RHSv, historical_kinematic_variables, rProcessInfo, 1);

            AssembleRHS(historical_spatial_rhs, historical_RHSv, historical_RHSF);

            noalias(rRHS) += (1.0 - theta) * historical_spatial_rhs;
        }

        VectorType current_values, previous_values;
        GetValuesVector(current_values, 0);
        GetValuesVector(previous_values, 1);
        noalias(rRHS) -= inverse_delta_time * prod(mass_matrix, current_values - previous_values);
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateAndAddKg(
    MatrixType& rK12, 
    const Matrix& rDW_DX, 
    const Vector& rStressVector, 
    const double weight
    ) const 
{
    KRATOS_TRY
    // This should be okay, check if it is correct to use A to derive all the entries of the matrix 
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);

    int self_index = this->GetNeighbourPosition(r_neighbours);
    const SizeType column_start = TDim * TDim * self_index;

    VectorType A = ZeroVector(TDim);
    SizeType F_index;

    for (IndexType i = 0; i < r_neighbours.size(); ++i){

        const SizeType row_start = TDim * i;

        if constexpr (TDim == 2){

            A[0] = rStressVector[0] * rDW_DX(i, 0) + rStressVector[2] * rDW_DX(i, 1);
            A[1] = rStressVector[2] * rDW_DX(i, 0) + rStressVector[1] * rDW_DX(i, 1);
            A *= weight;

        } else {

            A[0] = rStressVector[0] * rDW_DX(i, 0) + rStressVector[3] * rDW_DX(i, 1) + rStressVector[5] * rDW_DX(i, 2); 
            A[1] = rStressVector[3] * rDW_DX(i, 0) + rStressVector[1] * rDW_DX(i, 1) + rStressVector[4] * rDW_DX(i, 2); 
            A[2] = rStressVector[5] * rDW_DX(i, 0) + rStressVector[4] * rDW_DX(i, 1) + rStressVector[2] * rDW_DX(i, 2);
            A *= weight;

        }

        for (SizeType row = 0; row < TDim; ++row){
            for (SizeType col = 0; col < TDim; ++col){

                if constexpr (TDim == 2) {
                    constexpr SizeType indices[2][2] = {{0, 2}, {3, 1}};
                    F_index = indices[row][col];
                } else {
                    constexpr SizeType indices[3][3] = {{0, 3, 4}, {5, 1, 6}, {7, 8, 2}};
                    F_index = indices[row][col];
                }
                rK12(row_start + row, column_start + F_index) += A[col];
            }
        }
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateAndAddKm(
    MatrixType& rK12,
    const KinematicVariables& rThisKinematicVariables,
    const MatrixType& rConstitutiveMatrix,
    const double weight
    ) const 
{
    KRATOS_TRY
    // This is okay, only thing to check is the correcteness of the L matrix (also in case someday I change variable order)

    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const int self_index = this->GetNeighbourPosition(r_neighbours);

    MatrixType L;
    const MatrixType& rF = rThisKinematicVariables.F;
    if constexpr (TDim == 2) {
        L = ZeroMatrix(3, 4);
        L(0, 0) = rF(0, 0); L(0, 3) = rF(1, 0);
        L(1, 1) = rF(1, 1); L(1, 2) = rF(0, 1);
        L(2, 0) = rF(0, 1); L(2, 1) = rF(1, 0);
        L(2, 2) = rF(0, 0); L(2, 3) = rF(1, 1);
    } else {
        L = ZeroMatrix(6, 9);
        L(0, 0) = rF(0, 0); L(0, 5) = rF(1, 0); L(0, 7) = rF(2, 0);
        L(1, 3) = rF(0, 1); L(1, 1) = rF(1, 1); L(1, 8) = rF(2, 1);
        L(2, 4) = rF(0, 2); L(2, 6) = rF(1, 2); L(2, 2) = rF(2, 2);
        L(3, 0) = rF(0, 1); L(3, 3) = rF(0, 0); L(3, 5) = rF(1, 1);
        L(3, 1) = rF(1, 0); L(3, 7) = rF(2, 1); L(3, 8) = rF(2, 0);
        L(4, 3) = rF(0, 2); L(4, 4) = rF(0, 1); L(4, 1) = rF(1, 2);
        L(4, 6) = rF(1, 1); L(4, 8) = rF(2, 2); L(4, 2) = rF(2, 1);
        L(5, 0) = rF(0, 2); L(5, 4) = rF(0, 0); L(5, 5) = rF(1, 2);
        L(5, 6) = rF(1, 0); L(5, 7) = rF(2, 2); L(5, 2) = rF(2, 0);
    }

    const MatrixType material_tangent = weight * prod(trans(rThisKinematicVariables.B), Matrix(prod(rConstitutiveMatrix, L)));
    const SizeType column_start = self_index * TDim * TDim;
    noalias(project(rK12, range(0, rK12.size1()), range(column_start, column_start + TDim * TDim))) += material_tangent;

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateLinearMomentumResidualVector(
    VectorType& rRHSv,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rProcessInfo,
    const VectorType& rStressVector,
    const double weight,
    const int Step
)
{
    KRATOS_TRY

    const SizeType number_of_neigh = this->GetValue(NEIGHBOURS).size();
    
    // Only external forces contributions are taken into account at the moment 
    VectorType body_force(TDim); 
    SPHElementUtilities::GetLocalBodyForces(*this, body_force, Step);

    for (IndexType i = 0; i < number_of_neigh; ++i){
        const SizeType index = i * TDim;
        for (IndexType d = 0; d < TDim; ++d)
            rRHSv[index + d] += weight * rThisKinematicVariables.W[i] * body_force[d];
    }

    // Adding stress contribution to the residual vector
    noalias(rRHSv) -= weight * prod(trans(rThisKinematicVariables.B), rStressVector);

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateGeometricalTangentMatrix(
    MatrixType& rK21,
    const KinematicVariables& rThisKinematicVariables,
    const double weight
)
{
    KRATOS_TRY
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neigh = r_neighbours.size();

    const MatrixType E = IdentityMatrix(TDim);
    VectorType e(TDim), temp_vec(TDim * TDim);
    MatrixType temp(TDim, TDim);

    for (IndexType i = 0; i < number_of_neigh; ++i){
        for (IndexType d = 0; d < TDim; ++d){
            
            e = column(E, d); 
            temp = outer_prod(e, row(rThisKinematicVariables.DW_DX, i));
            temp_vec = SPHElementUtilities<TDim>::NonSymmetricTensorToVector(temp);
            
            const SizeType col_index = i * TDim + d;
            
            for (IndexType ii = 0; ii < number_of_neigh; ++ii){
                
                const SizeType row_start = TDim * TDim * ii;
                
                for (IndexType k = 0; k < TDim * TDim; ++k){
                    rK21(row_start + k, col_index) += weight * rThisKinematicVariables.W[ii] * temp_vec[k];
                }
            }
        }
    }

    KRATOS_CATCH("")   
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateGeometricalResidualVector(
    VectorType& rRHSF,
    KinematicVariables& rThisKinematicVariables,
    const double weight,
    const int Step
)
{
    KRATOS_TRY

    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neigh = r_neighbours.size();

    VectorType temp_residual(TDim * TDim), vel_aux(TDim);
    MatrixType temp = ZeroMatrix(TDim, TDim); 

    for (IndexType i = 0; i < number_of_neigh; ++i){
        
        const array_1d<double, 3>& velocity = r_neighbours[i]->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, Step);
        for (IndexType d = 0; d < TDim; ++d) vel_aux[d] = velocity[d];

        temp += outer_prod(vel_aux, row(rThisKinematicVariables.DW_DX, i));
    }

    temp_residual = SPHElementUtilities<TDim>::NonSymmetricTensorToVector(temp); 

    for (IndexType i = 0; i < number_of_neigh; ++i){
        noalias(project(rRHSF, range(TDim * TDim * i, TDim * TDim * (i + 1)))) += temp_residual * rThisKinematicVariables.W[i];
    }

    rRHSF *= weight;

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::AssembleLHS(
    MatrixType& rLHS,
    const MatrixType& rK11,
    const MatrixType& rK12,
    const MatrixType& rK21,
    const MatrixType& rK22
    )
{
    KRATOS_TRY
    const SizeType number_of_neighbours = this->GetValue(NEIGHBOURS).size();
    const SizeType first_equation_dofs = number_of_neighbours * TDim;

    for (IndexType i = 0; i < rK11.size1(); ++i)
        for (IndexType j = 0; j < rK11.size2(); ++j)
            rLHS(i, j) = rK11(i, j);
    
    for (IndexType i = 0; i < rK12.size1(); ++i)
        for (IndexType j = 0; j < rK12.size2(); ++j)
            rLHS(i, first_equation_dofs + j) = rK12(i, j);

    for (IndexType i = 0; i < rK21.size1(); ++i)
        for (IndexType j = 0; j < rK21.size2(); ++j)
            rLHS(first_equation_dofs + i, j) = rK21(i, j);

    for (IndexType i = 0; i < rK22.size1(); ++i)
        for (IndexType j = 0; j < rK22.size2(); ++j)
            rLHS(first_equation_dofs + i, first_equation_dofs + j) = rK22(i, j);

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::AssembleRHS(
    VectorType& rRHS,
    const VectorType& rRHSv,
    const VectorType& rRHSF
    )
{
    KRATOS_TRY
    const SizeType number_of_neighbours = this->GetValue(NEIGHBOURS).size();
    const SizeType first_equation_dofs = number_of_neighbours * TDim;
    
    for (IndexType i = 0; i < rRHSv.size(); ++i)
        rRHS[i] = rRHSv[i];
    
    for (IndexType i = 0; i < rRHSF.size(); ++i)
        rRHS[first_equation_dofs + i] = rRHSF[i];

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY
    const auto& r_geom = this->GetGeometry();
    const auto& r_prop = this->GetProperties();
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    const SizeType block_size = TDim + TDim * TDim;
    const SizeType mat_size = block_size * number_of_neighbours;
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize(mat_size, mat_size, false);
    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size); 

    KRATOS_ERROR_IF(r_prop.Has(DENSITY) == false) << "DENSITY not provided for element " << this->Id() << std::endl;

    const double density = r_prop[DENSITY];
    double thickness = 1.0;
    if constexpr (TDim == 2) thickness = r_prop.Has(THICKNESS) ? r_prop[THICKNESS] : 1.0;

    const auto& GPcoords = r_geom[0].GetInitialPosition();
    const double gauss_weight = r_geom[0].GetValue(VOLUME);

    double factor = density * thickness * gauss_weight;
    bool compute_lumped_mass_matrix = SPHElementUtilities::ComputeLumpedMassMatrix(r_prop, rProcessInfo);

    if (compute_lumped_mass_matrix){
        int this_id = this->Id();
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            if (r_neighbours[i]->Id() == this_id){
                const SizeType v_block = i * TDim;
                const SizeType F_block = i * TDim * TDim + TDim * number_of_neighbours;
                
                for (IndexType d = 0; d < TDim; ++d)
                    rMassMatrix(v_block + d, v_block + d) = factor;

                for (IndexType d = 0; d < TDim * TDim; ++d)
                    rMassMatrix(F_block + d, F_block + d) = gauss_weight;
            }
        }
    } else { // Consistent mass matrix
        const SizeType v_size = TDim * number_of_neighbours;
        const SizeType F_size = TDim * TDim * number_of_neighbours;

        MatrixType vBlock(v_size, TDim); vBlock.clear();
        MatrixType FBlock(F_size, TDim * TDim); FBlock.clear();
        
        VectorType X_AB_target(TDim);
        double kernel;

        for (IndexType i = 0; i < number_of_neighbours; ++i){

            const auto& JPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
            const double volume = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);

            for (IndexType d = 0; d < TDim; d++)
                X_AB_target[d] = GPcoords[d] - JPcoords[d];

            TKernelType::ComputeKernelValue(kernel, h, X_AB_target); 
            ComputeKernelCorrectionUtilities::ApplyKernelCorrection(*this, kernel);

            for (IndexType d = 0; d < TDim; ++d)
                vBlock(TDim * i + d, d) = kernel * volume;

            for (IndexType d = 0; d < TDim * TDim; ++d)
                FBlock(TDim * TDim * i + d, d) = kernel * volume;
        
        }

        noalias(project(rMassMatrix, range(0, v_size), range(0, v_size))) += factor * prod(vBlock, trans(vBlock));
        noalias(project(rMassMatrix, range(v_size, v_size + F_size), range(v_size, v_size + F_size))) += gauss_weight * prod(FBlock, trans(FBlock));
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateAndAddUpwindStabilizationResidual(
    VectorType& rRHSv,
    KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rProcessInfo, 
    int Step
)
{
    KRATOS_TRY
    const auto& r_geom = this->GetGeometry();
    const auto& r_props = this->GetProperties();
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);

    VectorType upwind_residual(TDim), velocity_jump(TDim); upwind_residual.clear();
    MatrixType stabilization_matrix(TDim, TDim);

    const int self_index = this->GetNeighbourPosition(r_neighbours);
    const array_1d<double, 3>& velocity_self = r_geom[0].FastGetSolutionStepValue(VELOCITY, Step);

    for (IndexType i = 0; i < r_neighbours.size(); ++i){
        
        if (i == self_index) continue;

        const array_1d<double, 3>& velocity_neigh = r_neighbours[i]->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, Step);

        for (IndexType d = 0; d < TDim; ++d)
            velocity_jump[d] = velocity_self[d] - velocity_neigh[d];
        
        CalculatePairUpwindStabilizationMatrix(stabilization_matrix, *r_neighbours[i], rThisKinematicVariables.F, rProcessInfo);

        upwind_residual += prod(stabilization_matrix, velocity_jump);
    }

    noalias(project(rRHSv, range(TDim * self_index, TDim * (self_index + 1)))) += upwind_residual;

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateAndAddUpwindStabilizationTangent(
    MatrixType& rK11,
    KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rProcessInfo
)
{
    KRATOS_TRY
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);

    MatrixType stabilization_matrix(TDim, TDim);

    const int self_index = this->GetNeighbourPosition(r_neighbours);

    for (IndexType i = 0; i < r_neighbours.size(); ++i){
        
        if (i == self_index) continue;
        
        CalculatePairUpwindStabilizationMatrix(stabilization_matrix, *r_neighbours[i], rThisKinematicVariables.F, rProcessInfo);
    
        // For simplicity, the tanget matrix doesn't consider the dependence of pressure and shear wave speeds on the deformation gradient 
        noalias(project(rK11, range(TDim * self_index, TDim * (self_index + 1)), range(TDim * i, TDim * (i + 1)))) += stabilization_matrix;
        noalias(project(rK11, range(TDim * self_index, TDim * (self_index + 1)), range(TDim * self_index, TDim * (self_index + 1)))) -= stabilization_matrix;
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculatePairUpwindStabilizationMatrix(
    MatrixType& rStabilizationMatrix,
    Element& rNeighbour,
    const MatrixType& rThisDeformationGradient,
    const ProcessInfo& rProcessInfo
)
{
    rStabilizationMatrix.clear();

    const auto& r_geom = this->GetGeometry();
    const auto& r_props = this->GetProperties();
    const auto& r_geom_neigh = rNeighbour.GetGeometry();

    VectorType X_AB_target(TDim), normal(TDim);
    for (IndexType d = 0; d < TDim; ++d)
        X_AB_target[d] = r_geom[0].GetInitialPosition()[d] - r_geom_neigh[0].GetInitialPosition()[d];
    
    const double norm_dist = norm_2(X_AB_target);

    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);
    double kernel;
    VectorType dkernel_b_at_a(TDim);

    TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
    TKernelType::ComputeKernelGradientValue(dkernel_b_at_a, h, X_AB_target);

    VectorType dkernel_a_at_b = dkernel_b_at_a;
    ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, kernel, dkernel_b_at_a);
    ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(rNeighbour, kernel, dkernel_a_at_b);

    const VectorType skew_kernel_gradient = 0.5 * (dkernel_b_at_a - dkernel_a_at_b);

    double pressure_wave_speed, shear_wave_speed;
    SPHElementUtilities::ComputeDeformationDependentWaveSpeed(pressure_wave_speed, shear_wave_speed, rThisDeformationGradient, r_props);

    const double density = r_props[DENSITY];
    const double volume_product = r_geom[0].GetValue(VOLUME) * r_geom_neigh[0].GetValue(VOLUME);
    const double coeff = rProcessInfo.GetValue(DISSIPATION_COEFFICIENT);

    normal = X_AB_target / norm_dist;
    rStabilizationMatrix = coeff * volume_product * density * norm_2(skew_kernel_gradient) *
        (pressure_wave_speed * outer_prod(normal, normal) + shear_wave_speed * (IdentityMatrix(TDim) - outer_prod(normal, normal)));

}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rProcessInfo,
    const int Step
    )
{
    KRATOS_TRY

    CalculateKernelsAndKernelGradients(rThisKinematicVariables.DW_DX, rThisKinematicVariables.W, rProcessInfo);
    AssembleDeformationGradient(rThisKinematicVariables.F, Step);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);

    const SizeType number_of_neighbours = this->GetValue(NEIGHBOURS).size();

    if constexpr (TDim == 2){
        SPHElementUtilities::Calculate2DB(rThisKinematicVariables.B, rThisKinematicVariables.F, rThisKinematicVariables.DW_DX, number_of_neighbours);
    } else {
        SPHElementUtilities::Calculate3DB(rThisKinematicVariables.B, rThisKinematicVariables.F, rThisKinematicVariables.DW_DX, number_of_neighbours);
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::CalculateKernelsAndKernelGradients(
    MatrixType& rDW_DX,
    VectorType& rW,
    const ProcessInfo& rProcessInfo
    )
{
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);
    const auto& IPcoords = this->GetGeometry()[0].GetInitialPosition();

    // Initialization of variables
    double kernel;
    VectorType dkernel(TDim), X_AB_target(TDim);

    for (IndexType i = 0; i < r_neighbours.size(); ++i){
        
        const auto& JPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
        const auto& JP_current_position = r_neighbours[i]->GetGeometry()[0].Coordinates();
        
        const double weight = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);

        for (IndexType d = 0; d < TDim; d++) 
            X_AB_target[d] = IPcoords[d] - JPcoords[d];

        TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
        TKernelType::ComputeKernelGradientValue(dkernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, kernel, dkernel);
            
        for (IndexType d = 0; d < TDim; d++) 
            rDW_DX(i, d) = weight * dkernel[d];

        rW[i] = weight * kernel;
    }
}

template<class TKernelType, std::size_t TDim>
void TotalLagrangianMixedStrainParticle<TKernelType, TDim>::AssembleDeformationGradient(
    MatrixType& rF,
    const int Step
)
{
    const auto& r_geom = this->GetGeometry();
    
    if constexpr (TDim == 2){
        rF(0,0) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XX, Step);
        rF(1,1) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YY, Step);
        rF(0,1) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XY, Step);
        rF(1,0) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YX, Step);
    } else {
        rF(0,0) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XX, Step);
        rF(1,1) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YY, Step);
        rF(2,2) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_ZZ, Step);
        rF(0,1) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XY, Step);
        rF(0,2) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XZ, Step);
        rF(1,0) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YX, Step);
        rF(1,2) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YZ, Step);
        rF(2,0) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_ZX, Step);
        rF(2,1) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_ZY, Step);
    }
}
    
template class TotalLagrangianMixedStrainParticle<CubicKernel2D, 2>;
template class TotalLagrangianMixedStrainParticle<CubicKernel3D, 3>;

}
