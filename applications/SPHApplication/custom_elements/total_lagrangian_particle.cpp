#include "custom_elements/total_lagrangian_particle.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::InitializeSolutionStep(const ProcessInfo& rProcessInfo)
{
    bool required = false;
    if (this->mThisConstitutiveLaw->RequiresInitializeMaterialResponse()) required = true;

    if (required){
        const auto& r_geom = this->GetGeometry();
        const auto& r_prop = this->GetProperties();
        const SizeType number_of_neighbours = this->GetValue(NEIGHBOURS).size();
        const SizeType dimension = r_geom.WorkingSpaceDimension();
        const SizeType strain_size = this->mThisConstitutiveLaw->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_neighbours);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        ConstitutiveLaw::Parameters Values(r_geom, r_prop, rProcessInfo);

        auto& ConstitutiveLawOptions = Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        Values.SetStrainVector(this_constitutive_variables.StrainVector);

        this->CalculateKinematicVariables(this_kinematic_variables, rProcessInfo);

        this->SetConstitutiveLawVariables(this_constitutive_variables, this_kinematic_variables, Values);

        this->mThisConstitutiveLaw->InitializeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);
    }
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::FinalizeSolutionStep(const ProcessInfo& rProcessInfo)
{
    bool required = false;
    if (this->mThisConstitutiveLaw->RequiresFinalizeMaterialResponse()) required = true;

    if (required) {
        auto& r_geom = this->GetGeometry();
        const auto& r_prop = this->GetProperties();
        const SizeType number_of_neighbours = this->GetValue(NEIGHBOURS).size();
        const SizeType dimension = r_geom.WorkingSpaceDimension();
        const SizeType strain_size = this->mThisConstitutiveLaw->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_neighbours);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        ConstitutiveLaw::Parameters Values(r_geom, r_prop, rProcessInfo);

        auto& ConstitutiveLawOptions = Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        Values.SetStrainVector(this_constitutive_variables.StrainVector);

        this->CalculateKinematicVariables(this_kinematic_variables, rProcessInfo);
        
        this->SetConstitutiveLawVariables(this_constitutive_variables, this_kinematic_variables, Values);

        this->mThisConstitutiveLaw->FinalizeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

        // This assign the value of PLASTIC_STRAIN_VECTOR to a non historical variable of the geometry
        VectorType temp(strain_size);
        this->mThisConstitutiveLaw->GetValue(PLASTIC_STRAIN_VECTOR, temp);
        r_geom[0].SetValue(PLASTIC_STRAIN_VECTOR, temp);
    }
}

template<class TKernelType>
Element::Pointer TotalLagrangianDisplacementParticle<TKernelType>::Clone( 
    IndexType NewId, 
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    TotalLagrangianDisplacementParticle<TKernelType>::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangianDisplacementParticle<TKernelType>>
        (NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    p_new_elem->SetConstitutiveLaw(BaseType::mThisConstitutiveLaw);

    return p_new_elem;

    KRATOS_CATCH("");
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateAll(
    MatrixType& rLHS, 
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY
    const auto&  r_geom = this->GetGeometry();
    const auto& r_props = this->GetProperties();
    const SizeType domain_size = r_geom.WorkingSpaceDimension();
    const SizeType strain_size = this->mThisConstitutiveLaw->GetStrainSize();

    const SizeType number_of_neigh = this->GetValue(NEIGHBOURS).size();
    const SizeType mat_size = domain_size * number_of_neigh;

    KinematicVariables this_kinematic_variables(strain_size, domain_size, number_of_neigh);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    if (CalculateStiffnessMatrixFlag){
        if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) rLHS.resize(mat_size, mat_size, false);
        noalias(rLHS) = ZeroMatrix(mat_size, mat_size);
    }

    if (CalculateResidualVectorFlag){
        if (rRHS.size() != mat_size) rRHS.resize(mat_size, false);
        noalias(rRHS) = ZeroVector(mat_size);
    }

    ConstitutiveLaw::Parameters Values(r_geom, r_props, rProcessInfo);
    auto& ConstitutiveLawOptions = Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if (CalculateStiffnessMatrixFlag) {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false); 
    }
    Values.SetStrainVector(this_constitutive_variables.StrainVector);

    double gauss_weight = r_geom[0].GetValue(VOLUME);
    VectorType body_force(domain_size);
    SPHElementUtilities::GetLocalBodyForces(*this, body_force);

    this->CalculateKinematicVariables(this_kinematic_variables, rProcessInfo);
    this->CalculateConstitutiveVariables(this_constitutive_variables, this_kinematic_variables, Values, ConstitutiveLaw::StressMeasure_PK2);

    const double thickness = (domain_size == 2 && r_props.Has(THICKNESS)) ? r_props[THICKNESS] : 1.0;
    gauss_weight *= thickness;

    if (CalculateStiffnessMatrixFlag){
        /* Geometric stiffness matrix */
        CalculateAndAddKg(rLHS, this_kinematic_variables.DW_DX, this_constitutive_variables.StressVector, gauss_weight);

        /* Material stiffness matrix */
        CalculateAndAddKm(rLHS, this_kinematic_variables.B, this_constitutive_variables.C, gauss_weight);  
    }

    if (CalculateResidualVectorFlag){
        CalculateAndAddResidualVector(rRHS, this_kinematic_variables, rProcessInfo, body_force, this_constitutive_variables.StressVector, gauss_weight);
    }

    if (rProcessInfo[PENALIZATION_COEFFICIENT] != 0.0){
        CalculateAndAddPenalization(rLHS, rRHS, this_kinematic_variables, rProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }
        
    
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateKinematicVariables(KinematicVariables& rThisKinematicVariables, const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    CalculateDeformationGradient(rThisKinematicVariables.F, rThisKinematicVariables.DW_DX, rThisKinematicVariables.W, rProcessInfo);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);

    if (this->GetGeometry().WorkingSpaceDimension() == 2){
        Calculate2DB(rThisKinematicVariables.B, rThisKinematicVariables.F, rThisKinematicVariables.DW_DX);
    } else {
        Calculate3DB(rThisKinematicVariables.B, rThisKinematicVariables.F, rThisKinematicVariables.DW_DX);
    }

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateDeformationGradient(
    MatrixType& rF, 
    MatrixType& rDW_DX,
    VectorType& rW,
    const ProcessInfo& rProcessInfo)
{
    // Initialization of variables
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    const int domain_size = this->GetGeometry().WorkingSpaceDimension();
    const auto& IPcoords = this->GetGeometry()[0].GetInitialPosition();

    // Initialization of variables
    double kernel;
    VectorType dkernel(domain_size), X_AB_target(domain_size), nodal_values(domain_size);

    rF = ZeroMatrix(domain_size, domain_size);
    const int self_index = GetNeighbourPosition(r_neighbours);

    for (IndexType i = 0; i < r_neighbours.size(); ++i){
        
        const auto& JPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
        const auto& JP_current_position = r_neighbours[i]->GetGeometry()[0].Coordinates();
        
        const double weight = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);

        for (IndexType d = 0; d < domain_size; d++) X_AB_target[d] = IPcoords[d] - JPcoords[d];

        TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
        TKernelType::ComputeKernelGradientValue(dkernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, kernel, dkernel);
        ComputeKernelCorrectionUtilities::ApplyIntegrationCorrection(*this, kernel, dkernel, i == self_index);
            
        for (IndexType d = 0; d < domain_size; d++){
            nodal_values[d] = JP_current_position[d];
            rDW_DX(i, d) = weight * dkernel[d];
        }

        rF += weight * outer_prod(nodal_values, dkernel);
        rW[i] = weight * kernel;
    }
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::Calculate2DB(MatrixType& rB, const MatrixType& rF, const MatrixType& rDW_DX)
{
    const SizeType number_of_neigh = this->GetValue(NEIGHBOURS).size();
    const SizeType domain_size = this->GetGeometry().WorkingSpaceDimension();

    for (IndexType i =0; i < number_of_neigh; ++i){
        const IndexType index = i * domain_size;
        rB(0, index + 0) = rF(0, 0) * rDW_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDW_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDW_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDW_DX(i, 1);
        rB(2, index + 0) = rF(0, 0) * rDW_DX(i, 1) + rF(0, 1) * rDW_DX(i, 0);
        rB(2, index + 1) = rF(1, 0) * rDW_DX(i, 1) + rF(1, 1) * rDW_DX(i, 0);
    }

}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::Calculate3DB(MatrixType& rB, const MatrixType& rF, const MatrixType& rDW_DX)
{
    const SizeType number_of_neigh = this->GetValue(NEIGHBOURS).size();
    const SizeType domain_size = this->GetGeometry().WorkingSpaceDimension();

    for (IndexType i =0; i < number_of_neigh; ++i){
        const IndexType index = i * domain_size;
        rB(0, index + 0) = rF(0, 0) * rDW_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDW_DX(i, 0);
        rB(0, index + 2) = rF(2, 0) * rDW_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDW_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDW_DX(i, 1);
        rB(1, index + 2) = rF(2, 1) * rDW_DX(i, 1);
        rB(2, index + 0) = rF(0, 2) * rDW_DX(i, 2);
        rB(2, index + 1) = rF(1, 2) * rDW_DX(i, 2);
        rB(2, index + 2) = rF(2, 2) * rDW_DX(i, 2);
        rB(3, index + 0) = rF(0, 0) * rDW_DX(i, 1) + rF(0, 1) * rDW_DX(i, 0);
        rB(3, index + 1) = rF(1, 0) * rDW_DX(i, 1) + rF(1, 1) * rDW_DX(i, 0);
        rB(3, index + 2) = rF(2, 0) * rDW_DX(i, 1) + rF(2, 1) * rDW_DX(i, 0);
        rB(4, index + 0) = rF(0, 1) * rDW_DX(i, 2) + rF(0, 2) * rDW_DX(i, 1);
        rB(4, index + 1) = rF(1, 1) * rDW_DX(i, 2) + rF(1, 2) * rDW_DX(i, 1);
        rB(4, index + 2) = rF(2, 1) * rDW_DX(i, 2) + rF(2, 2) * rDW_DX(i, 1);
        rB(5, index + 0) = rF(0, 2) * rDW_DX(i, 0) + rF(0, 0) * rDW_DX(i, 2);
        rB(5, index + 1) = rF(1, 2) * rDW_DX(i, 0) + rF(1, 0) * rDW_DX(i, 2);
        rB(5, index + 2) = rF(2, 2) * rDW_DX(i, 0) + rF(2, 0) * rDW_DX(i, 2);
    }

}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateConstitutiveVariables(
    ConstitutiveVariables& rThisConstitutiveVariables,
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
)
{
    SetConstitutiveLawVariables(rThisConstitutiveVariables, rThisKinematicVariables, rValues);
        
    this->mThisConstitutiveLaw->CalculateMaterialResponse(rValues, ThisStressMeasure);
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::SetConstitutiveLawVariables(
    ConstitutiveVariables& rThisConstitutiveVariables,
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveLaw::Parameters& rValues
)
{
    // Essential input parameters for the constitutive law
    rValues.SetDeterminantF(rThisKinematicVariables.detF);
    rValues.SetDeformationGradientF(rThisKinematicVariables.F);

    // Space in which the resulta shall be saved
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.C);
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector);
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateAndAddKg(
    MatrixType& rLHS, 
    const Matrix& DW_DX, 
    const Vector& stress_vector, 
    const double weight
    ) const 
{
    KRATOS_TRY
    
    const SizeType domain_size = this->GetGeometry().WorkingSpaceDimension();
    const MatrixType stress_tensor = weight * MathUtils<double>::StressVectorToTensor(stress_vector);
    MatrixType reduced_Kg(DW_DX.size1(), DW_DX.size1());
    MathUtils<double>::BDBtProductOperation(reduced_Kg, stress_tensor, DW_DX);
    MathUtils<double>::ExpandAndAddReducedMatrix(rLHS, reduced_Kg, domain_size);
    
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateAndAddKm(
    MatrixType& rLHS, 
    const Matrix& rB,
    const Matrix& rConstitutiveMatrix, 
    const double weight
    ) const 
{
    KRATOS_TRY

    noalias(rLHS) += weight * prod(trans(rB), Matrix(prod(rConstitutiveMatrix, rB)));

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateAndAddResidualVector(
    VectorType& rRHS,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rProcessInfo,
    const VectorType& rBodyForce,
    const VectorType& rStressVector,
    const double weigth 
    ) const
{
    KRATOS_TRY

    this->CalculateAndAddExternalForcesContribution(rThisKinematicVariables.W, rProcessInfo, rBodyForce, rRHS, weigth);
    noalias(rRHS) -= weigth * prod(trans(rThisKinematicVariables.B), rStressVector);

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateAndAddPenalization(
    MatrixType& rLHS,
    VectorType& rRHS,
    KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY
    
    const auto& r_props = this->GetProperties();
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neigh = r_neighbours.size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    // Initialization of variables
    double norm_dist, kernel, norm_dkernel, penalization_factor;
    VectorType X_AB_target(dimension), normal(dimension), dkernel(dimension), displacement_jump(dimension);
    MatrixType AcousticTensor(dimension, dimension);

    const int self_index = GetNeighbourPosition(r_neighbours);

    const double alpha = rProcessInfo.GetValue(PENALIZATION_COEFFICIENT);
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    const auto& IPcoords = this->GetGeometry()[0].GetInitialPosition();
    const double weight1 = this->GetGeometry()[0].GetValue(VOLUME);

    for (IndexType i = 0; i < number_of_neigh; ++i){

        const double weight2 = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);
        const auto& JPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();

        for (IndexType d = 0; d < dimension; ++d) X_AB_target[d] = IPcoords[d] - JPcoords[d];

        TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
        TKernelType::ComputeKernelGradientValue(dkernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, kernel, dkernel);

        norm_dist = norm_2(X_AB_target);
        if (norm_dist < 1e-15) continue; 

        normal = X_AB_target / norm_dist; 
        
        penalization_factor = alpha * weight1 * weight2 * norm_2(dkernel) / norm_dist;

        SPHElementUtilities::ComputeLinearElasticAcousticTensor(AcousticTensor, normal, r_props);

        if (CalculateStiffnessMatrixFlag){
            noalias(project(rLHS, range(self_index * dimension, (self_index + 1) * dimension), range(i * dimension, (i + 1) * dimension))) -= penalization_factor * AcousticTensor;
            noalias(project(rLHS, range(self_index * dimension, (self_index + 1) * dimension), range(self_index * dimension, (self_index + 1) * dimension))) += penalization_factor * AcousticTensor;
        }

        if (CalculateResidualVectorFlag){
            SPHElementUtilities::ComputeParticleJump(displacement_jump, *this, *r_neighbours[i], X_AB_target, rProcessInfo); 
            noalias(project(rRHS, range(self_index * dimension, (self_index + 1) * dimension))) += penalization_factor * prod(AcousticTensor, displacement_jump);
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rProcessInfo)
{
    BaseType::CalculateDampingMatrix(rDampingMatrix, rProcessInfo);
    if (rProcessInfo[DISSIPATION_COEFFICIENT] != 0.0) CalculateAndAddDissipation(rDampingMatrix, rProcessInfo);

}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateAndAddDissipation(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rProcessInfo
)
{
    KRATOS_TRY
    const auto& r_props = this->GetProperties();
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neigh = r_neighbours.size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    // Initialization of variables
    double pressure_wave_speed, shear_wave_speed, norm_dist, kernel, 
        norm_dkernel, dissipation_factor;
    VectorType X_AB_target(dimension), normal(dimension), dkernel(dimension), jump(dimension);
    MatrixType DissipationMatrix(dimension, dimension);

    const int self_index = GetNeighbourPosition(r_neighbours);

    const double alpha = rProcessInfo.GetValue(DISSIPATION_COEFFICIENT);
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    const auto& IPcoords = this->GetGeometry()[0].GetInitialPosition();
    const double weight1 = this->GetGeometry()[0].GetValue(VOLUME);
    const double density = r_props[DENSITY];

    SPHElementUtilities::ComputeWaveSpeed(pressure_wave_speed, shear_wave_speed, r_props);

    for (IndexType i = 0; i < number_of_neigh; ++i){

        const double weight2 = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);
        const auto& JPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();

        for (IndexType d = 0; d < dimension; ++d) X_AB_target[d] = IPcoords[d] - JPcoords[d];

        TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
        TKernelType::ComputeKernelGradientValue(dkernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, kernel, dkernel);

        norm_dist = norm_2(X_AB_target);
        if (norm_dist < 1e-15) continue; 

        normal = X_AB_target / norm_dist; 
        
        dissipation_factor  = alpha * density * weight1 * weight2 * norm_2(dkernel);

        DissipationMatrix = pressure_wave_speed * outer_prod(normal, normal) + shear_wave_speed * (IdentityMatrix(dimension) - outer_prod(normal, normal));

        noalias(project(rDampingMatrix, range(self_index * dimension, (self_index + 1) * dimension), range(i * dimension, (i + 1) * dimension))) -= dissipation_factor * DissipationMatrix;
        noalias(project(rDampingMatrix, range(self_index * dimension, (self_index + 1) * dimension), range(self_index * dimension, (self_index + 1) * dimension))) += dissipation_factor * DissipationMatrix;
        
    }

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rOutput,
    const ProcessInfo& rProcessInfo
)
{
    if (rVariable == F_DEFORMATION_GRADIENT){
        // Initialization of variables
        const auto& r_neighbours = this->GetValue(NEIGHBOURS);
        const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

        const int domain_size = this->GetGeometry().WorkingSpaceDimension();
        const auto& IPcoords = this->GetGeometry()[0].GetInitialPosition();

        // Initialization of variables
        double kernel;
        VectorType dkernel(domain_size), X_AB_target(domain_size), nodal_values(domain_size);

        MatrixType rF = ZeroMatrix(domain_size, domain_size);
        const int self_index = GetNeighbourPosition(r_neighbours);

        for (IndexType i = 0; i < r_neighbours.size(); ++i){

            const auto& JPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
            const auto& JP_current_position = r_neighbours[i]->GetGeometry()[0].Coordinates();

            const double weight = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);

            for (IndexType d = 0; d < domain_size; d++) X_AB_target[d] = IPcoords[d] - JPcoords[d];

            TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
            TKernelType::ComputeKernelGradientValue(dkernel, h, X_AB_target);
            ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, kernel, dkernel);
            ComputeKernelCorrectionUtilities::ApplyIntegrationCorrection(*this, kernel, dkernel, i == self_index);

            for (IndexType d = 0; d < domain_size; d++){
                nodal_values[d] = JP_current_position[d];
            }

            rF += weight * outer_prod(nodal_values, dkernel);
        }

        if (rOutput.size() != 1)
        rOutput.resize(1);

        rOutput[0] = rF;
    }
}

template class TotalLagrangianDisplacementParticle<CubicKernel2D>;
template class TotalLagrangianDisplacementParticle<CubicKernel3D>;

}

