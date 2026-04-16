#include "custom_elements/total_lagrangian_particle.h"

namespace Kratos
{

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    InitializeMaterial();

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::InitializeMaterial()
{
    KRATOS_TRY

    if (GetProperties()[CONSTITUTIVE_LAW] != nullptr){
        const auto& r_geom = GetGeometry();
        const auto& r_prop = GetProperties();

        mThisConstitutiveLaw = r_prop[CONSTITUTIVE_LAW]->Clone();
        mThisConstitutiveLaw->InitializeMaterial(r_prop, r_geom, Vector());

    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID" << this->Id() << std::endl;
    }

    KRATOS_CATCH("")
}

template<class TKernelType>
Element::Pointer TotalLagrangianDisplacementParticle<TKernelType>::Clone( 
    IndexType NewId, 
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    TotalLagrangianDisplacementParticle<TKernelType>::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangianDisplacementParticle<TKernelType>>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    p_new_elem->SetConstitutiveLaw(mThisConstitutiveLaw);

    return p_new_elem;

    KRATOS_CATCH("");
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::GetNodalValuesVector(VectorType& rNodalValues) const
{
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    if (rNodalValues.size() != dimension * number_of_neighbours){
        rNodalValues.resize(dimension * number_of_neighbours, false);
    }

    IndexType index = 0;

    for (IndexType i = 0; i < number_of_neighbours; ++i){

        const auto& r_geom = r_neighbours[i]->GetGeometry();
        const auto& r_displ = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT);
        rNodalValues[index++] = r_displ[0];
        rNodalValues[index++] = r_displ[1];
        if (dimension == 3){
            rNodalValues[index++] = r_displ[2];
        }
    }
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    IndexType local_index = 0;

    if (rResult.size() != dimension * number_of_neighbours)
        rResult.resize(dimension * number_of_neighbours, false);

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        const auto& r_geom = r_neighbours[i]->GetGeometry();

        const IndexType xpos = r_geom[0].GetDofPosition(DISPLACEMENT_X);
        
        rResult[local_index++] = r_geom[0].GetDof(DISPLACEMENT_X, xpos).EquationId();
        rResult[local_index++] = r_geom[0].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
        if (dimension == 3){
            rResult[local_index++] = r_geom[0].GetDof(DISPLACEMENT_Z, xpos + 2).EquationId();
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    IndexType index = 0;

    if(rElementalDofList.size() != dimension * number_of_neighbours)
        rElementalDofList.resize(dimension * number_of_neighbours);

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        const auto& r_geom = r_neighbours[i]->GetGeometry();

        rElementalDofList[index++] = r_geom[0].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index++] = r_geom[0].pGetDof(DISPLACEMENT_Y);
        if (dimension == 3){
            rElementalDofList[index++] = r_geom[0].pGetDof(DISPLACEMENT_Z);
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::GetFirstDerivativesVector(VectorType& rValues, int step) const
{
    KRATOS_TRY
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = dimension * number_of_neighbours;

    if (rValues.size() != mat_size) rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        
        const auto& r_geom = r_neighbours[i]->GetGeometry();

        const array_1d<double, 3>& velocity = r_geom[0].FastGetSolutionStepValue(VELOCITY, step);
        const SizeType index = i * dimension;
        for (unsigned int k = 0; k < dimension; ++k){
            rValues[index + k] = velocity[k];
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::GetSecondDerivativesVector(VectorType& rValues, int step) const
{
    KRATOS_TRY
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = dimension * number_of_neighbours;

    if (rValues.size() != mat_size) rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        
        const auto& r_geom = r_neighbours[i]->GetGeometry();

        const array_1d<double, 3>& acceleration = r_geom[0].FastGetSolutionStepValue(ACCELERATION, step);
        const SizeType index = i * dimension;
        for (unsigned int k = 0; k < dimension; ++k){
            rValues[index + k] = acceleration[k];
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
)
{
    KRATOS_TRY
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = true;
    
    CalculateAll(rLHS, rRHS, rProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
)
{
    KRATOS_TRY
    
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = false;
    VectorType RHS;
    
    CalculateAll(rLHS, RHS, rProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
)
{
    KRATOS_TRY
    const bool CalculateStiffnessMatrixFlag = false;
    const bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();
     
    CalculateAll(temp, rRHS, rProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    
    KRATOS_CATCH("")
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
    const auto&  r_geom = GetGeometry();
    const auto& r_props = GetProperties();
    const SizeType domain_size = r_geom.WorkingSpaceDimension();
    const SizeType strain_size = mThisConstitutiveLaw->GetStrainSize();
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    // The neighbours of the particle are like the nodes belonging to the element in standard FEM
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neigh = r_neighbours.size();
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

    VectorType nodal_values(mat_size), dkernel(domain_size), X_AB_target(domain_size), body_force(domain_size);
    double kernel;
    GetNodalValuesVector(nodal_values);

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

    const auto& GPcoords = r_geom[0].GetInitialPosition();
    double gauss_weight = r_geom[0].GetValue(VOLUME);
    SPHElementUtilities::GetLocalBodyForces(*this, body_force);

    CalculateKinematicVariables(this_kinematic_variables, rProcessInfo);
    CalculateConstitutiveVariables(this_constitutive_variables, this_kinematic_variables, Values, ConstitutiveLaw::StressMeasure_PK2);

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
    
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateKinematicVariables(KinematicVariables& rThisKinematicVariables, const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    CalculateDeformationGradient(rThisKinematicVariables.F, rThisKinematicVariables.DW_DX, rThisKinematicVariables.W, rProcessInfo);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
    
    CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.F, rThisKinematicVariables.DW_DX);

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateDeformationGradient(
    MatrixType& rF, 
    MatrixType& rDW_DX,
    VectorType rW,
    const ProcessInfo& rProcessInfo)
{
    // Initialization of variables
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    const int domain_size = GetGeometry().WorkingSpaceDimension();
    const auto& IPcoords = GetGeometry()[0].GetInitialPosition();

    // Initialization of variables
    double kernel;
    VectorType dkernel(domain_size), X_AB_target(domain_size), 
    nodal_values(domain_size * r_neighbours.size()), local_nodal_values(domain_size);

    GetNodalValuesVector(nodal_values);
    rF = IdentityMatrix(domain_size);

    for (IndexType i = 0; i < r_neighbours.size(); ++i){
        const auto& JPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
        const double weight = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);

        for (IndexType d = 0; d < domain_size; d++) X_AB_target[d] = IPcoords[d] - JPcoords[d];

        TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
        TKernelType::ComputeKernelGradientValue(dkernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, kernel, dkernel);
            
        GlobalSizeVector(local_nodal_values, nodal_values, i);

        rF += weight * outer_prod(local_nodal_values, dkernel);
        for (IndexType d = 0; d < domain_size; d++) rDW_DX(i, d) = weight * dkernel[d];
        rW[i] = weight * kernel;

    }
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateB(MatrixType& rB, const MatrixType& rF, const MatrixType& rDW_DX)
{
    const SizeType number_of_neigh = GetValue(NEIGHBOURS).size();
    const SizeType domain_size = GetGeometry().WorkingSpaceDimension();

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
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateConstitutiveVariables(
    ConstitutiveVariables& rThisConstitutiveVariables,
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
)
{
    SetConstitutiveLawVariables(rThisConstitutiveVariables, rThisKinematicVariables, rValues);

    mThisConstitutiveLaw->CalculateMaterialResponse(rValues, ThisStressMeasure);
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
    
    const SizeType domain_size = GetGeometry().WorkingSpaceDimension();
    const MatrixType stress_tensor = MathUtils<double>::StressVectorToTensor(stress_vector);
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
    MatrixType temp = prod(rConstitutiveMatrix, rB);
    noalias(rLHS) += weight * prod(trans(rB), temp);
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateAndAddResidualVector(
    VectorType& rRHS,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rProcessInfo,
    const VectorType& rBodyForce,
    const Vector& rStressVector,
    const double weigth 
    ) const
{
    KRATOS_TRY

    CalculateAndAddExternelForcesContribution(rThisKinematicVariables.W, rProcessInfo, rBodyForce, rRHS, weigth);
    noalias(rRHS) -= weigth * prod(trans(rThisKinematicVariables.B), rStressVector);

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateAndAddExternelForcesContribution(
    const VectorType& rW,
    const ProcessInfo& rProcessInfo,
    const VectorType& rBodyForce,
    VectorType& rRHS,
    const double weight
    ) const
{
    const SizeType number_of_neigh = GetValue(NEIGHBOURS).size();
    const SizeType domain_size = GetGeometry().WorkingSpaceDimension();

    for (IndexType i = 0; i < number_of_neigh; ++i){
        const SizeType index = i * domain_size;

        for (IndexType d = 0; d < domain_size; ++d){
            rRHS[index + d] += weight * rW[i] * rBodyForce[d];
        }
    }
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neigh = r_neighbours.size();
    
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType mat_size = dimension * number_of_neigh;
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size) rMassMatrix.resize(mat_size, mat_size, false);
    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);

    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY has to be provided for the calculation of the Mass Matrix"<<std::endl;
    
    const double density = r_prop[DENSITY];
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;
    const auto& GPcoords = r_geom[0].GetInitialPosition();
    const double gauss_weight = r_geom[0].GetValue(VOLUME);

    MatrixType J0(dimension, dimension), MassMatrix(mat_size, dimension);
    noalias(MassMatrix) = ZeroMatrix(mat_size, dimension);
    VectorType X_AB_target(dimension);
    double kernel, temp;

    for (IndexType i = 0; i < number_of_neigh; ++i){
        
        const auto& IPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
        const double volume = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);
        
        for (IndexType d = 0; d < dimension; ++d){
            X_AB_target[d] = GPcoords[d] - IPcoords[d];
        }

        TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelCorrection(*this, kernel);

        noalias(J0) = IdentityMatrix(dimension);
        temp = volume * kernel;
        J0 *= temp;

        noalias(project(MassMatrix, range(dimension * i, dimension * (i +1)), range(0, dimension))) += J0;

    }

    double factor = gauss_weight * density * thickness;

    noalias(rMassMatrix) = prod(MassMatrix, trans(MassMatrix));
    rMassMatrix *= factor;

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rProcessInfo)
{
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType mat_size = GetGeometry().WorkingSpaceDimension() * r_neighbours.size();

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(*this, rDampingMatrix, rProcessInfo, mat_size);
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable, 
    std::vector<Vector>& rOutput, 
    const ProcessInfo& rProcessInfo
)
{
    if (rVariable == SPH_KERNEL_GRADIENT){
        const auto& r_neighbours = this->GetValue(NEIGHBOURS);
        const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

        const int domain_size = GetGeometry().WorkingSpaceDimension();
        const auto& IPcoords = GetGeometry()[0].GetInitialPosition();

        IndexType index = 0;
        rOutput.resize(r_neighbours.size());

        for (auto& JP : r_neighbours){
            const auto& JPcoords = JP->GetGeometry()[0].GetInitialPosition();

            VectorType X_AB_target(domain_size);
            for (IndexType d = 0; d < domain_size; d++){
                X_AB_target[d] = IPcoords[d] - JPcoords[d];
            }

            VectorType dkernel(domain_size);
            TKernelType::ComputeKernelGradientValue(dkernel, h, X_AB_target);
            rOutput[index] = dkernel;
            ++index;
        }
    }
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, 
    std::vector<double>& rOutput, 
    const ProcessInfo& rProcessInfo
)
{
    if (rVariable == SPH_KERNEL){
        const auto& r_neighbours = this->GetValue(NEIGHBOURS);
        const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

        const int domain_size = GetGeometry().WorkingSpaceDimension();
        const auto& IPcoords = GetGeometry()[0].GetInitialPosition();

        IndexType index = 0;
        rOutput.resize(r_neighbours.size());

        for (auto& JP : r_neighbours){
            const auto& JPcoords = JP->GetGeometry()[0].GetInitialPosition();

            VectorType X_AB_target(domain_size);
            for (IndexType d = 0; d < domain_size; d++){
                X_AB_target[d] = IPcoords[d] - JPcoords[d];
            }

            double kernel;
            TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
            rOutput[index] = kernel;
            ++index;
        }
    }
}


template<class TKernelType>
int TotalLagrangianDisplacementParticle<TKernelType>::Check(const ProcessInfo& rCurrentProcessInfo) const 
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1)<<"Particle element with invalid Id"<< this->Id() << std::endl;
    KRATOS_ERROR_IF(this->GetGeometry().size() != 1) << "Particle element" << this->Id() << "must have exactly 1 node, found " << this->GetGeometry().size() << std::endl;

    return 0; 
    KRATOS_CATCH("")
}

template class TotalLagrangianDisplacementParticle<CubicKernel2D>;

}

