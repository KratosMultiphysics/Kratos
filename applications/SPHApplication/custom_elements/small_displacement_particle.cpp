// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

// Application includes
#include "custom_elements/small_displacement_particle.h"

namespace Kratos
{

template<class TKernelType>
void SmallDisplacementParticle<TKernelType>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    InitializeMaterial();

    KRATOS_CATCH("")
}

template<class TKernelType>
void SmallDisplacementParticle<TKernelType>::InitializeMaterial()
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
Element::Pointer SmallDisplacementParticle<TKernelType>::Clone( 
    IndexType NewId, 
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    SmallDisplacementParticle<TKernelType>::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementParticle<TKernelType>>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    p_new_elem->SetConstitutiveLaw(mThisConstitutiveLaw);

    return p_new_elem;

    KRATOS_CATCH("");
}

template<class TKernelType>
void SmallDisplacementParticle<TKernelType>::GetNodalValuesVector(VectorType& rNodalValues) const
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
void SmallDisplacementParticle<TKernelType>::EquationIdVector(
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
void SmallDisplacementParticle<TKernelType>::GetDofList(
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
void SmallDisplacementParticle<TKernelType>::GetFirstDerivativesVector(VectorType& rValues, int step) const
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
void SmallDisplacementParticle<TKernelType>::GetSecondDerivativesVector(VectorType& rValues, int step) const
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
void SmallDisplacementParticle<TKernelType>::GetShapeFunctionDerivatives(
    MatrixType& rB,
    VectorType& rCGCK,
    const double volume
)
{
    rB.clear();
    rB(0,0) = rCGCK[0];
    rB(1,1) = rCGCK[1];
    rB(2,0) = rCGCK[1];
    rB(2,1) = rCGCK[0];
    rB *= volume;

}

template<class TKernelType>
void SmallDisplacementParticle<TKernelType>::CalculateStrainVector(
    VectorType& rStrainVector, 
    VectorType& rCGCK, 
    const VectorType& rNodalValues,
    const double volume
)
{
    SizeType strain_size = rStrainVector.size();
    SizeType domain_size = rCGCK.size();
    MatrixType B(strain_size, domain_size);

    rStrainVector.clear();
    GetShapeFunctionDerivatives(B, rCGCK, volume);
    noalias(rStrainVector) = prod(B, rNodalValues);
}

template<class TKernelType>
void SmallDisplacementParticle<TKernelType>::CalculateLocalSystem(
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
void SmallDisplacementParticle<TKernelType>::CalculateLeftHandSide(
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
void SmallDisplacementParticle<TKernelType>::CalculateRightHandSide(
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
void SmallDisplacementParticle<TKernelType>::CalculateAll(
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
    
    if (CalculateStiffnessMatrixFlag){
        if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) rLHS.resize(mat_size, mat_size, false);
        noalias(rLHS) = ZeroMatrix(mat_size, mat_size);
    }

    if (CalculateResidualVectorFlag){
        if (rRHS.size() != mat_size) rRHS.resize(mat_size, false);
        noalias(rRHS) = ZeroVector(mat_size);
    }

    VectorType nodal_values(mat_size), dkernel(domain_size), X_AB_target(domain_size);
    GetNodalValuesVector(nodal_values);
    double kernel;

    // Constitutive law 
    ConstitutiveLaw::Parameters cl_values(r_geom, r_props, rProcessInfo);
    auto& r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true); //To compute the strain outside the CL 

    // Constitutive law initialization
    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);
    cl_values.SetDeformationGradientF(IdentityMatrix(2)); // To maintain compatibility with Kratos CLs
    cl_values.SetDeterminantF(1.0); // To maintain compatibility with Kratos CLs

    // Initilization of other things
    MatrixType B(strain_size, domain_size), temp(strain_size, domain_size);
    VectorType f_local(domain_size), local_nodal_values(domain_size), body_force(domain_size);

    MatrixType B_global(domain_size * number_of_neigh, strain_size), CB_global(strain_size, domain_size * number_of_neigh),
        strain_storage(strain_size, number_of_neigh);
    noalias(B_global) = ZeroMatrix(domain_size * number_of_neigh, strain_size);
    noalias(CB_global) = ZeroMatrix(strain_size, domain_size * number_of_neigh);
    noalias(strain_storage) = ZeroMatrix(strain_size, number_of_neigh);

    const auto& GPcoords = r_geom[0].GetInitialPosition();
    const double gauss_weight = r_geom[0].GetValue(VOLUME);

    for (IndexType i = 0; i < number_of_neigh; ++i){
        
        const auto& IPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
        double volume = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);
            
        for (IndexType d = 0; d < domain_size; d++){
            X_AB_target[d] = GPcoords[d] - IPcoords[d];
        }

        TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
        TKernelType::ComputeKernelGradientValue(dkernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, kernel, dkernel);

        if (CalculateResidualVectorFlag){
            SPHElementUtilities::GetLocalBodyForces(*this, body_force);
            noalias(project(rRHS, range(domain_size * i, domain_size * (i + 1)))) += gauss_weight * volume * body_force * kernel;
        }

        GlobalSizeVector(local_nodal_values, nodal_values, i);
        strain_vector.clear();
        CalculateStrainVector(strain_vector, dkernel, local_nodal_values, volume);

        mThisConstitutiveLaw->CalculateMaterialResponseCauchy(cl_values);
        const VectorType& r_stress_vector = cl_values.GetStressVector();
        const MatrixType& r_constitutive_matrix = cl_values.GetConstitutiveMatrix();
        
        GetShapeFunctionDerivatives(B, dkernel, volume);

        noalias(project(B_global, range(i * domain_size, (i + 1) * domain_size), range(0, strain_size))) += trans(B);
        
        temp = prod(r_constitutive_matrix, B);
        noalias(project(CB_global, range(0, strain_size), range(i * domain_size, (i + 1) * domain_size))) += temp;

        column(strain_storage, i) = r_stress_vector;
    }

    if (CalculateStiffnessMatrixFlag){
        noalias(rLHS) = prod(B_global, CB_global);
        rLHS *= gauss_weight;
    }

    if (CalculateResidualVectorFlag){
        VectorType ones_vector(number_of_neigh, 1.0);
        VectorType internal_forces = prod(prod(B_global, strain_storage), ones_vector);
        internal_forces *= gauss_weight;
        noalias(rRHS) -= internal_forces;
    }

    KRATOS_CATCH("")
}

template<class TKernelType>
void SmallDisplacementParticle<TKernelType>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rProcessInfo)
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
void SmallDisplacementParticle<TKernelType>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rProcessInfo)
{
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType mat_size = GetGeometry().WorkingSpaceDimension() * r_neighbours.size();

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(*this, rDampingMatrix, rProcessInfo, mat_size);
}

template<class TKernelType>
void SmallDisplacementParticle<TKernelType>::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rProcessInfo)
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
void SmallDisplacementParticle<TKernelType>::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rProcessInfo)
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
int SmallDisplacementParticle<TKernelType>::Check(const ProcessInfo& rCurrentProcessInfo) const 
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1)<<"Particle element with invalid Id"<< this->Id() << std::endl;
    KRATOS_ERROR_IF(this->GetGeometry().size() != 1) << "Particle element" << this->Id() << "must have exactly 1 node, found " << this->GetGeometry().size() << std::endl;

    return 0; 
    KRATOS_CATCH("")
}

template class SmallDisplacementParticle<CubicKernel2D>;

}





