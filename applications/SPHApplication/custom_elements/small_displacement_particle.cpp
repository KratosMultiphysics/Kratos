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
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
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
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
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
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
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
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
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
/*
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
    const auto&  r_geom = this->GetGeometry();
    const auto& r_props = this->GetProperties();
    const SizeType domain_size = r_geom.WorkingSpaceDimension();
    const SizeType strain_size = this->mThisConstitutiveLaw->GetStrainSize();
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    // The neighbours of the particle are like the nodes belonging to the element in standard FEM
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType mat_size = domain_size * r_neighbours.size();
    
    if (CalculateStiffnessMatrixFlag){
        if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) rLHS.resize(mat_size, mat_size, false);
        noalias(rLHS) = ZeroMatrix(mat_size, mat_size);
    }

    if (CalculateResidualVectorFlag){
        if (rRHS.size() != mat_size) rRHS.resize(mat_size, false);
        noalias(rRHS) = ZeroVector(mat_size);
    }

    VectorType nodal_values(mat_size);
    this->GetNodalValuesVector(nodal_values);

    // Kernel variables initialization 
    double ip_kernel, jp_kernel;
    VectorType ip_dkernel(domain_size), jp_dkernel(domain_size), X_AB_target(domain_size);

    // Constitutive law 
    ConstitutiveLaw::Parameters cl_values(r_geom, r_props, rProcessInfo);
    auto& r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Constitutive law initialization
    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    // Initilization of other things
    MatrixType K_loc_temp(domain_size, domain_size), B_b(strain_size, domain_size), 
        B_a(strain_size, domain_size), temp(strain_size, domain_size);
    VectorType f_local(domain_size), local_nodal_values(domain_size);

    const auto& GPcoords = r_geom[0].GetInitialPosition();
    const double volume_gp = r_geom[0].GetValue(VOLUME);

    IndexType ip_index = 0;
    
    for (auto& IP : r_neighbours){

        const auto& geom_ip = IP->GetGeometry();
        const auto& IPcoords = geom_ip[0].GetInitialPosition();
        double volume_ip = geom_ip[0].GetValue(VOLUME);
            
        for (IndexType d = 0; d < domain_size; d++){
            X_AB_target[d] = GPcoords[d] - IPcoords[d];
        }

        TKernelType::ComputeKernelValue(ip_kernel, h, X_AB_target);
        TKernelType::ComputeKernelGradientValue(ip_dkernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, ip_kernel, ip_dkernel);

        if (CalculateResidualVectorFlag){
            array_1d<double, 2> body_force = this->GetLocalBodyForces();
            f_local = volume_gp * volume_ip * body_force * ip_kernel; 
            noalias(project(rRHS, range(domain_size * ip_index, domain_size * (ip_index + 1)))) += f_local;
        }

        this->GlobalSizeVector(local_nodal_values, nodal_values, ip_index);
        strain_vector.clear();
        this->CalculateStrainVector(strain_vector, ip_dkernel, local_nodal_values, volume_ip);

        this->mThisConstitutiveLaw->CalculateMaterialResponseCauchy(cl_values);
        const VectorType& r_stress_vector = cl_values.GetStressVector();
        const MatrixType& r_constitutive_matrix = cl_values.GetConstitutiveMatrix();

        // Local system computation and assembly
        this->GetShapeFunctionDerivatives(B_a, ip_dkernel, volume_ip); 

        IndexType jp_index = 0;

        for (auto& JP : r_neighbours){

            const auto& geom_jp = JP->GetGeometry();
            const auto& JPcoords = geom_jp[0].GetInitialPosition();
            double volume_jp = geom_jp[0].GetValue(VOLUME);
                
            for (IndexType d = 0; d < domain_size; d++){
                X_AB_target[d] = GPcoords[d] - JPcoords[d];
            }

            TKernelType::ComputeKernelValue(jp_kernel, h, X_AB_target);
            TKernelType::ComputeKernelGradientValue(jp_dkernel, h, X_AB_target);
            ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, jp_kernel, jp_dkernel);

            // Local system computation and assembly
            this->GetShapeFunctionDerivatives(B_b, jp_dkernel, volume_jp);

            f_local = prod(trans(B_b), r_stress_vector);
            f_local *= volume_gp;

            temp = prod(r_constitutive_matrix, B_a);
            K_loc_temp = prod(trans(B_b), temp);
            K_loc_temp *= volume_gp;
                
            if (CalculateStiffnessMatrixFlag){
                noalias(project(rLHS, range(domain_size * jp_index, domain_size * (jp_index + 1)), range(domain_size * ip_index, domain_size * (ip_index + 1)))) += K_loc_temp;
            }
            
            if (CalculateResidualVectorFlag){
                noalias(project(rRHS, range(domain_size * jp_index, domain_size * (jp_index + 1)))) -= f_local;
            }

            jp_index++;
        }
        ip_index++;
    }
    KRATOS_CATCH("")
} */


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
    const SizeType strain_size = this->mThisConstitutiveLaw->GetStrainSize();
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    // The neighbours of the particle are like the nodes belonging to the element in standard FEM
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
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

    // Precomputation of the corrected kernels and kernels gradients  

    const auto& GPcoords = r_geom[0].GetInitialPosition();
    const double volume_gp = r_geom[0].GetValue(VOLUME);

    // Kernel variables initialization 
    double ip_kernel;
    VectorType ip_dkernel(domain_size), X_AB_target(domain_size);
    std::vector<VectorType> dkernel_storage(number_of_neigh, VectorType(domain_size));
    VectorType kernel_storage(number_of_neigh), volume_storage(number_of_neigh);
    

    for (IndexType index = 0; index < number_of_neigh; ++index){
        const auto& geom_ip = r_neighbours[index]->GetGeometry();
        const auto& IPcoords = geom_ip[0].GetInitialPosition();
        double volume = geom_ip[0].GetValue(VOLUME);
            
        for (IndexType d = 0; d < domain_size; d++){
            X_AB_target[d] = GPcoords[d] - IPcoords[d];
        }

        TKernelType::ComputeKernelValue(ip_kernel, h, X_AB_target);
        TKernelType::ComputeKernelGradientValue(ip_dkernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, ip_kernel, ip_dkernel);

        kernel_storage[index] = ip_kernel;
        dkernel_storage[index] = ip_dkernel;
        volume_storage[index] = volume;
    }

    VectorType nodal_values(mat_size);
    GetNodalValuesVector(nodal_values);

    // Constitutive law 
    ConstitutiveLaw::Parameters cl_values(r_geom, r_props, rProcessInfo);
    auto& r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Constitutive law initialization
    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    // DA VEDERE SE SERVE
    MatrixType K_loc_temp(domain_size, domain_size), B_b(strain_size, domain_size), 
        B_a(strain_size, domain_size), temp(strain_size, domain_size);
    VectorType f_local(domain_size), local_nodal_values(domain_size);
    //////////
    VectorType body_force(domain_size);

    for (IndexType ip = 0; ip < number_of_neigh; ++ip){
        
        if (CalculateResidualVectorFlag){
            SPHElementUtilities::GetLocalBodyForces(*this, body_force);
            f_local = volume_gp * volume_storage[ip] * body_force * kernel_storage[ip]; 
            noalias(project(rRHS, range(domain_size * ip, domain_size * (ip + 1)))) += f_local;
        }

        GlobalSizeVector(local_nodal_values, nodal_values, ip);
        strain_vector.clear();
        CalculateStrainVector(strain_vector, dkernel_storage[ip], local_nodal_values, volume_storage[ip]);

        mThisConstitutiveLaw->CalculateMaterialResponseCauchy(cl_values);
        const VectorType& r_stress_vector = cl_values.GetStressVector();
        const MatrixType& r_constitutive_matrix = cl_values.GetConstitutiveMatrix();

        // Local system computation and assembly
        GetShapeFunctionDerivatives(B_a, dkernel_storage[ip], volume_storage[ip]); 
        
        for (IndexType jp = 0; jp < number_of_neigh; ++jp){

            // Local system computation and assembly
            GetShapeFunctionDerivatives(B_b, dkernel_storage[jp], volume_storage[jp]);

            f_local = prod(trans(B_b), r_stress_vector);
            f_local *= volume_gp;

            temp = prod(r_constitutive_matrix, B_a);
            K_loc_temp = prod(trans(B_b), temp);
            K_loc_temp *= volume_gp;
            
            if (CalculateStiffnessMatrixFlag){
                noalias(project(rLHS, range(domain_size * jp, domain_size * (jp + 1)), range(domain_size * ip, domain_size * (ip + 1)))) += K_loc_temp;
            }
            
            if (CalculateResidualVectorFlag){
                noalias(project(rRHS, range(domain_size * jp, domain_size * (jp + 1)))) -= f_local;
            }

        }
    }

    KRATOS_CATCH("")
}

template<class TKernelType>
void SmallDisplacementParticle<TKernelType>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
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

    MatrixType J0(dimension, dimension);
    noalias(J0) = ZeroMatrix(dimension, dimension);
    VectorType kernel_storage(number_of_neigh), X_AB_target(dimension);
    double kernel;

    for (IndexType i = 0; i < number_of_neigh; ++i){
        
        const auto& IPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
        const double volume = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);
        
        for (IndexType d = 0; d < dimension; ++d){
            X_AB_target[d] = GPcoords[d] - IPcoords[d];
        }

        TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelCorrection(*this, kernel);

        kernel_storage[i] = volume * kernel;
    }

    double temp; 
    double factor = gauss_weight * density * thickness;

    for (IndexType i = 0; i < number_of_neigh; ++i){
        for (IndexType j = 0; j < number_of_neigh; ++j){   
            temp = factor * kernel_storage[i] * kernel_storage[j];
            noalias(J0) = IdentityMatrix(dimension);
            J0 *= temp;
            noalias(project(rMassMatrix, range(dimension * j, dimension * (j + 1)), range(dimension * i, dimension * (i + 1)))) += J0;
        }
    }

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





