#include "custom_elements/total_lagrangian_mixed_strain_particle.h"

#include "constitutive_laws_application_variables.h"

namespace Kratos
{

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rCurrentProcessInfo[IS_RESTARTED]){
        InitializeMaterial();
    }

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::InitializeMaterial()
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
void TotalLagrangianMixedStrainParticle<TKernelType>::InitializeSolutionStep(const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    bool required = false;
    if (mThisConstitutiveLaw->RequiresInitializeMaterialResponse()) required = true;

    if (required){
        const auto& r_geom = GetGeometry();
        const auto& r_prop = GetProperties();
        const SizeType dimension = r_geom.WorkingSpaceDimension();
        const SizeType strain_size = mThisConstitutiveLaw->GetStrainSize();
        const SizeType number_of_neighbours = GetValue(NEIGHBOURS).size();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_neighbours);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        ConstitutiveLaw::Parameters cl_values(r_geom, r_prop, rProcessInfo);
        // Set constitutive law flags
        Flags& cl_options = cl_values.GetOptions(); 
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        CalculateKinematicVariables(this_kinematic_variables, rProcessInfo);

        SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values);
        cl_values.SetStrainVector(this_constitutive_variables.StrainVector);

        // Call the constitutive law to update material variables 
        mThisConstitutiveLaw->InitializeMaterialResponse(cl_values, ConstitutiveLaw::StressMeasure_PK2);
    }

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::FinalizeSolutionStep(const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    bool required = false;
    if (mThisConstitutiveLaw->RequiresInitializeMaterialResponse()) required = true;

    if (required){
        const auto& r_geom = GetGeometry();
        const auto& r_prop = GetProperties();
        const SizeType dimension = r_geom.WorkingSpaceDimension();
        const SizeType strain_size = mThisConstitutiveLaw->GetStrainSize();
        const SizeType number_of_neighbours = GetValue(NEIGHBOURS).size();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_neighbours);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        ConstitutiveLaw::Parameters cl_values(r_geom, r_prop, rProcessInfo);
        // Set constitutive law flags
        Flags& cl_options = cl_values.GetOptions(); 
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        CalculateKinematicVariables(this_kinematic_variables, rProcessInfo);

        SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values);
        cl_values.SetStrainVector(this_constitutive_variables.StrainVector);

        // Call the constitutive law to update material variables 
        mThisConstitutiveLaw->FinalizeMaterialResponse(cl_values, ConstitutiveLaw::StressMeasure_PK2);
    }

    // Forse qui devi aggiornare posizione

    KRATOS_CATCH("")
}


template<class TKernelType>
Element::Pointer TotalLagrangianMixedStrainParticle<TKernelType>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
) const
{
    KRATOS_TRY
    
    TotalLagrangianMixedStrainParticle<TKernelType>::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangianMixedStrainParticle<TKernelType>>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    p_new_elem->SetConstitutiveLaw(mThisConstitutiveLaw);

    return p_new_elem;

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    const SizeType dofs_per_node = dimension + (dimension * dimension);

    if (rResult.size() != dofs_per_node * number_of_neighbours)
        rResult.resize(dofs_per_node * number_of_neighbours, false);

    if (dimension == 2){
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * dimension;
            const SizeType F_block = i * dimension * dimension + dimension * number_of_neighbours;

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
            const SizeType v_block = i * dimension;
            const SizeType F_block = i * dimension * dimension + dimension * number_of_neighbours;

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

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    const SizeType dofs_per_node = dimension + (dimension * dimension);

    if(rElementalDofList.size() != dofs_per_node * number_of_neighbours)
        rElementalDofList.resize(dofs_per_node * number_of_neighbours);
    
    if (dimension == 2){
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * dimension;
            const SizeType F_block = i * dimension * dimension + dimension * number_of_neighbours;

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
            const SizeType v_block = i * dimension;
            const SizeType F_block = i * dimension * dimension + dimension * number_of_neighbours;

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

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_TRY
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_per_node = dimension + (dimension * dimension);
    const SizeType mat_size = dofs_per_node * number_of_neighbours;

    if (rValues.size() != mat_size) 
        rValues.resize(mat_size, false);

    if (dimension == 2){
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * dimension;
            const SizeType F_block = i * dimension * dimension + dimension * number_of_neighbours;

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
            const SizeType v_block = i * dimension;
            const SizeType F_block = i * dimension * dimension + dimension * number_of_neighbours;

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

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::GetFirstDerivativesVector(VectorType& rValues, int step) const
{
    KRATOS_TRY
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_per_node = dimension + (dimension * dimension);
    const SizeType mat_size = dofs_per_node * number_of_neighbours;

    if (rValues.size() != mat_size) rValues.resize(mat_size, false);

    if (dimension == 2){
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            const auto& r_geom = r_neighbours[i]->GetGeometry();
            const SizeType v_block = i * dimension;
            const SizeType F_block = i * dimension * dimension + dimension * number_of_neighbours;

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
            const SizeType v_block = i * dimension;
            const SizeType F_block = i * dimension * dimension + dimension * number_of_neighbours;

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

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateLocalSystem(
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
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateLeftHandSide(
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
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateRightHandSide(
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
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateAll(
    MatrixType& rLHS, 
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const auto& r_geom = GetGeometry(); 
    const auto& r_prop = GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType strain_size = mThisConstitutiveLaw->GetStrainSize();
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType block_size = dimension + dimension * dimension; 
    const SizeType mat_size = block_size * number_of_neighbours;

    if (CalculateStiffnessMatrixFlag){
        if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) rLHS.resize(mat_size, mat_size, false);
        noalias(rLHS) = ZeroMatrix(mat_size, mat_size);
    }

    if (CalculateResidualVectorFlag){
        if (rRHS.size() != mat_size) rRHS.resize(mat_size, false);
        noalias(rRHS) = ZeroVector(mat_size);
    }

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_neighbours);
    ConstitutiveVariables this_constitutive_variables(strain_size);

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
    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, ConstitutiveLaw::StressMeasure_PK2);

    if (CalculateStiffnessMatrixFlag){
        // Initializing the blocks
        MatrixType K11(dimension * number_of_neighbours, dimension * number_of_neighbours); K11.clear(); 
        MatrixType K12(dimension * number_of_neighbours, dimension * dimension * number_of_neighbours); K12.clear();
        MatrixType K21(dimension * dimension * number_of_neighbours, dimension * number_of_neighbours); K21.clear();
        MatrixType K22(dimension * dimension * number_of_neighbours, dimension * dimension * number_of_neighbours); K22.clear();

        /* Geometric stiffness matrix */
        //CalculateAndAddKg(K12, this_kinematic_variables.DW_DX, this_constitutive_variables.StressVector, gauss_weight);

        /* Material stiffness matrix */
        //CalculateAndAddKm(K12, this_kinematic_variables.B, this_constitutive_variables.C, gauss_weight); 

        AssembleLHS(rLHS, K11, K12, K21, K22);
    }

    if (CalculateResidualVectorFlag){
        VectorType RHSv(dimension * number_of_neighbours); RHSv.clear();
        VectorType RHSF(dimension * dimension * number_of_neighbours); RHSF.clear();

        CalculateLinearMomentumResidualVector(RHSv, this_kinematic_variables, rProcessInfo, this_constitutive_variables.StressVector, gauss_weight);

        CalculateGeometricalResidualVector(RHSF, this_kinematic_variables, gauss_weight);

        AssembleRHS(rRHS, RHSv, RHSF);
    }

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateAndAddKg(
    MatrixType& rK, 
    const Matrix& DW_DX, 
    const Vector& stress_vector, 
    const double weight
    ) const 
{
    KRATOS_TRY
    const SizeType domain_size = GetGeometry().WorkingSpaceDimension();
    const MatrixType stress_tensor = weight * MathUtils<double>::StressVectorToTensor(stress_vector);
    MatrixType reduced_Kg(DW_DX.size1(), DW_DX.size1());
    MathUtils<double>::BDBtProductOperation(reduced_Kg, stress_tensor, DW_DX);
    MathUtils<double>::ExpandAndAddReducedMatrix(rK, reduced_Kg, domain_size);
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateAndAddKm(
    MatrixType& rK, 
    const Matrix& rB,
    const Matrix& rConstitutiveMatrix, 
    const double weight
    ) const 
{
    KRATOS_TRY
    noalias(rK) += weight * prod(trans(rB), Matrix(prod(rConstitutiveMatrix, rB)));
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateLinearMomentumResidualVector(
    VectorType& rRHSv,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rProcessInfo,
    const VectorType& rStressVector,
    const double weight
)
{
    KRATOS_TRY

    const SizeType number_of_neigh = GetValue(NEIGHBOURS).size();
    const SizeType domain_size = GetGeometry().WorkingSpaceDimension();
    
    // Only external forces contributions are taken into account at the moment 
    VectorType body_force(domain_size); SPHElementUtilities::GetLocalBodyForces(*this, body_force);

    for (IndexType i = 0; i < number_of_neigh; ++i){
        const SizeType index = i * domain_size;
        for (IndexType d = 0; d < domain_size; ++d)
            rRHSv[index + d] += weight * rThisKinematicVariables.W[i] * body_force[d];
    }

    // Adding stress contribution to the residual vector
    noalias(rRHSv) -= weight * prod(trans(rThisKinematicVariables.B), rStressVector);

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateGeometricalResidualVector(
    VectorType& rRHSF,
    KinematicVariables& rThisKinematicVariables,
    const double weight
)
{
    KRATOS_TRY

    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neigh = r_neighbours.size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    VectorType temp_residual(dimension * dimension), vel_aux(dimension);
    MatrixType temp = ZeroMatrix(dimension, dimension); 

    for (IndexType i = 0; i < number_of_neigh; ++i){
        
        const array_1d<double, 3>& velocity = r_neighbours[i]->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        for (IndexType d = 0; d < dimension; ++d) vel_aux[d] = velocity[d];

        temp += outer_prod(vel_aux, row(rThisKinematicVariables.DW_DX, i));
    }

    temp_residual = SPHElementUtilities::NonSymmetricTensorToVector(temp); 

    for (IndexType i = 0; i < number_of_neigh; ++i){
        noalias(project(rRHSF, range(dimension * dimension * i, dimension * dimension * (i + 1)))) += temp_residual * rThisKinematicVariables.W[i];
    }

    rRHSF *= weight;

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateLeftHandSideK21Block(
    MatrixType& rK21,
    const KinematicVariables& rThisKinematicVariables,
    const double weight
)
{
    KRATOS_TRY

    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neigh = r_neighbours.size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    const MatrixType E = IdentityMatrix(dimension);
    VectorType e(dimension), temp_vec(dimension * dimension);
    MatrixType temp(dimension, dimension);

    for (IndexType i = 0; i < number_of_neigh; ++i){
        for (IndexType d = 0; d < dimension; ++d){
            
            e = column(E, d); 
            temp = outer_prod(e, row(rThisKinematicVariables.DW_DX, i));
            temp_vec = SPHElementUtilities::NonSymmetricTensorToVector(temp);
            
            const SizeType col_index = i * dimension + d;
            
            for (IndexType ii = 0; ii < number_of_neigh; ++ii){
                
                const SizeType row_start = dimension * dimension * ii;
                
                for (IndexType k = 0; k < dimension * dimension; ++k){
                    rK21(row_start + k, col_index) += weight * rThisKinematicVariables.W[ii] * temp_vec[k];
                }
            }
        }
    }

    KRATOS_CATCH("")   
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::AssembleLHS(
    MatrixType& rLHS,
    const MatrixType& rK11,
    const MatrixType& rK21,
    const MatrixType& rK12,
    const MatrixType& rK22
    )
{
    KRATOS_TRY

    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_neighbours = GetValue(NEIGHBOURS).size();

    const SizeType first_equation_dofs = number_of_neighbours * dimension;

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

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::AssembleRHS(
    VectorType& rRHS,
    const VectorType& rRHSv,
    const VectorType& rRHSF
    )
{
    KRATOS_TRY

    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_neighbours = GetValue(NEIGHBOURS).size();

    const SizeType first_equation_dofs = number_of_neighbours * dimension;
    
    for (IndexType i = 0; i < rRHSv.size(); ++i)
        rRHS[i] = rRHSv[i];
    
    for (IndexType i = 0; i < rRHSF.size(); ++i)
        rRHS[first_equation_dofs + i] = rRHSF[i];

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY
    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType block_size = dimension + dimension * dimension;
    const SizeType mat_size = block_size * number_of_neighbours;
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize(mat_size, mat_size, false);
    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size); 

    KRATOS_ERROR_IF(r_prop.Has(DENSITY) == false) << "DENSITY not provided for element " << this->Id() << std::endl;

    const double density = r_prop[DENSITY];
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    const auto& GPcoords = r_geom[0].GetInitialPosition();
    const double gauss_weight = r_geom[0].GetValue(VOLUME);

    double factor = density * thickness * gauss_weight;
    bool compute_lumped_mass_matrix = SPHElementUtilities::ComputeLumpedMassMatrix(r_prop, rProcessInfo);

    if (compute_lumped_mass_matrix){
        int this_id = this->Id();
        for (IndexType i = 0; i < number_of_neighbours; ++i){
            if (r_neighbours[i]->Id() == this_id){
                const SizeType v_block = i * dimension;
                const SizeType F_block = i * dimension * dimension + dimension * number_of_neighbours;
                
                for (IndexType d = 0; d < dimension; ++d)
                    rMassMatrix(v_block + d, v_block + d) = factor;

                for (IndexType d = 0; d < dimension * dimension; ++d)
                    rMassMatrix(F_block + d, F_block + d) = gauss_weight;
            }
        }
    } else { // Consistent mass matrix
        
        MatrixType vBlock(mat_size, block_size); vBlock.clear();
        MatrixType FBlock(mat_size, block_size); FBlock.clear();
        
        VectorType X_AB_target(dimension);
        double kernel;

        for (IndexType i = 0; i < number_of_neighbours; ++i){

            const auto& JPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
            const double volume = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);

            for (IndexType d = 0; d < dimension; d++){
                X_AB_target[d] = GPcoords[d] - JPcoords[d];
            }

            TKernelType::ComputeKernelValue(kernel, h, X_AB_target); 
            ComputeKernelCorrectionUtilities::ApplyKernelCorrection(*this, kernel);

            for (IndexType d = 0; d < dimension; ++d){
                vBlock(dimension * i + d, d) = kernel * volume;
            }

            for (IndexType d = 0; d < dimension * dimension; ++d){
                FBlock(dimension * dimension * i + d, d) = kernel * volume;
            }
        }

        // Aseembling the mass matrix
        int v_size = dimension * number_of_neighbours;
        int F_size = dimension * dimension * number_of_neighbours;

        noalias(project(rMassMatrix, range(0, v_size), range(0, v_size))) += factor * prod(vBlock, trans(vBlock));
        noalias(project(rMassMatrix, range(v_size, v_size + F_size), range(v_size, v_size + F_size))) += gauss_weight * prod(FBlock, trans(FBlock));
    }
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues
    ) const
{
    // Essential input parameters for the constitutive law
    rValues.SetDeterminantF(rThisKinematicVariables.detF);
    rValues.SetDeformationGradientF(rThisKinematicVariables.F);

    // Space in which the resulta shall be saved
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.C);
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector);
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues);

    // Actually do the computations in the ConstitutiveLaw
    mThisConstitutiveLaw->CalculateMaterialResponse(rValues, ThisStressMeasure);
}
    
template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY

    CalculateKernelsAndKernelGradients(rThisKinematicVariables.DW_DX, rThisKinematicVariables.W, rProcessInfo);
    AssembleDeformationGradient(rThisKinematicVariables.F);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);

    const SizeType number_of_neighbours = GetValue(NEIGHBOURS).size();

    if (GetGeometry().WorkingSpaceDimension() == 2){
        SPHElementUtilities::Calculate2DB(rThisKinematicVariables.B, rThisKinematicVariables.F, rThisKinematicVariables.DW_DX, number_of_neighbours);
    } else {
        SPHElementUtilities::Calculate3DB(rThisKinematicVariables.B, rThisKinematicVariables.F, rThisKinematicVariables.DW_DX, number_of_neighbours);
    }

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateKernelsAndKernelGradients(
    MatrixType& rDW_DX,
    VectorType& rW,
    const ProcessInfo& rProcessInfo
    )
{
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);
    const int domain_size = GetGeometry().WorkingSpaceDimension();
    const auto& IPcoords = GetGeometry()[0].GetInitialPosition();

    // Initialization of variables
    double kernel;
    VectorType dkernel(domain_size), X_AB_target(domain_size);

    for (IndexType i = 0; i < r_neighbours.size(); ++i){
        
        const auto& JPcoords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
        const auto& JP_current_position = r_neighbours[i]->GetGeometry()[0].Coordinates();
        
        const double weight = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);

        for (IndexType d = 0; d < domain_size; d++) 
            X_AB_target[d] = IPcoords[d] - JPcoords[d];

        TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
        TKernelType::ComputeKernelGradientValue(dkernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, kernel, dkernel);
            
        for (IndexType d = 0; d < domain_size; d++) 
            rDW_DX(i, d) = weight * dkernel[d];

        rW[i] = weight * kernel;
    }
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::AssembleDeformationGradient(
    MatrixType& rF
)
{
    const auto& r_geom = GetGeometry();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    
    if (dimension == 2){
        rF(0,0) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XX);
        rF(1,1) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YY);
        rF(0,1) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XY);
        rF(1,0) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YX);
    } else {
        rF(0,0) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XX);
        rF(1,1) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YY);
        rF(2,2) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_ZZ); 
        rF(0,1) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XY);
        rF(0,2) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_XZ);
        rF(1,0) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YX);
        rF(1,2) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_YZ);
        rF(2,0) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_ZX);
        rF(2,1) = r_geom[0].FastGetSolutionStepValue(DEFORMATION_GRADIENT_ZY);
    }
}

template<class TKernelType>
int TotalLagrangianMixedStrainParticle<TKernelType>::Check(const ProcessInfo& rProcessInfo) const 
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1)<<"Particle element with invalid Id"<< this->Id() << std::endl;
    KRATOS_ERROR_IF(this->GetGeometry().size() != 1) << "Particle element" << this->Id() << "must have exactly 1 node, found " << this->GetGeometry().size() << std::endl;

    const auto& r_prop = GetProperties();

    if (r_prop[CONSTITUTIVE_LAW] != nullptr) mThisConstitutiveLaw->Check( r_prop, GetGeometry(), rProcessInfo );

    return 0; 
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateOnIntegrationPoints(
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
    } else if (mThisConstitutiveLaw->Has(rVariable)){  // At the moment this funtion is never called because the Point2D/Point3D geometries does not have an integration point
        const SizeType strain_size = mThisConstitutiveLaw->GetStrainSize();
        rOutput.resize(1, ZeroVector(strain_size));
        mThisConstitutiveLaw->GetValue(rVariable, rOutput[0]);
    }
}

template<class TKernelType>
void TotalLagrangianMixedStrainParticle<TKernelType>::CalculateOnIntegrationPoints(
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
    
template class TotalLagrangianMixedStrainParticle<CubicKernel2D>; 
template class TotalLagrangianMixedStrainParticle<CubicKernel3D>;

}