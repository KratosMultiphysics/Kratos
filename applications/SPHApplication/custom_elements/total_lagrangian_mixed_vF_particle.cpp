#include "custom_elements/total_lagrangian_mixed_vF_particle.h"

#include "constitutive_laws_application_variables.h"

namespace Kratos
{

template<class TKernelType>
void TotalLagrangianMixedvFParticle<TKernelType>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    InitializeMaterial();

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedvFParticle<TKernelType>::InitializeMaterial()
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
void TotalLagrangianMixedvFParticle<TKernelType>::InitializeSolutionStep(const ProcessInfo& rProcessInfo)
{
    bool required = false;
    if (mThisConstitutiveLaw->RequiresInitializeMaterialResponse()) required = true;

    if (required){
        KRATOS_ERROR << "InitializeSolutionStep not implemented"
    }
}

template<class TKernelType>
void TotalLagrangianDisplacementParticle<TKernelType>::FinalizeSolutionStep(const ProcessInfo& rProcessInfo)
{
    bool required = false;
    if (mThisConstitutiveLaw->RequiresFinalizeMaterialResponse()) required = true;

    if (required) {
        KRATOS_ERROR << "FinalizeSolutionStep not implemented"
    }
}


template<class TKernelType>
Element::Pointer TotalLagrangianMixedvFParticle<TKernelType>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
) const
{
    KRATOS_TRY
    
    TotalLagrangianMixedvFParticle<TKernelType>::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangianMixedvFParticle<TKernelType>>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    p_new_elem->SetConstitutiveLaw(mThisConstitutiveLaw);

    return p_new_elem;

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedvFParticle<TKernelType>::EquationIdVector(
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

    IndexType local_index = 0;

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        const auto& r_geom = r_neighbours[i]->GetGeometry();

        const IndexType xpos = r_geom[0].GetDofPosition(DISPLACEMENT_X);
        const IndexType Fpos = r_geom[0].GetDofPosition(DEFORMATION_GRADIENT_XX);
        
        rResult[local_index++] = r_geom[0].GetDof(DISPLACEMENT_X, xpos).EquationId();
        rResult[local_index++] = r_geom[0].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
        if (dimension == 3)
            rResult[local_index++] = r_geom[0].GetDof(DISPLACEMENT_Z, xpos + 2).EquationId();

        rResult[local_index++] = r_geom[0].GetDof(DEFORMATION_GRADIENT_XX, Fpos).EquationId();
        rResult[local_index++] = r_geom[0].GetDof(DEFORMATION_GRADIENT_YY, Fpos + 1).EquationId();
        rResult[local_index++] = r_geom[0].GetDof(DEFORMATION_GRADIENT_XY, Fpos + 2).EquationId();
        rResult[local_index++] = r_geom[0].GetDof(DEFORMATION_GRADIENT_YX, Fpos + 3).EquationId();
        if (dimension == 3){
            rResult[local_index++] = r_geom[0].GetDof(DEFORMATION_GRADIENT_ZZ, Fpos + 4).EquationId();
            rResult[local_index++] = r_geom[0].GetDof(DEFORMATION_GRADIENT_YZ, Fpos + 5).EquationId();
            rResult[local_index++] = r_geom[0].GetDof(DEFORMATION_GRADIENT_ZY, Fpos + 6).EquationId();
            rResult[local_index++] = r_geom[0].GetDof(DEFORMATION_GRADIENT_XZ, Fpos + 7).EquationId();
            rResult[local_index++] = r_geom[0].GetDof(DEFORMATION_GRADIENT_ZX, Fpos + 8).EquationId();
        }
    }

    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedvFParticle<TKernelType>::GetDofList(
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
    
    IndexType local_index = 0;

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        const auto& r_geom = r_neighbours[i]->GetGeometry();

        rElementalDofList[local_index++] = r_geom[0].pGetDof(DISPLACEMENT_X);
        rElementalDofList[local_index++] = r_geom[0].pGetDof(DISPLACEMENT_Y);
        if (dimension == 3)
            rElementalDofList[local_index++] = r_geom[0].pGetDof(DISPLACEMENT_Z);

        rElementalDofList[local_index++] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_XX);
        rElementalDofList[local_index++] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_YY);
        rElementalDofList[local_index++] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_XY);
        rElementalDofList[local_index++] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_YX);
        if (dimension == 3){
            rElementalDofList[local_index++] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_ZZ);
            rElementalDofList[local_index++] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_YZ);
            rElementalDofList[local_index++] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_ZY);
            rElementalDofList[local_index++] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_XZ);
            rElementalDofList[local_index++] = r_geom[0].pGetDof(DEFORMATION_GRADIENT_ZX);
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType>
void TotalLagrangianMixedvFParticle<TKernelType>::GetFirstDerivativesVector(VectorType& rValues, int step) const
{
    KRATOS_TRY
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_per_node = dimension + (dimension * dimension);
    const SizeType mat_size = dofs_per_node * number_of_neighbours;

    if (rValues.size() != mat_size) rValues.resize(mat_size, false);

    IndexType local_index = 0; 

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        const auto& r_geom = r_neighbours[i]->GetGeometry();
        const array_1d<double, 3>& acceleration = r_geom[0].FastGetSolutionStepValue(ACCELERATION, step);
        
        rValues[local_index++] = acceleration[0];
        rValues[local_index++] = acceleration[1];
        if (dimension == 3)
            rValues[local_index++] = acceleration[2];
        
        rValues[local_index++] = r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_XX, step);
        rValues[local_index++] = r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_YY, step);
        rValues[local_index++] = r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_XY, step);
        rValues[local_index++] = r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_YX, step);
        if (dimension == 3){
            rValues[local_index++] = r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_ZZ, step);
            rValues[local_index++] = r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_YZ, step);
            rValues[local_index++] = r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_ZY, step);
            rValues[local_index++] = r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_XZ, step);
            rValues[local_index++] = r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT_DOT_ZX, step);
        }
    }
    KRATOS_CATCH("")
}

template<class TKernelType>
int TotalLagrangianMixedvFParticle<TKernelType>::Check(const ProcessInfo& rProcessInfo) const 
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1)<<"Particle element with invalid Id"<< this->Id() << std::endl;
    KRATOS_ERROR_IF(this->GetGeometry().size() != 1) << "Particle element" << this->Id() << "must have exactly 1 node, found " << this->GetGeometry().size() << std::endl;

    const auto& r_prop = GetProperties();

    if (r_prop[CONSTITUTIVE_LAW] != nullptr) mThisConstitutiveLaw->Check( r_prop, GetGeometry(), rProcessInfo );

    return 0; 
    KRATOS_CATCH("")
}
    
template class TotalLagrangianMixedvFParticle<CubicKernel2D>; 
template class TotalLagrangianMixedvFParticle<CubicKernel3D>;

}