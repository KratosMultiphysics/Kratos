#pragma once

#include "custom_elements/small_displacement_particle_2D.h"

namespace Kratos
{

template<class TKernelType>
SmallDisplacementParticle2D<TKernelType>::CalculateAll(
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
    const SizeType dofs_per_node = GetDoFsPerNode();
    const SizeType strain_size = mThisConstitutiveLaw->GetStrainSize();
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    // The neighbours of the particle are like the nodes belonging to the element in standard FEM
    const auto& r_neighbours = this->GetValue(NEIGHBOURS);
    const SizeType mat_size = dofs_per_node * r_neighbours.size();
    
    if (CalculateStiffnessMatrixFlag){
        if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) rLHS.resize(mat_size, mat_size, false);
        noalias(rLHS) = ZeroMatrix(mat_size, mat_size);
    }

    if (CalculateResidualVectorFlag){
        if (rRHS.size() != mat_size) rRHS.resize(mat_size, false);
        noalias(rRHS) = ZeroVector(mat_size);
    }

    VectorType nodal_values(mat_size);
    GetNodalValuesVector(nodal_values);

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
    VectorType f_local(domain_size), local_nodal_values(dofs_per_node);

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
            array_1d<double, 2> body_force = GetLocalBodyForces();
            f_local = volume_gp * volume_ip * body_force * ip_kernel; 
        }

        noalias(project(rRHS, range(dofs_per_node * ip_index, dofs_per_node * (ip_index + 1)))) += f_local;

        VectorType local_nodal_values(dofs_per_node);
        GlobalSizeVector(local_nodal_values, nodal_values, ip_index);
        strain_vector.clear();
        CalculateStrainVector(strain_vector, ip_dkernel, local_nodal_values, volume_ip);

        mThisConstitutiveLaw->CalculateMaterialResponseCauchy(cl_values);
        const VectorType& r_stress_vector = cl_values.GetStressVector();
        const MatrixType& r_constitutive_matrix = cl_values.GetConstitutiveMatrix();

        // Local system computation and assembly
        GetShapeFunctionDerivatives(B_a, ip_dkernel, volume_ip); 

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
            GetShapeFunctionDerivatives(B_b, jp_dkernel, volume_jp);

            f_local = prod(trans(B_b), r_stress_vector);
            f_local *= volume_gp;

            temp = prod(r_constitutive_matrix, B_a);
            K_loc_temp = prod(trans(B_b), temp);
            K_loc_temp *= volume_gp;
                
            if (CalculateStiffnessMatrixFlag){
                noalias(project(rLHS, range(dofs_per_node * jp_index, dofs_per_node * (jp_index + 1)), range(dofs_per_node * ip_index, dofs_per_node * (ip_index + 1)))) += K_loc_temp;
            }
            
            if (CalculateResidualVectorFlag){
                noalias(project(rRHS, range(dofs_per_node * jp_index, dofs_per_node * (jp_index + 1)))) -= f_local;
            }

            jp_index++;
        }
        ip_index++;
    }
    KRATOS_CATCH("")
}

template class SmallDisplacementParticle2D<CubicKernel2D>;


}
