#pragma once

#include "custom_utilities/compute_kernel_correction_utilities.h"

namespace Kratos
{

void ComputeKernelCorrectionUtilities::ComputeWeightedSums(ModelPart& rThisModelPart)
{
    auto& rElem = rThisModelPart.Elements();
    const SizeType domain_size = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    for (auto IP = rElem.begin(); IP != rElem.end(); ++IP){

        const auto& r_neighbours = IP->GetValue(NEIGHBOURS);

        double vw_kernel_aux = 0.0;
        VectorType vw_dkernel_aux(domain_size);
        noalias(vw_dkernel_aux) = ZeroVector(domain_size);

        std::vector<double> kernel;
        std::vector<Vector> dkernel;

        IP->CalculateOnIntegrationPoints(SPH_KERNEL, kernel, rThisModelPart.GetProcessInfo());
        IP->CalculateOnIntegrationPoints(SPH_KERNEL_GRADIENT, dkernel, rThisModelPart.GetProcessInfo());

        for (IndexType index = 0; index < r_neighbours.size(); index++ ){
            
            const auto& JP = r_neighbours[index];
            const double volume = JP->GetGeometry()[0].GetValue(VOLUME);
            
            vw_kernel_aux += volume * kernel[index];
            noalias(vw_dkernel_aux) += volume * dkernel[index];
        }

        IP->SetValue(VW_KERNEL, vw_kernel_aux);
        IP->SetValue(VW_DKERNEL, vw_dkernel_aux);
    }
}

void ComputeKernelCorrectionUtilities::ComputeGradientCorrection(ModelPart& rThisModelPart)
{   
    auto& rElem = rThisModelPart.Elements();
    const SizeType domain_size = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    for (auto IP = rElem.begin(); IP != rElem.end(); ++IP){

        const auto& r_neighbours = IP->GetValue(NEIGHBOURS);
        const auto& r_geom = IP->GetGeometry();

        Matrix gradient_L_aux(domain_size, domain_size);
        noalias(gradient_L_aux) = ZeroMatrix(domain_size, domain_size);

        std::vector<double> kernel;
        std::vector<Vector> dkernel;

        IP->CalculateOnIntegrationPoints(SPH_KERNEL, kernel, rThisModelPart.GetProcessInfo());
        IP->CalculateOnIntegrationPoints(SPH_KERNEL_GRADIENT, dkernel, rThisModelPart.GetProcessInfo());

        const double vw_kernel = IP->GetValue(VW_KERNEL);
        const Vector& vw_dkernel = IP->GetValue(VW_DKERNEL); 
        const auto& IPcoords = r_geom[0].Coordinates();

        for (IndexType index = 0; index < r_neighbours.size(); index++){

            const auto& JP = r_neighbours[index];
            const auto& r_geom_neigh = JP->GetGeometry();

            Vector X_AB_target(domain_size);
            const auto& JPcoords = r_geom_neigh[0].Coordinates();
            for (IndexType d = 0; d < domain_size; d++){
                X_AB_target[d] = IPcoords[d] - JPcoords[d];
            }
            
            const double volume = r_geom_neigh[0].GetValue(VOLUME);

            Vector dckernel = dkernel[index] / vw_kernel - kernel[index] * vw_dkernel / (vw_kernel * vw_kernel);
            noalias(gradient_L_aux) += volume * outer_prod(dckernel, - X_AB_target);

        }

        Matrix inv_gradient_L(domain_size, domain_size);
        double det_L = MathUtils<double>::Det(gradient_L_aux); 

        MathUtils<double>::InvertMatrix(gradient_L_aux, inv_gradient_L, det_L);
        IP->SetValue(GRADIENT_CORRECTION, inv_gradient_L);

    }
}

void ComputeKernelCorrectionUtilities::ComputeIntegrationCorrection(
    ModelPart& rThisModelPart,
    Parameters& rThisParameters,
    unsigned int& rIterations)
{
    KRATOS_TRY
    /**
     * The integration correction is computed iteratively. As described in the reference this correspond to 
     * compute iteratively a vector of correction coefficients gamma_i for each particle i. The iterative process
     * consists in solving an equation of the following form:
     * gamma_i^(k+1) = A gamma_i^(k) + rhs_i
     */

    auto& rElem = rThisModelPart.Elements();
    const SizeType number_of_particles = rElem.size();
    const SizeType domain_size = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    const double relative_tolerance = 1.0e-08;
    const double absolute_tolerance = 1.0e-10;
    const unsigned int max_iterations = 25000;

    for (auto& r_particle : rElem){
        r_particle.GetGeometry()[0].SetValue(INTEGRATION_CORRECTION_VARIABLE, ZeroVector(domain_size));
    }

    // Initialization of the gamma vector.
    std::vector<Vector> new_correction(number_of_particles, ZeroVector(domain_size));

    // Iterative procedure until convergence.
    for (rIterations = 0; rIterations <= max_iterations; rIterations++){

        SizeType particle_index = 0;
        
        for (auto IP = rElem.begin(); IP != rElem.end(); ++IP){

            const auto& r_neighbours = IP->GetValue(NEIGHBOURS);

            std::vector<double> kernel;
            std::vector<Vector> dkernel;
            IP->CalculateOnIntegrationPoints(SPH_KERNEL, kernel, rThisModelPart.GetProcessInfo());
            IP->CalculateOnIntegrationPoints(SPH_KERNEL_GRADIENT, dkernel, rThisModelPart.GetProcessInfo());

            Vector rhs(domain_size); rhs.clear();

            // Compute the right-hand side rhs_i of the equation.
            for (IndexType j = 0; j < r_neighbours.size(); ++j){

                ApplyKernelGradientCorrectionInverted(*r_neighbours[j], kernel[j], dkernel[j]);

                const auto& r_geom = r_neighbours[j]->GetGeometry();
                const double volume = r_geom[0].GetValue(VOLUME);
                const Vector& r_boundary_normal_area = r_geom[0].GetValue(BOUNDARY_NORMAL_AREA);

                noalias(rhs) += r_boundary_normal_area * kernel[j] - volume * dkernel[j];

                const Vector& r_old_correction = r_geom[0].GetValue(INTEGRATION_CORRECTION_VARIABLE);

                noalias(rhs) += volume * kernel[j] * r_old_correction;
            }

            new_correction[particle_index++].swap(rhs);
        }

        double maximum_change = 0.0; double maximum_correction = 0.0;
        
        particle_index = 0;
        for (auto IP = rElem.begin(); IP != rElem.end(); ++IP){
            
            const Vector& r_old_correction = IP->GetGeometry()[0].GetValue(INTEGRATION_CORRECTION_VARIABLE);

            const double correction_change = norm_2(new_correction[particle_index] - r_old_correction);
            const double correction_norm = norm_2(new_correction[particle_index]);

            maximum_change = std::max(maximum_change, correction_change);
            maximum_correction = std::max(maximum_correction, correction_norm);
            
            ++particle_index;
        }

        particle_index = 0;
        for (auto IP = rElem.begin(); IP != rElem.end(); ++IP){
            IP->GetGeometry()[0].SetValue(INTEGRATION_CORRECTION_VARIABLE, new_correction[particle_index++]);
        }

        if (maximum_change <= absolute_tolerance || maximum_change <= relative_tolerance * maximum_correction){
            break;
        }

    }

    KRATOS_CATCH("")
}

void ComputeKernelCorrectionUtilities::ApplyKernelCorrection(Element& IP, double& kernel_target)
{
    kernel_target /= IP.GetValue(VW_KERNEL);
} 

void ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(
    Element& rIntegrationParticle,
    double& rKernel,
    VectorType& rKernelGradient)
{
    rKernel /= rIntegrationParticle.GetValue(VW_KERNEL);

    const Vector corrected_kernel_gradient = rKernelGradient / rIntegrationParticle.GetValue(VW_KERNEL) - rKernel * rIntegrationParticle.GetValue(VW_DKERNEL) / rIntegrationParticle.GetValue(VW_KERNEL);
    noalias(rKernelGradient) = prod(rIntegrationParticle.GetValue(GRADIENT_CORRECTION), corrected_kernel_gradient);
}

void ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrectionInverted(
    Element& rNeighbouringParticle,
    double& rKernel,
    VectorType& rKernelGradient)
{
    rKernelGradient *= -1.0;
    ApplyKernelGradientCorrection(rNeighbouringParticle, rKernel, rKernelGradient); 
}

void ComputeKernelCorrectionUtilities::ApplyIntegrationCorrection(
    Element& rIntegrationParticle,
    double& rKernel,
    VectorType& rKernelGradient,
    bool IsParticleItself)
{
    const Vector& r_gamma = rIntegrationParticle.GetGeometry()[0].GetValue(INTEGRATION_CORRECTION_VARIABLE);
    const double volume = rIntegrationParticle.GetGeometry()[0].GetValue(VOLUME);

    noalias(rKernelGradient) -= r_gamma * rKernel;
    if (IsParticleItself) noalias(rKernelGradient) += r_gamma / volume;
}

bool ComputeKernelCorrectionUtilities::VerifyKernelCorrection(ModelPart& rThisModelPart, Parameters& rThisParameters)
{
    KRATOS_TRY 

    auto& rElem = rThisModelPart.Elements();
    const SizeType domain_size = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
    
    const double tol = rThisParameters["tol"].GetDouble();

    // Initializing the controls 
    Matrix I(domain_size, domain_size);
    noalias(I) = IdentityMatrix(domain_size);

    for (auto& IP : rElem){

        const auto& r_neighbours = IP.GetValue(NEIGHBOURS);
        const auto& IPcoords = IP.GetGeometry()[0].Coordinates();

        // Initializing the controls 
        double control1 = 0.0;
        Matrix control3(domain_size, domain_size);
        noalias(control3) = ZeroMatrix(domain_size, domain_size);
        Vector control2(domain_size);
        noalias(control2) = ZeroVector(domain_size);

        std::vector<double> kernel;
        std::vector<Vector> dkernel;

        IP.CalculateOnIntegrationPoints(SPH_KERNEL, kernel, rThisModelPart.GetProcessInfo());
        IP.CalculateOnIntegrationPoints(SPH_KERNEL_GRADIENT, dkernel, rThisModelPart.GetProcessInfo());

        for (IndexType index = 0; index < r_neighbours.size(); index++){
            
            auto& JP = r_neighbours[index];

            Vector X_AB_target(domain_size);
            const auto& JPcoords = JP->GetGeometry()[0].Coordinates();
            for (IndexType d = 0; d < domain_size; d++){
                X_AB_target[d] = IPcoords[d] - JPcoords[d];
            }

            ApplyKernelGradientCorrection(IP, kernel[index], dkernel[index]);

            const double volume = JP->GetGeometry()[0].GetValue(VOLUME);

            control1 += volume * kernel[index];
            noalias(control2) += volume * dkernel[index];
            noalias(control3) += volume * outer_prod(dkernel[index], - X_AB_target);
            
        }

        if (std::abs(control1 - 1.0) > tol){
            KRATOS_WARNING("ComputeKernelCorrections")<<"Zeroth order check failed"<<std::endl;
            return false;
        }

        if (norm_2(control2) > tol){
            KRATOS_WARNING("ComputeKernelCorrections")<<"First order gradient check failed"<<std::endl;
            return false;
        }

        if(norm_frobenius(control3 - I) > tol){
            KRATOS_WARNING("ComputeKernelCorrections")<<"Linear completeness check failed"<<std::endl;
            return false;
        }

    }

    return true;

    KRATOS_CATCH("")
}

bool ComputeKernelCorrectionUtilities::VerifyIntegrationCorrection(ModelPart& rThisModelPart, Parameters& rThisParameters)
{
    KRATOS_TRY
    auto& r_elements = rThisModelPart.Elements();
    const SizeType domain_size = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double tol = rThisParameters["tol"].GetDouble();

    bool is_particle_itself; 

    for (auto rElem = r_elements.begin(); rElem != r_elements.end(); ++rElem){

        const auto& r_neighbours = rElem->GetValue(NEIGHBOURS);
        Vector residual = ZeroVector(domain_size);

        std::vector<double> kernel;
        std::vector<Vector> kernel_gradient;
        rElem->CalculateOnIntegrationPoints(SPH_KERNEL, kernel, rThisModelPart.GetProcessInfo());
        rElem->CalculateOnIntegrationPoints(SPH_KERNEL_GRADIENT, kernel_gradient, rThisModelPart.GetProcessInfo());

        for (IndexType neighbour_index = 0; neighbour_index < r_neighbours.size(); ++neighbour_index) {
            Element& r_neighbour = *r_neighbours[neighbour_index];
            const auto& r_neighbour_node = r_neighbour.GetGeometry()[0];
            const double volume = r_neighbour_node.GetValue(VOLUME);

            is_particle_itself = (rElem->Id() == r_neighbour.Id());

            ApplyKernelGradientCorrectionInverted(r_neighbour, kernel[neighbour_index], kernel_gradient[neighbour_index]);
            ApplyIntegrationCorrection(r_neighbour, kernel[neighbour_index], kernel_gradient[neighbour_index], is_particle_itself);

            noalias(residual) += volume * kernel_gradient[neighbour_index] - r_neighbour_node.GetValue(BOUNDARY_NORMAL_AREA) * kernel[neighbour_index];
        }

        KRATOS_WATCH(norm_2(residual));

        if (norm_2(residual) > tol){
            KRATOS_WARNING("ComputeKernelCorrections")<<"Integration correction check failed"<<std::endl;
            return false;
        }
    }

    return true;

    KRATOS_CATCH("")
}

}