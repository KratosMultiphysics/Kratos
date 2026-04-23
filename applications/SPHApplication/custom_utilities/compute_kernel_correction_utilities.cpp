// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once

#include "custom_utilities/compute_kernel_correction_utilities.h"


namespace Kratos
{

void ComputeKernelCorrectionUtilities::ComputeWeightedSums(ModelPart& rThisModelPart)
{
    auto& rElem = rThisModelPart.Elements();
    const double h = rThisModelPart.GetProcessInfo()[SMOOTHING_LENGTH];
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
    const double h = rThisModelPart.GetProcessInfo()[SMOOTHING_LENGTH];
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
        double det_L = 0.0;

        MathUtils<double>::InvertMatrix(gradient_L_aux, inv_gradient_L, det_L);
        IP->SetValue(GRADIENT_CORRECTION, inv_gradient_L);

    }
}

bool ComputeKernelCorrectionUtilities::VerifyKernelCorrection(ModelPart& rThisModelPart, Parameters& rThisParameters)
{
    KRATOS_TRY 

    auto& rElem = rThisModelPart.Elements();
    const double h = rThisModelPart.GetProcessInfo()[SMOOTHING_LENGTH];
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

            ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(IP, kernel[index], dkernel[index]);

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

void ComputeKernelCorrectionUtilities::ApplyKernelCorrection(Element& IP, double& kernel_target)
{
    // kernel becomes corrected kernel
    kernel_target /= IP.GetValue(VW_KERNEL);
} 

void ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(Element& IP, double& kernel_target, VectorType& dkernel_target)
{
    // kernel becomes corrected kernel 
    kernel_target /= IP.GetValue(VW_KERNEL);
    
    // kernel gradient becomes corrected kernel gradient
    VectorType dckernel = dkernel_target / IP.GetValue(VW_KERNEL) - kernel_target * IP.GetValue(VW_DKERNEL) / IP.GetValue(VW_KERNEL);
    noalias(dkernel_target) = prod(IP.GetValue(GRADIENT_CORRECTION), dckernel);
} 

}