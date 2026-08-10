// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#include "custom_processes/compute_kernel_correction_process.h"
#include "custom_utilities/compute_kernel_correction_utilities.h"
#include "sph_application_variables.h" 


namespace Kratos
{
    void ComputeKernelCorrectionProcess::Execute()
    {
        KRATOS_TRY

        bool flag = mrThisParameters["controls"].GetBool();
        unsigned int iter;

        ComputeKernelCorrectionUtilities::ComputeWeightedSums(mrThisModelPart);
        ComputeKernelCorrectionUtilities::ComputeGradientCorrection(mrThisModelPart);
        //ComputeKernelCorrectionUtilities::ComputeIntegrationCorrection(mrThisModelPart, mrThisParameters, iter);
        
        if (flag == true){
            bool correction_flag = ComputeKernelCorrectionUtilities::VerifyKernelCorrection(mrThisModelPart, mrThisParameters);
            KRATOS_INFO("ComputeKernelCorrectionProcess")<<"Performing verification of kernel and kernel gradient corrections..."<<std::endl
                <<"Kernel correction verification completed. Result: "<<(correction_flag ? "successful" : "failed")<<std::endl;
            KRATOS_INFO("ComputeKernelCorrectionProcess")<<"The integration correction process was executed"<<std::endl
                <<"Number of iterations for convergence: "<< iter <<std::endl;
        }

        KRATOS_INFO("ComputeKernelCorrectionProcess")<<"The kernel correction process was executed"<<std::endl;

        KRATOS_CATCH("")
    }

    void ComputeKernelCorrectionProcess::ExecuteInitialize()
    {
        this->Execute();
    }
}