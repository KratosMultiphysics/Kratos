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

        ComputeKernelCorrectionUtilities::ComputeWeightedSums(mrThisModelPart);

        ComputeKernelCorrectionUtilities::ComputeGradientCorrection(mrThisModelPart);

        ComputeKernelCorrectionUtilities::ComputeIntegrationCorrectionVector(mrThisModelPart);

        KRATOS_CATCH("")
    }

    void ComputeKernelCorrectionProcess::ExecuteInitialize()
    {
        this->Execute();

        if (mrThisParameters["controls"].GetBool()){
            bool correction_flag = ComputeKernelCorrectionUtilities::VerifyKernelCorrection(mrThisModelPart, mrThisParameters);
            KRATOS_INFO("ComputeKernelCorrectionProcess")<<"Performing verification of kernel and kernel gradient corrections..."<<std::endl
                <<"Kernel correction verification completed. Result: "<<(correction_flag ? "successful" : "failed")<<std::endl;
        }
        
        KRATOS_INFO("ComputeKernelCorrectionProcess")<<"The kernel correction process was executed"<<std::endl;
    }
}