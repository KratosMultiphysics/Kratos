#include "custom_processes/compute_kernel_correction_process.h"
#include "custom_utilities/compute_kernel_correction_utilities.h"
#include "sph_application_variables.h" 

namespace Kratos
{
    void ComputeKernelCorrectionProcess::Execute()
    {
        KRATOS_TRY

        const bool flag = mrThisParameters["controls"].GetBool();
        ComputeKernelCorrectionUtilities::ComputeWeightedSums(mrThisModelPart);
        ComputeKernelCorrectionUtilities::ComputeGradientCorrection(mrThisModelPart);
        ComputeKernelCorrectionUtilities::ComputeIntegrationCorrection(mrThisModelPart);
        
        if (flag == true){
            bool correction_flag = ComputeKernelCorrectionUtilities::VerifyKernelCorrection(mrThisModelPart, mrThisParameters);
            KRATOS_INFO("ComputeKernelCorrectionProcess")<<"Kernel correction verification completed. Result: "<<(correction_flag ? "successful" : "failed")<<std::endl;
            
            const bool integration_correction_flag = ComputeKernelCorrectionUtilities::VerifyIntegrationCorrection(mrThisModelPart, mrThisParameters);
            KRATOS_INFO("ComputeKernelCorrectionProcess") << "Integration correction verification completed. Result: "
                << (integration_correction_flag ? "successful" : "failed") << std::endl;
        }

        KRATOS_INFO("ComputeKernelCorrectionProcess")<<"The kernel correction process was executed"<<std::endl;

        KRATOS_CATCH("")
    }

    void ComputeKernelCorrectionProcess::ExecuteInitialize()
    {
        this->Execute();
    }
}
