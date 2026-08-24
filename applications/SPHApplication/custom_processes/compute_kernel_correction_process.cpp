#include "custom_processes/compute_kernel_correction_process.h"
#include "custom_utilities/compute_kernel_correction_utilities.h"
#include "sph_application_variables.h" 

// TO DO: The integration correction is not working properly, the structure is implemented but the implementation needs to be fixed. Keep the lines commented for now. 

namespace Kratos
{
    void ComputeKernelCorrectionProcess::Execute()
    {
        KRATOS_TRY

        const bool flag = mrThisParameters["controls"].GetBool();
        unsigned int iter;

        ComputeKernelCorrectionUtilities::ComputeWeightedSums(mrThisModelPart);
        ComputeKernelCorrectionUtilities::ComputeGradientCorrection(mrThisModelPart);
        //ComputeKernelCorrectionUtilities::ComputeIntegrationCorrection(mrThisModelPart, mrThisParameters, iter);
        
        if (flag == true){
            bool correction_flag = ComputeKernelCorrectionUtilities::VerifyKernelCorrection(mrThisModelPart, mrThisParameters);
            KRATOS_INFO("ComputeKernelCorrectionProcess")<<"Kernel correction verification completed. Result: "<<(correction_flag ? "successful" : "failed")<<std::endl;
            
            //bool integration_correction_flag = ComputeKernelCorrectionUtilities::VerifyIntegrationCorrection(mrThisModelPart, mrThisParameters);
            //KRATOS_INFO("ComputeKernelCorrectionProcess")<<"The integration correction process was executed: "
            //    <<(integration_correction_flag ? "successful" : "failed")<<" in "<< iter << " iterations" <<std::endl;
        }

        KRATOS_INFO("ComputeKernelCorrectionProcess")<<"The kernel correction process was executed"<<std::endl;

        KRATOS_CATCH("")
    }

    void ComputeKernelCorrectionProcess::ExecuteInitialize()
    {
        this->Execute();
    }
}