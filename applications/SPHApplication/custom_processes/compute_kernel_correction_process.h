// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once

#include "processes/process.h"
#include "includes/model_part.h"

/**
 * @class 
 * @brief 
 */

namespace Kratos
{
class KRATOS_API(SPH_APPLICATION) ComputeKernelCorrectionProcess
    : public Process
{
public:
    
    KRATOS_CLASS_POINTER_DEFINITION(ComputeKernelCorrectionProcess);

    ComputeKernelCorrectionProcess(ModelPart& rThisModelPart, Parameters rThisParameters)  
        : mrThisModelPart(rThisModelPart), mrThisParameters(rThisParameters)
    {
    }

    void Execute() override;

    void ExecuteInitialize() override;

protected:

private:
    ModelPart& mrThisModelPart;
    Parameters mrThisParameters;
};
}