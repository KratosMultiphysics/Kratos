// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once

#include "processes/process.h"
#include "includes/model_part.h"

/**
 * @class ComputeKernelCorrectionProcess
 * @brief This class computes the kernel correction for the SPH method to ensure
 * zeroth and first-order consistency in both domain and boundaries.
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

    /**
     * @details In a Total Lagrangian framework, the kernel corrections are compluted just once at trhe beginning of the simulation.
     */
    void ExecuteInitialize() override;

protected:

private:
    ModelPart& mrThisModelPart;
    Parameters mrThisParameters;
};
}