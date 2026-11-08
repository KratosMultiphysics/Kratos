// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once

#include "processes/process.h"
#include "includes/model_part.h"
#include "sph_application_variables.h"

/**
 * @class ComputeKernelCorrectionProcess
 * @brief This class assign the values for the stabilizations coefficients, 
 * and stores them into ProcessInfo 
 */

namespace Kratos
{

//class KRATOS_API(SPH_APPLICATION) AssignStabilizationCoefficientsProcess

class AssignStabilizationCoefficientsProcess
    : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(AssignStabilizationCoefficientsProcess);
    
    AssignStabilizationCoefficientsProcess(ModelPart& rThisModelPart, Parameters rThisParameters)  
        : mrThisModelPart(rThisModelPart), mrThisParameters(rThisParameters)
    {
    }

    void Execute() override
    {
        mrThisModelPart.GetProcessInfo().SetValue(PENALIZATION_COEFFICIENT, mrThisParameters["penalization_coeff"].GetDouble());
        mrThisModelPart.GetProcessInfo().SetValue(DISSIPATION_COEFFICIENT, mrThisParameters["dissipation_coeff"].GetDouble());
    }

    void ExecuteInitialize() override
    {
        this->Execute();
    }

protected:

private:
    ModelPart& mrThisModelPart;
    Parameters mrThisParameters;
};


}