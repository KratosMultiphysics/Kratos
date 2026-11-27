// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once

#include "processes/process.h"
#include "includes/model_part.h"
#include "sph_application_variables.h"

/**
 * @class AssignStabilizationCoefficientsProcess
 * @brief This class assign the values for the stabilizations coefficients and stores them into ProcessInfo 
 */

namespace Kratos
{

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
        /**
         * @class TotaLagrangianDisplacementParticle and SmallDisplacementParticle 
         * @brief This coefficient is used to penalize relative current position of the particles that are inconsistent with the local deformation gradient. 
         */
        mrThisModelPart.GetProcessInfo().SetValue(PENALIZATION_COEFFICIENT, mrThisParameters["penalization_coeff"].GetDouble());

        /**
         * @class TotaLagrangianDisplacementParticle, SmallDisplacementParticle and TotaLagrangianMixedStrainParticle
         * @brief This coefficient is used to penalize shocks or strong discontinuities by introducing an additional localised dissipation. 
         */
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