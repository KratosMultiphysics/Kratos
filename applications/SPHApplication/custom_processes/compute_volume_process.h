// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once

#include "processes/process.h"
#include "includes/model_part.h"
#include "sph_application_variables.h"
#include "geometries/geometry.h"
#include "custom_utilities/compute_volume_utilities.h"

/**
 * @class ComputeVolumeProcess
 * @brief This process computes and assign the correct volume to each particle  
 */

namespace Kratos
{
class KRATOS_API(SPH_APPLICATION) ComputeVolumeProcess 
    : public Process
{
public:
    
    using VectorType = Vector;
    using SizeType = size_t;

    KRATOS_CLASS_POINTER_DEFINITION(ComputeVolumeProcess);

    ComputeVolumeProcess(ModelPart& rThisModelPart, Parameters rThisParameters)  
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