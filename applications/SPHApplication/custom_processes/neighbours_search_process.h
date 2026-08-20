//  ____  ____  _   _                   _ _           _   _             
// / ___||  _ \| | | | __ _ _ __  _ __ | (_) ___ __ _| |_(_) ___  _ __  
// \___ \| |_) | |_| |/ _` | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_ \ 
//  ___) |  __/|  _  | (_| | |_) | |_) | | | (_| (_| | |_| | (_) | | | |
// |____/|_|   |_| |_|\__,_| .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                         |_|   |_|                                    

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once

#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_utilities/neighbours_search_utilities.h"
#include "sph_application_variables.h" 

/**
 * @class NeighboursSearchProcess
 * @brief This process computes the neighbours of particle inside a radius twice the smoothing length 
 */

namespace Kratos
{
class KRATOS_API(SPH_APPLICATION) NeighboursSearchProcess
    : public Process
{
public:

    using SizeType = std::size_t;
    using ElementType = Element::Pointer;

    KRATOS_CLASS_POINTER_DEFINITION(NeighboursSearchProcess);

    NeighboursSearchProcess(ModelPart& rThisModelPart, Parameters rThisParameters)  
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
