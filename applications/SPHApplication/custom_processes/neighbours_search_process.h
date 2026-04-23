// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once

#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_utilities/neighbours_search_utilities.h"
#include "sph_application_variables.h" 
#include "spatial_containers/bins_static.h" // when sph is eulerian you should use bins dynamic (?)

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