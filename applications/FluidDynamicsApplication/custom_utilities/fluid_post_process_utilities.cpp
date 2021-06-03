//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Miguel Angel Celigueta
//                  Aditya Ghantasala
//

#include "fluid_auxiliary_utilities.h"
#include "fluid_post_process_utilities.h"


namespace Kratos {

    double FluidPostProcessUtilities::CalculateFlow(const ModelPart& rModelPart)
    {
        KRATOS_WARNING("CalculateFlow") << "\'CalculateFlow\' is deprecated. Please use \'CalculateFlowRate\' from the \'FluidAuxiliaryUtilities\'." << std::endl;
        return FluidAuxiliaryUtilities::CalculateFlowRate(rModelPart);
    }

} // namespace Kratos.
