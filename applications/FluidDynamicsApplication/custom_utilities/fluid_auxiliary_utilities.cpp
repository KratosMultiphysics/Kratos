//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Suneth Warnakulasuriya
//

// System includes


// External includes


// Project includes
#include "utilities/parallel_utilities.h"

// Application includes
#include "fluid_auxiliary_utilities.h"

namespace Kratos
{

double FluidAuxiliaryUtilities::CalculateFluidVolume(ModelPart& rModelPart)
{
    KRATOS_ERROR_IF(rModelPart.NumberOfElements() == 0) << "There are no elements in the provided model part. Fluid volume cannot be computed." << std::endl;

    double fluid_volume = 0.0;
    fluid_volume = block_for_each<SumReduction<double>>(rModelPart.Elements(), [](Element& rElement){
        const auto& r_geom = rElement.GetGeometry();
        return r_geom.DomainSize();
    });
    rModelPart.GetCommunicator().GetDataCommunicator().SumAll(fluid_volume);

    return fluid_volume;
}

double FluidAuxiliaryUtilities::CalculateFluidPositiveVolume(ModelPart& rModelPart)
{
    // KRATOS_ERROR_IF(rModelPart.NumberOfElements() == 0) << "There are no elements in the provided model part. Fluid volume cannot be computed." << std::endl;

    // double fluid_volume = 0.0;
    // fluid_volume = block_for_each<SumReduction<double>>(rModelPart.Elements(), [](Element& rElement){
    //     const auto& r_geom = rElement.GetGeometry();
    //     return r_geom.DomainSize();
    // });
    // rModelPart.GetCommunicator().GetDataCommunicator().SumAll(fluid_volume);

    // return fluid_volume;
}

double FluidAuxiliaryUtilities::CalculateFluidNegativeVolume(ModelPart& rModelPart)
{
    // KRATOS_ERROR_IF(rModelPart.NumberOfElements() == 0) << "There are no elements in the provided model part. Fluid volume cannot be computed." << std::endl;

    // double fluid_volume = 0.0;
    // fluid_volume = block_for_each<SumReduction<double>>(rModelPart.Elements(), [](Element& rElement){
    //     const auto& r_geom = rElement.GetGeometry();
    //     return r_geom.DomainSize();
    // });
    // rModelPart.GetCommunicator().GetDataCommunicator().SumAll(fluid_volume);

    // return fluid_volume;
}

} // namespace Kratos
