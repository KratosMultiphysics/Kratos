//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah
//
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

// Application includes
#include "custom_utilities/parabolic_profile_utilities.h"
#include "custom_utilities/wss_statistics_utilities.h"

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // WSS statistics utilities
    py::class_<WssStatisticsUtilities>(m,"WssStatisticsUtilities")
        .def_static("InitializeWSSVariables", &WssStatisticsUtilities::InitializeWSSVariables)
        .def_static("CalculateWSS", &WssStatisticsUtilities::CalculateWSS)
        .def_static("CalculateOSI", &WssStatisticsUtilities::CalculateOSI)
        ;

    // Parabolic profile utilities
    py::class_<ParabolicProfileUtilities>(m,"ParabolicProfileUtilities")
        .def_static("CalculateWallParallelDistance", &ParabolicProfileUtilities::CalculateWallParallelDistance)
        .def_static("ImposeParabolicInlet", [](ModelPart& rModelPart, const double rMaxParabolaValue){
            ParabolicProfileUtilities::ImposeParabolicInlet(rModelPart, rMaxParabolaValue);})
        .def_static("ImposeParabolicInlet", [](ModelPart& rModelPart, const GenericFunctionUtility::Pointer rMaxParabolaValue){
            ParabolicProfileUtilities::ImposeParabolicInlet(rModelPart, rMaxParabolaValue);})
        .def_static("FreeParabolicInlet", &ParabolicProfileUtilities::FreeParabolicInlet)
        ;
}

}  // namespace Python.

} // Namespace Kratos
