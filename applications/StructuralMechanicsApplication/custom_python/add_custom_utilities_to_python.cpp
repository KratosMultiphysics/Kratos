// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/volume_calculation_under_plane_utility.h"
#include "custom_utilities/vtk_output.hpp"

//Utilities
#include "custom_utilities/formfinding_io_utility.h"

namespace Kratos
{
namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<FormfindingIOUtility>(m,"FormfindingIOUtility")
    .def(init<ModelPart&, const Parameters>())
    .def("PrintModelPart",&FormfindingIOUtility::PrintModelPart)
    .def("ReadPrestressData",&FormfindingIOUtility::ReadPrestressData )
    .def("PrintPrestressData",&FormfindingIOUtility::PrintPrestressData )
    ;

     class_<VolumeCalculationUnderPlaneUtility>(m,"VolumeCalculationUnderPlaneUtility")
        .def(init<Vector, double, Vector>())
        //.def("CalculateVolumeEnclosedByClosedSurface", &VolumeCalcUsingSurfaceUtility::CalculateVolumeEnclosedByClosedSurface)
        .def("CalculateVolume", &VolumeCalculationUnderPlaneUtility::CalculateVolume)
        .def("UpdatePositionOfPlaneBasedOnTargetVolume", &VolumeCalculationUnderPlaneUtility::UpdatePositionOfPlaneBasedOnTargetVolume)
        .def("SetPlaneParameters", &VolumeCalculationUnderPlaneUtility::SetPlaneParameters)
        .def("GetIntersectedArea", &VolumeCalculationUnderPlaneUtility::GetIntersectedArea);
        
    class_<VtkOutput>(m,"VtkOutput")
        .def(init<ModelPart &, std::string, Parameters>())
        .def("PrintOutput", &VtkOutput::PrintOutput)
        .def("PrintOutput", &VtkOutput::PrintOutputSubModelPart);
}

}  // namespace Python.

} // Namespace Kratos


     