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
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/volume_calculation_under_plane_utility.h"
#include "custom_utilities/vtk_output.hpp"

//Utilities

namespace Kratos
{
namespace Python
{

void AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    //     typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    //     typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    //     typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
     class_<VolumeCalculationUnderPlaneUtility, boost::noncopyable>("VolumeCalculationUnderPlaneUtility", init<ModelPart &>())
        //.def("CalculateVolumeEnclosedByClosedSurface", &VolumeCalcUsingSurfaceUtility::CalculateVolumeEnclosedByClosedSurface)
        .def("CalculateVolume", &VolumeCalculationUnderPlaneUtility::CalculateVolume)
        .def("UpdatePositionOfPlaneBasedOnTargetVolume", &VolumeCalculationUnderPlaneUtility::UpdatePositionOfPlaneBasedOnTargetVolume)
        .def("GetIntersectedArea", &VolumeCalculationUnderPlaneUtility::GetIntersectedArea); 

    class_<VtkOutput, boost::noncopyable>("VtkOutput", init<ModelPart &, std::string, Parameters>())
        .def("PrintOutput", &VtkOutput::PrintOutput)
        .def("PrintOutput", &VtkOutput::PrintOutputSubModelPart);
}

} // namespace Python.

} // Namespace Kratos
