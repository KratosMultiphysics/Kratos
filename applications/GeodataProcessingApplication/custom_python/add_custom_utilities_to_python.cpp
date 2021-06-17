//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nicola Germano
//                   Simon Wenczowski
//


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_utilities/hello.h"

#include "custom_utilities/cleaning_utilities.h"
#include "custom_utilities/extrusion_height_utilities.h"
#include "custom_utilities/fill_cfd_modelpart_utilities.h"
#include "custom_utilities/building_utilities.h"
// #include "custom_utilities/nicola.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

/*     typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType; */

    // for testing --- delete ---
    py::class_<Hello>(m,"Hello")
        .def(py::init<>())
        .def("Greet",&Hello::Greet)
        .def("test_nicola", &Hello::test_nicola)
        .def("CheckIfInternal", &Hello::CheckIfInternal)
        ;

    // to set the extrusion height of the ground terrain
    py::class_<ExtrusionHeightUtilities>(m,"ExtrusionHeightUtilities")
        .def(py::init< ModelPart& >())
        .def("SetExtrusionHeight", &ExtrusionHeightUtilities::SetExtrusionHeight)
        .def("SmoothExtrusionHeight", &ExtrusionHeightUtilities::SmoothExtrusionHeight)
        ;

    // to keep the model part in a good quality state
    py::class_<CleaningUtilities>(m,"CleaningUtilities")
        .def(py::init< ModelPart& >())
        .def("CleanIsolatedNodes", &CleaningUtilities::CleanIsolatedNodes)
        .def("CleanConditions", &CleaningUtilities::CleanConditions)
        .def("CleanConditionsAngles", &CleaningUtilities::CleanConditionsAngles)
        .def("HardCopyBeforeSurfaceDiscretization", &CleaningUtilities::HardCopyBeforeSurfaceDiscretization)
        .def("HardCopyAfterSurfaceDiscretization", &CleaningUtilities::HardCopyAfterSurfaceDiscretization)
        .def("FillBottom", &CleaningUtilities::FillBottom)  // [NG] CHECK WHETHER TO MOVE IT TO ANOTHER CLASS
        ;
    
    // to fill cfd model part
    py::class_<FillCfdModelpartUtilities>(m,"FillCfdModelpartUtilities")
        .def(py::init< ModelPart& >())
        .def("FillModelPart", &FillCfdModelpartUtilities::FillModelPart)
        .def("FillPartsFluid", &FillCfdModelpartUtilities::FillPartsFluid)
        .def("FillInlet", &FillCfdModelpartUtilities::FillInlet)
        .def("FillOutlet", &FillCfdModelpartUtilities::FillOutlet)
        .def("FillSlip", &FillCfdModelpartUtilities::FillSlip)
        .def("FillNoslip", &FillCfdModelpartUtilities::FillNoslip)
        ;
    

    // some operations on buildings
    py::class_<BuildingUtilities>(m,"BuildingUtilities")
        .def(py::init< ModelPart& >())
        .def("CheckOverlapElement", &BuildingUtilities::CheckOverlapElement)
        .def("DeleteNotValidElements", &BuildingUtilities::DeleteNotValidElements)
        .def("CheckIfInternal", &BuildingUtilities::CheckIfInternal)
        // .def("SplitBuilding", &BuildingUtilities::SplitBuilding)
        ;
}

} // namespace Python.

} // Namespace Kratos
