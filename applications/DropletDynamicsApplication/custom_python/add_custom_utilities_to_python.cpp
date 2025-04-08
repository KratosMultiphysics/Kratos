//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//


// System includes

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>  // This is needed for py::array_t

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_utilities/contact_angle_evaluator.h"

#include "custom_utilities/intersection_points_utility.h"  // Include for IntersectionPointsUtility
#include "droplet_dynamics_application_variables.h"  // Include for IntersectionPointData


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    py::class_<ContactAngleEvaluator, ContactAngleEvaluator::Pointer, Process>(m,"ContactAngleEvaluatorProcess")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters& >());

    // Register intersection points data and utility
    py::class_<IntersectionPointData>(m, "IntersectionPointData")
        .def(py::init<>())
        .def_readwrite("elementId", &IntersectionPointData::elementId)
        .def_readwrite("pointId", &IntersectionPointData::pointId)
        .def_property("coordinates",
            [](IntersectionPointData& self) { return py::array_t<double>(3, &self.coordinates[0]); },
            [](IntersectionPointData& self, py::array_t<double> arr) {
                for (int i = 0; i < 3; i++) self.coordinates[i] = arr.at(i);
            });
    
    py::class_<KratosDropletDynamics::IntersectionPointsUtility>(m, "IntersectionPointsUtility")
        .def_static("CollectElementIntersectionPoints", &KratosDropletDynamics::IntersectionPointsUtility::CollectElementIntersectionPoints)
        .def_static("ClearIntersectionPoints", &KratosDropletDynamics::IntersectionPointsUtility::ClearIntersectionPoints)
        .def_static("GetIntersectionPoints", &KratosDropletDynamics::IntersectionPointsUtility::GetIntersectionPoints, py::return_value_policy::reference)
        .def_static("SaveIntersectionPointsToFile", &KratosDropletDynamics::IntersectionPointsUtility::SaveIntersectionPointsToFile)
        // .def_static("AddIntersectionPoint", &KratosDropletDynamics::IntersectionPointsUtility::AddIntersectionPoint)
        .def_static("ExtractIntersectionPointsFromSplitter", &KratosDropletDynamics::IntersectionPointsUtility::ExtractIntersectionPointsFromSplitter)
        .def_static("DiagnosticOutput", &KratosDropletDynamics::IntersectionPointsUtility::DiagnosticOutput)
        .def_static("ProcessIntersectionPointsAndFitCurves", &KratosDropletDynamics::IntersectionPointsUtility::ProcessIntersectionPointsAndFitCurves)
        .def_static("ProcessIntersectionPointsAndFitCurvesparabola", &KratosDropletDynamics::IntersectionPointsUtility::ProcessIntersectionPointsAndFitCurvesparabola);

}

} // namespace Python.
} // Namespace Kratos
