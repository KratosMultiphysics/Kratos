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
#include "custom_utilities/intersection_points_utility.h" // Include for IntersectionPointsUtility
/////////////
#include "custom_utilities/intersection_points_utility.cpp" // Include for IntersectionPointsUtility
/////////////
#include "custom_utilities/intersection_points_container.h" // For IntersectionPointData and InterfaceAverageData
#include "droplet_dynamics_application_variables.h"

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
        .def_static("ProcessIntersectionPointsAndFitCurvesparabola", &KratosDropletDynamics::IntersectionPointsUtility::ProcessIntersectionPointsAndFitCurvesparabola)
        .def_static("ProcessIntersectionPointsAndFitGeneralConic", &KratosDropletDynamics::IntersectionPointsUtility::ProcessIntersectionPointsAndFitGeneralConic);
    
    // Register interface averages data and utility
    py::class_<InterfaceAverageData>(m, "InterfaceAverageData")
        .def(py::init<>())
        .def_readwrite("elementId", &InterfaceAverageData::elementId)
        .def_readwrite("numberOfPoints", &InterfaceAverageData::numberOfPoints)
        .def_readwrite("interfaceArea", &InterfaceAverageData::interfaceArea)
        .def_property("averageCoordinates",
            [](InterfaceAverageData& self) { 
                return py::array_t<double>(3, &self.averageCoordinates[0]); 
            },
            [](InterfaceAverageData& self, py::array_t<double> arr) {
                for (int i = 0; i < 3; i++) 
                    self.averageCoordinates[i] = arr.at(i);
            })
        .def_property("averageNormal",
            [](InterfaceAverageData& self) { 
                return py::array_t<double>(3, &self.averageNormal[0]); 
            },
            [](InterfaceAverageData& self, py::array_t<double> arr) {
                for (int i = 0; i < 3; i++) 
                    self.averageNormal[i] = arr.at(i);
            });
    
    py::class_<KratosDropletDynamics::InterfaceAveragesUtility>(m, "InterfaceAveragesUtility")
        .def_static("CollectElementInterfaceAverages", &KratosDropletDynamics::InterfaceAveragesUtility::CollectElementInterfaceAverages)
        .def_static("ComputeModelPartInterfaceAverages", &KratosDropletDynamics::InterfaceAveragesUtility::ComputeModelPartInterfaceAverages)
        .def_static("ClearInterfaceAverages", &KratosDropletDynamics::InterfaceAveragesUtility::ClearInterfaceAverages)
        .def_static("GetInterfaceAverages", &KratosDropletDynamics::InterfaceAveragesUtility::GetInterfaceAverages, py::return_value_policy::reference)
        .def_static("SaveInterfaceAveragesToFile", &KratosDropletDynamics::InterfaceAveragesUtility::SaveInterfaceAveragesToFile)
        .def_static("DiagnosticOutput", &KratosDropletDynamics::InterfaceAveragesUtility::DiagnosticOutput);
        // .def_static("ApplyInterfaceAveragesToModelPart", &KratosDropletDynamics::InterfaceAveragesUtility::ApplyInterfaceAveragesToModelPart);
        
    // Register the normal fitting functions
    // m.def("FitLinearNormal", &KratosDropletDynamics::FitLinearNormal, 
    //       py::arg("rModelPart"), 
    //       py::arg("rInterfaceAverages"), 
    //       py::arg("ElementId"),
    //       "Fit a normal vector linearly based on 3 points (element + 2 neighbors)");
m.def("FitLinearNormal", &KratosDropletDynamics::FitLinearNormal,
    py::arg("rModelPart"),
    py::arg("rInterfaceAverages"),
    py::arg("ElementId"),
    py::arg("a0"),
    py::arg("a1"),
    py::arg("a2"),
    py::arg("b0"),
    py::arg("b1"),
    py::arg("b2"),
    "Fit a normal vector using exactly three points (target element and two neighbors)");

          
    // m.def("ApplyFittedNormalsToModelPart", &KratosDropletDynamics::ApplyFittedNormalsToModelPart,
    //       py::arg("rModelPart"),
    //       py::arg("rInterfaceAverages"),
    //       py::arg("StoreOriginalNormal") = true,
    //       "Apply fitted normals to all interface elements in the model part");
          
    m.def("SaveFittedNormalsToFile", &KratosDropletDynamics::SaveFittedNormalsToFile,
          py::arg("rModelPart"),
          py::arg("rInterfaceAverages"),
          py::arg("Filename"),
          "Save fitted normals to a file for visualization and debugging");
}

} // namespace Python.
} // Namespace Kratos
