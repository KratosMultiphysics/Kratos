// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Manuel Messmer
//

// System includes

// External includes


// Project includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/pertube_geometry.h"
#include "custom_processes/pertube_geometry_small_correlation_length.h"


namespace Kratos {
namespace Python {

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// Processes
    py::class_<PertubeGeometryProcess, PertubeGeometryProcess::Pointer, Process>(m,"PertubeGeometryProcess")
        .def(py::init<ModelPart&,double>())
        .def("CreateEigenvectors", &PertubeGeometryProcess::CreateEigenvectors)
        .def("AssembleEigenvectors", &PertubeGeometryProcess::AssembleEigenvectors)
        .def("Average", &PertubeGeometryProcess::Average)
        ;
    
    py::class_<PertubeGeometrySmallCorrelationLength, PertubeGeometrySmallCorrelationLength::Pointer, Process>(m,"PertubeGeometrySmallCorrelationLength")
        .def(py::init<ModelPart&,double>()) 
        .def("CreateEigenvectors", &PertubeGeometrySmallCorrelationLength::CreateEigenvectors)
        .def("AssembleEigenvectors", &PertubeGeometrySmallCorrelationLength::AssembleEigenvectors)
        .def("Average", &PertubeGeometrySmallCorrelationLength::Average)
        ;

}

}  // namespace Python.
} // Namespace Kratos

