//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "custom_processes/output_quadrature_domain_process.h"
#include "custom_processes/output_eigen_values_process.h"
#include "custom_processes/nitsche_stabilization_model_part_process.h"
#include "custom_processes/map_nurbs_volume_results_to_embedded_geometry_process.h"
#include "custom_processes/assign_integration_points_to_background_elements_process.h"

#include "iga_application_variables.h"


namespace Kratos {
namespace Python {

void AddCustomProcessesToPython(
    pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<OutputQuadratureDomainProcess, OutputQuadratureDomainProcess::Pointer, Process>(m, "OutputQuadratureDomainProcess")
        .def(py::init<Model&, Parameters >())
        ;

    py::class_<OutputEigenValuesProcess, OutputEigenValuesProcess::Pointer, Process>(m, "OutputEigenValuesProcess")
        .def(py::init<Model&, Parameters >())
        ;

    py::class_<NitscheStabilizationModelPartProcess, NitscheStabilizationModelPartProcess::Pointer, Process>(m, "NitscheStabilizationModelPartProcess")
        .def(py::init<ModelPart& >())
        ;

    py::class_<MapNurbsVolumeResultsToEmbeddedGeometryProcess, MapNurbsVolumeResultsToEmbeddedGeometryProcess::Pointer, Process>(m, "MapNurbsVolumeResultsToEmbeddedGeometryProcess")
        .def(py::init<Model&, Parameters >())
        .def("MapVariables", &MapNurbsVolumeResultsToEmbeddedGeometryProcess::MapVariables )
        ;

    py::class_<AssignIntegrationPointsToBackgroundElementsProcess, AssignIntegrationPointsToBackgroundElementsProcess::Pointer, Process>(m, "AssignIntegrationPointsToBackgroundElementsProcess")
        .def(py::init<Model&, Parameters >())
        ;



}

} // namespace Python
} // Namespace Kratos
