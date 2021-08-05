// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

// External includes 
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "processes/process.h"

// Application includes
#include "custom_python/add_custom_utilities_to_python.h"

// Utilities
#include "custom_utilities/structure_response_function_utilities.h"
#include "custom_utilities/topology_filtering_utilities.h"
#include "custom_utilities/topology_extractor_utilities.h"
#include "custom_utilities/topology_smoothing_utilities.h"
#include "custom_utilities/topology_updating_utilities.h"
#include "custom_utilities/io_utilities.h"


namespace Kratos
{

namespace Python
{

void AddCustomUtilitiesToPython(pybind11::module& m)

{
	namespace py = pybind11;



	// =============================================================================================================================================
	// Utility Classes
	// =============================================================================================================================================

	py::class_<StructureResponseFunctionUtilities>(m, "StructureResponseFunctionUtilities")
    .def(py::init<ModelPart& >())
	.def("ComputeStrainEnergy", &StructureResponseFunctionUtilities::ComputeStrainEnergy)
	.def("ComputeVolumeFraction", &StructureResponseFunctionUtilities::ComputeVolumeFraction)
	;

	py::class_<TopologyFilteringUtilities >(m, "TopologyFilteringUtilities")
    .def(py::init<ModelPart&, const double, const int>())
	.def("ApplyFilterSensitivity", &TopologyFilteringUtilities::ApplyFilterSensitivity)
	.def("ApplyFilterDensity", &TopologyFilteringUtilities::ApplyFilterDensity)
	;

	py::class_<TopologyExtractorUtilities >(m, "TopologyExtractorUtilities")
    .def(py::init<>())
	.def("ExtractVolumeMesh", &TopologyExtractorUtilities::ExtractVolumeMesh)
	.def("ExtractSurfaceMesh", &TopologyExtractorUtilities::ExtractSurfaceMesh)
	;

	py::class_<TopologySmoothingUtilities >(m, "TopologySmoothingUtilities" )
    .def(py::init<>())
	.def("SmoothMesh", &TopologySmoothingUtilities::SmoothMesh)
	;

	py::class_<TopologyUpdatingUtilities >(m, "TopologyUpdatingUtilities")
	.def(py::init<ModelPart&>())
    .def("UpdateDensitiesUsingOCMethod", &TopologyUpdatingUtilities::UpdateDensitiesUsingOCMethod)
	.def("UpdateDensitiesUsingMMAMethod", &TopologyUpdatingUtilities::UpdateDensitiesUsingMMAMethod)
/* 	.def("InitializeMMA", &TopologyUpdatingUtilities::InitializeMMA) */
	;

	py::class_<IOUtilities >(m, "IOUtilities" )
	.def(py::init<>())
    .def("SaveOptimizationResults", &IOUtilities::SaveOptimizationResults)
	.def("WriteSurfaceAsSTLFile", &IOUtilities::WriteSurfaceAsSTLFile)
	;
}

}  // namespace Python.

} // Namespace Kratos
