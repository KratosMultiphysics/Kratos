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
#include <boost/python.hpp>

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


void AddCustomUtilitiesToPython()
{
	using namespace boost::python;

	// =============================================================================================================================================
	// Utility Classes
	// =============================================================================================================================================

	class_<StructureResponseFunctionUtilities, bases<Process> >("StructureResponseFunctionUtilities", init<ModelPart&>())
    .def("ComputeStrainEnergy", &StructureResponseFunctionUtilities::ComputeStrainEnergy)
	.def("ComputeVolumeFraction", &StructureResponseFunctionUtilities::ComputeVolumeFraction)
	;

	class_<TopologyFilteringUtilities, bases<Process> >("TopologyFilteringUtilities", init<ModelPart&, const double, const int>())
    .def("ApplyFilter", &TopologyFilteringUtilities::ApplyFilter)
	;

	class_<TopologyExtractorUtilities, bases<Process> >("TopologyExtractorUtilities", init<>())
    .def("ExtractVolumeMesh", &TopologyExtractorUtilities::ExtractVolumeMesh)
	.def("ExtractSurfaceMesh", &TopologyExtractorUtilities::ExtractSurfaceMesh)
	;

	class_<TopologySmoothingUtilities, bases<Process> >("TopologySmoothingUtilities", init<>())
    .def("SmoothMesh", &TopologySmoothingUtilities::SmoothMesh)
	;

	class_<TopologyUpdatingUtilities, bases<Process> >("TopologyUpdatingUtilities", init<ModelPart&>())
    .def("UpdateDensitiesUsingOCMethod", &TopologyUpdatingUtilities::UpdateDensitiesUsingOCMethod)
	;

	class_<IOUtilities, bases<Process> >("IOUtilities", init<>())
    .def("SaveOptimizationResults", &IOUtilities::SaveOptimizationResults)
	.def("WriteSurfaceAsSTLFile", &IOUtilities::WriteSurfaceAsSTLFile)
	;
}

}  // namespace Python.

} // Namespace Kratos

