// ==============================================================================
/*
 KratosTopologyOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosTopology                        $
//   Last modified by:	  $Author:   daniel.baumgaertner@tum.de $
// 						  $Co-Author: Octaviano Malfavón Farías $
//   Date:                $Date:                    August 2016 $
//   Revision:            $Revision:                        0.0 $
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

