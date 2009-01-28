//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-08-21 14:11:10 $
//   Revision:            $Revision: 1.3 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/FSI_utils.h"
#include "custom_utilities/aitken_utils.h"

namespace Kratos
{
	
namespace Python
{
	
  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;

	  class_<FSIUtils>("FSIUtils", init<>())
//		.def("FSIUtils",&FSIUtils::GenerateCouplingElements)
		.def("CheckPressureConvergence",&FSIUtils::CheckPressureConvergence)
		.def("StructuralPressurePrediction",&FSIUtils::StructuralPressurePrediction)
		;

	  class_<AitkenUtils>("AitkenUtils", init<>())
		.def("ComputeAitkenFactor",&AitkenUtils::ComputeAitkenFactor)
		.def("ComputeRelaxedDisplacement",&AitkenUtils::ComputeRelaxedDisplacement)
		;
  }
	
}  // namespace Python.

} // Namespace Kratos

