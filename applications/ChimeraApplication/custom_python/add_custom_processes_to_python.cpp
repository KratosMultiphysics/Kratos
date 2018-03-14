//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
// ==============================================================================

// System includes

// External includes
#include <boost/python.hpp>
// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/custom_hole_cutting_process.h"
#include "custom_processes/apply_chimera_process.h"
#include "custom_processes/custom_calculate_signed_distance_process.h"
#include "custom_processes/calculate_signed_distance_to_2d_condition_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"
#include "custom_processes/rotate_region_process.h"
//#include "custom_processes/apply_multi_point_constraints_process.h"
namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython()
{
	using namespace boost::python;
		
  
	/*
	 * ApplyChimeraProcess for 2d and 3d
	 */
	class_<ApplyChimeraProcess<2>,bases<Process> >("ApplyChimeraProcess2d", init< ModelPart&, Parameters >())
			.def("ApplyMpcConstraint", &ApplyChimeraProcess<2>::ApplyMpcConstraint)
			.def("FormulateChimera2D", &ApplyChimeraProcess<2>::FormulateChimera)
			.def("SetOverlapDistance",&ApplyChimeraProcess<2>::SetOverlapDistance)
			.def("CalculateNodalAreaAndNodalMass",&ApplyChimeraProcess<2>::CalculateNodalAreaAndNodalMass)
			.def("SetType",&ApplyChimeraProcess<2>::SetType);


	class_<ApplyChimeraProcess<3>,bases<Process> >("ApplyChimeraProcess3d", init< ModelPart&, Parameters >())
			.def("ApplyMpcConstraint", &ApplyChimeraProcess<3>::ApplyMpcConstraint)		
			.def("FormulateChimera3D", &ApplyChimeraProcess<3>::FormulateChimera)
			.def("SetOverlapDistance",&ApplyChimeraProcess<3>::SetOverlapDistance)
			.def("CalculateNodalAreaAndNodalMass",&ApplyChimeraProcess<3>::CalculateNodalAreaAndNodalMass)

			.def("SetType",&ApplyChimeraProcess<3>::SetType);

    class_< RotateRegionProcess, bases<Process>, boost::noncopyable >
    ("RotateRegionProcess",init<ModelPart&, Parameters& >())
	.def("SetCentreOfRotation", &RotateRegionProcess::SetCentreOfRotation);

	class_<CustomCalculateSignedDistanceProcess<2>>("SignedDistanceProcess2d", init< >())
			.def("CalculateSignedDistance", &CustomCalculateSignedDistanceProcess<2>::CalculateSignedDistance);
	


	class_<CustomCalculateSignedDistanceProcess<3>>("SignedDistanceProcess3d", init<>())
			.def("CalculateSignedDistance", &CustomCalculateSignedDistanceProcess<3>::CalculateSignedDistance)	;	



}



} // namespace Python.

} // Namespace Kratos

