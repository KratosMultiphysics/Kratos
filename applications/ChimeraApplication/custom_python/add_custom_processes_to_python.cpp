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

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"


//#include "custom_processes/spalart_allmaras_turbulence_model_for_chimera.h"
#include "custom_processes/custom_hole_cutting_process.h"
#include "custom_processes/apply_chimera_process_monolithic.h"
#include "custom_processes/apply_chimera_process_FractionalStep.h"
#include "custom_processes/custom_calculate_signed_distance_process.h"
#include "custom_processes/calculate_signed_distance_to_2d_condition_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"
#include "custom_processes/rotate_region_process.h"
namespace Kratos
{

namespace Python
{
using namespace pybind11;

void AddCustomProcessesToPython(pybind11::module& m)
{

	class_<ApplyChimeraProcessMonolithic<2>, ApplyChimeraProcessMonolithic<2>::Pointer, Process>(m, "ApplyChimeraProcessMonolithic2d")
		.def(init< ModelPart&, Parameters >())
		.def("ApplyMpcConstraint", &ApplyChimeraProcessMonolithic<2>::ApplyMpcConstraint)
		.def("FormulateChimera2D", &ApplyChimeraProcessMonolithic<2>::FormulateChimera)
		.def("SetOverlapDistance",&ApplyChimeraProcessMonolithic<2>::SetOverlapDistance)
		;

	class_<ApplyChimeraProcessMonolithic<3>, ApplyChimeraProcessMonolithic<3>::Pointer, Process>(m, "ApplyChimeraProcessMonolithic3d")
		.def(init< ModelPart&, Parameters >())
		.def("ApplyMpcConstraint", &ApplyChimeraProcessMonolithic<3>::ApplyMpcConstraint)
		.def("FormulateChimera2D", &ApplyChimeraProcessMonolithic<3>::FormulateChimera)
		.def("SetOverlapDistance",&ApplyChimeraProcessMonolithic<3>::SetOverlapDistance)
		;

	class_<ApplyChimeraProcessFractionalStep<2>, ApplyChimeraProcessFractionalStep<2>::Pointer, Process>(m, "ApplyChimeraProcessFractionalStep2d")
		.def(init< ModelPart&, Parameters >())
		.def("FormulateChimera2D", &ApplyChimeraProcessFractionalStep<2>::FormulateChimera)
		.def("SetOverlapDistance",&ApplyChimeraProcessFractionalStep<2>::SetOverlapDistance)
		.def("CalculateNodalAreaAndNodalMass",&ApplyChimeraProcessFractionalStep<2>::CalculateNodalAreaAndNodalMass)
		;

	class_<ApplyChimeraProcessFractionalStep<3>, ApplyChimeraProcessFractionalStep<3>::Pointer, Process>(m, "ApplyChimeraProcessFractionalStep3d")
		.def(init< ModelPart&, Parameters >())
		.def("FormulateChimera2D", &ApplyChimeraProcessFractionalStep<3>::FormulateChimera)
		.def("SetOverlapDistance",&ApplyChimeraProcessFractionalStep<3>::SetOverlapDistance)
		.def("CalculateNodalAreaAndNodalMass",&ApplyChimeraProcessFractionalStep<3>::CalculateNodalAreaAndNodalMass)
		;


    class_< RotateRegionProcess,RotateRegionProcess::Pointer, Process >(m, "RotateRegionProcess")
		.def(init< ModelPart&, Parameters >())
		.def("SetCentreOfRotation", &RotateRegionProcess::SetCentreOfRotation)
		.def("ChangeAngularVelocity", &RotateRegionProcess::ChangeAngularVelocity)
		;

	class_<CustomCalculateSignedDistanceProcess<2>>(m, "SignedDistanceProcess2d")
		.def(init<>())
		.def("CalculateSignedDistance", &CustomCalculateSignedDistanceProcess<2>::CalculateSignedDistance)
		;

	class_<CustomCalculateSignedDistanceProcess<3>>(m, "SignedDistanceProcess3d")
		.def(init<>())
		.def("CalculateSignedDistance", &CustomCalculateSignedDistanceProcess<3>::CalculateSignedDistance)
		;

}



} // namespace Python.

} // Namespace Kratos

