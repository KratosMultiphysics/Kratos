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
#include "custom_processes/apply_chimera_process_Monolithic.h"
#include "custom_processes/apply_chimera_process_FractionalStep.h"
#include "custom_processes/custom_calculate_signed_distance_process.h"
#include "custom_processes/calculate_signed_distance_to_2d_condition_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"
#include "custom_processes/rotate_region_process.h"
//#include "custom_processes/apply_multi_point_constraints_process.h"
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
		.def("CalculateNodalAreaAndNodalMass",&ApplyChimeraProcessMonolithic<2>::CalculateNodalAreaAndNodalMass)
		.def("SetType",&ApplyChimeraProcessMonolithic<2>::SetType)
		;

	class_<ApplyChimeraProcessMonolithic<3>, ApplyChimeraProcessMonolithic<3>::Pointer, Process>(m, "ApplyChimeraProcessMonolithic3d")
		.def(init< ModelPart&, Parameters >())
		.def("ApplyMpcConstraint", &ApplyChimeraProcessMonolithic<3>::ApplyMpcConstraint)
		.def("FormulateChimera2D", &ApplyChimeraProcessMonolithic<3>::FormulateChimera)
		.def("SetOverlapDistance",&ApplyChimeraProcessMonolithic<3>::SetOverlapDistance)
		.def("CalculateNodalAreaAndNodalMass",&ApplyChimeraProcessMonolithic<3>::CalculateNodalAreaAndNodalMass)
		.def("SetType",&ApplyChimeraProcessMonolithic<3>::SetType)
		;

	class_<ApplyChimeraProcessFractionalStep<2>, ApplyChimeraProcessFractionalStep<2>::Pointer, Process>(m, "ApplyChimeraProcessFractionalStep2d")
		.def(init< ModelPart&, Parameters >())
		.def("ApplyMpcConstraint", &ApplyChimeraProcessFractionalStep<2>::ApplyMpcConstraint)
		.def("FormulateChimera2D", &ApplyChimeraProcessFractionalStep<2>::FormulateChimera)
		.def("SetOverlapDistance",&ApplyChimeraProcessFractionalStep<2>::SetOverlapDistance)
		.def("CalculateNodalAreaAndNodalMass",&ApplyChimeraProcessFractionalStep<2>::CalculateNodalAreaAndNodalMass)
		.def("SetType",&ApplyChimeraProcessFractionalStep<2>::SetType)
		;

	class_<ApplyChimeraProcessFractionalStep<3>, ApplyChimeraProcessFractionalStep<3>::Pointer, Process>(m, "ApplyChimeraProcessFractionalStep3d")
		.def(init< ModelPart&, Parameters >())
		.def("ApplyMpcConstraint", &ApplyChimeraProcessFractionalStep<3>::ApplyMpcConstraint)
		.def("FormulateChimera2D", &ApplyChimeraProcessFractionalStep<3>::FormulateChimera)
		.def("SetOverlapDistance",&ApplyChimeraProcessFractionalStep<3>::SetOverlapDistance)
		.def("CalculateNodalAreaAndNodalMass",&ApplyChimeraProcessFractionalStep<3>::CalculateNodalAreaAndNodalMass)
		.def("SetType",&ApplyChimeraProcessFractionalStep<3>::SetType)
		;

	/* class_<ApplyChimeraProcess<3>,bases<Process> >("ApplyChimeraProcess3d", init< ModelPart&, Parameters >())
			.def("ApplyMpcConstraint", &ApplyChimeraProcess<3>::ApplyMpcConstraint)
			.def("FormulateChimera3D", &ApplyChimeraProcess<3>::FormulateChimera)
			.def("SetOverlapDistance",&ApplyChimeraProcess<3>::SetOverlapDistance)
			.def("CalculateNodalAreaAndNodalMass",&ApplyChimeraProcess<3>::CalculateNodalAreaAndNodalMass)

			.def("SetType",&ApplyChimeraProcess<3>::SetType);
 */

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

	/* typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


   class_< SpalartAllmarasTurbulenceModelForChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases<Process>, boost::noncopyable >
    ("SpalartAllmarasTurbulenceModelForChimera", init < ModelPart&, LinearSolverType::Pointer, std::size_t, double, std::size_t, bool, std::size_t>())
    .def("ActivateDES", &SpalartAllmarasTurbulenceModelForChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >::ActivateDES)
    .def("AdaptForFractionalStep", &SpalartAllmarasTurbulenceModelForChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >::AdaptForFractionalStep)
    ; */

}



} // namespace Python.

} // Namespace Kratos

