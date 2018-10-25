// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes


// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "structural_mechanics_application_variables.h"

//Processes
#include "custom_processes/prism_neighbours_process.h"
#include "custom_processes/postprocess_eigenvalues_process.h"
#include "custom_processes/total_structural_mass_process.h"
#include "custom_processes/compute_center_of_gravity_process.h"
#include "custom_processes/shell_to_solid_shell_process.h"
#include "custom_processes/solid_shell_thickness_compute_process.h"
#include "custom_processes/spr_error_process.h"
#include "custom_processes/impose_rigid_movement_process.h"

namespace Kratos
{
namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

    /// Processes
    class_<PostprocessEigenvaluesProcess, PostprocessEigenvaluesProcess::Pointer, Process>(m,"PostprocessEigenvaluesProcess")
        .def(init<ModelPart&, Parameters>());

    class_<TotalStructuralMassProcess, TotalStructuralMassProcess::Pointer, Process>(m,"TotalStructuralMassProcess")
        .def(init<ModelPart&>())
        .def_static("CalculateElementMass", &TotalStructuralMassProcess::CalculateElementMass);
        ;

    class_<ComputeCenterOfGravityProcess, ComputeCenterOfGravityProcess::Pointer, Process>(m,"ComputeCenterOfGravityProcess")
        .def(init<ModelPart&>())
        ;


    class_<SolidShellThickComputeProcess, SolidShellThickComputeProcess::Pointer, Process>(m,"SolidShellThickComputeProcess")
        .def(init<ModelPart&>())
        ;

    class_<PrismNeighboursProcess, PrismNeighboursProcess::Pointer, Process>(m, "PrismNeighboursProcess")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, const bool >())
        ;

    class_<ShellToSolidShellProcess<3>, ShellToSolidShellProcess<3>::Pointer, Process>(m, "TriangleShellToSolidShellProcess")
        .def(init<ModelPart&>())
        .def(init< ModelPart&, Parameters >())
        ;

    class_<ShellToSolidShellProcess<4>, ShellToSolidShellProcess<4>::Pointer, Process>(m, "QuadrilateralShellToSolidShellProcess")
        .def(init<ModelPart&>())
        .def(init< ModelPart&, Parameters >())
        ;

    //SPR_ERROR
    class_<SPRErrorProcess<2>, SPRErrorProcess<2>::Pointer, Process >(m, "SPRErrorProcess2D")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        ;

    class_<SPRErrorProcess<3>, SPRErrorProcess<3>::Pointer, Process >(m, "SPRErrorProcess3D")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        ;

    class_<ImposeRigidMovementProcess, ImposeRigidMovementProcess::Pointer, Process>(m, "ImposeRigidMovementProcess")
        .def(init<ModelPart&>())
        .def(init< ModelPart&, Parameters >())
        ;
}

}  // namespace Python.

} // Namespace Kratos

