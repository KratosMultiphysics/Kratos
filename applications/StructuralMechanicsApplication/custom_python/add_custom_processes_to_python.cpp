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

namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// Processes
    py::class_<PostprocessEigenvaluesProcess, PostprocessEigenvaluesProcess::Pointer, Process>(m,"PostprocessEigenvaluesProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<TotalStructuralMassProcess, TotalStructuralMassProcess::Pointer, Process>(m,"TotalStructuralMassProcess")
        .def(py::init<ModelPart&>())
        .def_static("CalculateElementMass", &TotalStructuralMassProcess::CalculateElementMass);
        ;

    py::class_<SolidShellThickComputeProcess, SolidShellThickComputeProcess::Pointer, Process>(m,"SolidShellThickComputeProcess")
        .def(py::init<ModelPart&>())
        ;


/*     class_<CableNetMpcProcess, CableNetMpcProcess::Pointer, ApplyMultipointConstraintsProcess>(m,"CableNetMpcProcess")
        .def(init<ModelPart&,Parameters&>())
        .def("CoupleModelParts", &CableNetMpcProcess::CoupleModelParts)
        ;  */
    py::class_<ComputeCenterOfGravityProcess, ComputeCenterOfGravityProcess::Pointer, Process>(m,"ComputeCenterOfGravityProcess")
        .def(py::init<ModelPart&>())
        ;

    py::class_<PrismNeighboursProcess, PrismNeighboursProcess::Pointer, Process>(m, "PrismNeighboursProcess")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, const bool >())

        ;

    py::class_<ShellToSolidShellProcess<3>, ShellToSolidShellProcess<3>::Pointer, Process>(m, "TriangleShellToSolidShellProcess")
        .def(py::init<ModelPart&>())
        .def(py::init< ModelPart&, Parameters >())
        ;

    py::class_<ShellToSolidShellProcess<4>, ShellToSolidShellProcess<4>::Pointer, Process>(m, "QuadrilateralShellToSolidShellProcess")
        .def(py::init<ModelPart&>())
        .def(py::init< ModelPart&, Parameters >())
        ;

    //SPR_ERROR
    py::class_<SPRErrorProcess<2>, SPRErrorProcess<2>::Pointer, Process >(m, "SPRErrorProcess2D")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, Parameters>())
        ;

    py::class_<SPRErrorProcess<3>, SPRErrorProcess<3>::Pointer, Process >(m, "SPRErrorProcess3D")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, Parameters>())
        ;

    py::class_<ImposeRigidMovementProcess, ImposeRigidMovementProcess::Pointer, Process>(m, "ImposeRigidMovementProcess")
        .def(py::init<ModelPart&>())
        .def(py::init< ModelPart&, Parameters >())
        ;
}

}  // namespace Python.
} // Namespace Kratos

