// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "structural_mechanics_application_variables.h"

//Processes
#include "custom_processes/prism_neighbours_process.h"
#include "custom_processes/postprocess_eigenvalues_process.h"
#include "custom_processes/total_structural_mass_process.h"
#include "custom_processes/compute_center_of_gravity_process.h"
#include "custom_processes/compute_mass_moment_of_inertia_process.h"
#include "custom_processes/shell_to_solid_shell_process.h"
#include "custom_processes/solid_shell_thickness_compute_process.h"
#include "custom_processes/spr_error_process.h"
#include "custom_processes/assign_nodal_elements_to_nodes_process.h"
#include "custom_processes/impose_rigid_movement_process.h"
#include "custom_processes/impose_z_strain_process.h"
#include "custom_processes/distribute_load_on_surface_process.h"
#include "custom_processes/set_moving_load_process.h"
#include "custom_processes/set_cartesian_local_axes_process.h"
#include "custom_processes/set_cylindrical_local_axes_process.h"
#include "custom_processes/set_spherical_local_axes_process.h"
#include "custom_processes/set_automated_initial_variable_process.h"

namespace Kratos::Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// Processes
    py::class_<PostprocessEigenvaluesProcess, PostprocessEigenvaluesProcess::Pointer, Process>(m,"PostprocessEigenvaluesProcess")
        .def(py::init<Model&, Parameters>())
        ;

    py::class_<TotalStructuralMassProcess, TotalStructuralMassProcess::Pointer, Process>(m,"TotalStructuralMassProcess")
        .def(py::init<ModelPart&>())
        .def_static("CalculateElementMass", &TotalStructuralMassProcess::CalculateElementMass);
        ;

    py::class_<ComputeCenterOfGravityProcess, ComputeCenterOfGravityProcess::Pointer, Process>(m,"ComputeCenterOfGravityProcess")
        .def(py::init<ModelPart&>())
        ;

    py::class_<ComputeMassMomentOfInertiaProcess, ComputeMassMomentOfInertiaProcess::Pointer, Process>(m,"ComputeMassMomentOfInertiaProcess")
        .def(py::init<ModelPart&, const Point&, const Point&>())
        ;

    py::class_<SolidShellThickComputeProcess, SolidShellThickComputeProcess::Pointer, Process>(m,"SolidShellThickComputeProcess")
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

    py::class_<AssignNodalElementsToNodesProcess, AssignNodalElementsToNodesProcess::Pointer, Process>(m, "AssignNodalElementsToNodesProcess")
        .def(py::init<Model&, Parameters>())
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, Parameters >())
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

    py::class_<ImposeZStrainProcess, ImposeZStrainProcess::Pointer, Process>(m, "ImposeZStrainProcess")
        .def(py::init< ModelPart&, Parameters >())
        ;

    py::class_<DistributeLoadOnSurfaceProcess, DistributeLoadOnSurfaceProcess::Pointer, Process>(m,"DistributeLoadOnSurfaceProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<SetMovingLoadProcess, SetMovingLoadProcess::Pointer, Process>(m, "SetMovingLoadProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<SetCartesianLocalAxesProcess, SetCartesianLocalAxesProcess::Pointer, Process>(m,"SetCartesianLocalAxesProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<SetCylindricalLocalAxesProcess, SetCylindricalLocalAxesProcess::Pointer, Process>(m,"SetCylindricalLocalAxesProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<SetSphericalLocalAxesProcess, SetSphericalLocalAxesProcess::Pointer, Process>(m,"SetSphericalLocalAxesProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<SetAutomatedInitialVariableProcess, SetAutomatedInitialVariableProcess::Pointer, Process>(m,"SetAutomatedInitialVariableProcess")
        .def(py::init<ModelPart&, Parameters>());
}

}  // namespace Kratos::Python

