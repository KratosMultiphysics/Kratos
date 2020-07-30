//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/elemental_refining_criteria_process.h"
#include "custom_processes/apply_perturbation_function_process.h"
#include "custom_processes/apply_sinusoidal_function_process.h"
#include "custom_processes/rough_porous_layer_wetting_model.h"
#include "custom_processes/negative_height_wetting_model.h"
#include "custom_processes/id_renumbering_process.h"
#include "custom_processes/compute_velocity_process.h"
#include "custom_processes/move_shallow_particles_process.h"


namespace Kratos
{

namespace Python
{

    void  AddCustomProcessesToPython(pybind11::module& m)
    {
        namespace py = pybind11;

        typedef VariableComponent<VectorComponentAdaptor<array_1d<double,3>>> VariableComponentType;

        py::class_<ElementalRefiningCriteriaProcess, ElementalRefiningCriteriaProcess::Pointer, Process>
        (m, "ElementalRefiningCriteriaProcess")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, Parameters>())
        .def(py::init<ModelPart&, const Variable<double>&, double, bool>())
        ;

        typedef ApplyPerturbationFunctionProcess<Variable<double>> ApplyPerturbationScalarFunctionProcess;
        py::class_<ApplyPerturbationScalarFunctionProcess, ApplyPerturbationScalarFunctionProcess::Pointer, Process>
        (m, "ApplyPerturbationFunctionToScalar")
        .def(py::init<ModelPart&, Node<3>::Pointer, Variable<double>&, Parameters&>())
        .def(py::init<ModelPart&, ModelPart::NodesContainerType&, Variable<double>&, Parameters&>())
        ;

        typedef ApplyPerturbationFunctionProcess<VariableComponentType> ApplyPerturbationComponentFunctionProcess;
        py::class_<ApplyPerturbationComponentFunctionProcess, ApplyPerturbationComponentFunctionProcess::Pointer, Process>
        (m, "ApplyPerturbationFunctionToComponent")
        .def(py::init<ModelPart&, Node<3>::Pointer, VariableComponentType&, Parameters&>())
        .def(py::init<ModelPart&, ModelPart::NodesContainerType&, VariableComponentType&, Parameters&>())
        ;

        typedef ApplySinusoidalFunctionProcess<Variable<double>> ApplySinusoidalScalarFunctionProcess;
        py::class_<ApplySinusoidalScalarFunctionProcess, ApplySinusoidalScalarFunctionProcess::Pointer, Process>
        (m, "ApplySinusoidalFunctionToScalar")
        .def(py::init<ModelPart&, Variable<double>&, Parameters&>())
        ;

        typedef ApplySinusoidalFunctionProcess<VariableComponentType> ApplySinusoidalComponentFunctionProcess;
        py::class_<ApplySinusoidalComponentFunctionProcess, ApplySinusoidalComponentFunctionProcess::Pointer, Process>
        (m, "ApplySinusoidalFunctionToComponent")
        .def(py::init<ModelPart&, VariableComponentType&, Parameters&>())
        ;

        typedef ApplySinusoidalFunctionProcess<Variable<array_1d<double,3>>> ApplySinusoidalVectorFunctionProcess;
        py::class_<ApplySinusoidalVectorFunctionProcess, ApplySinusoidalVectorFunctionProcess::Pointer, Process>
        (m, "ApplySinusoidalFunctionToVector")
        .def(py::init<ModelPart&, Variable<array_1d<double,3>>&, Parameters&>())
        ;

        py::class_<RoughPorousLayerWettingModel, RoughPorousLayerWettingModel::Pointer, Process>
        (m, "RoughPorousLayerWettingModel")
        .def(py::init<ModelPart&, Parameters>())
        .def(py::init<ModelPart&, double, double>())
        ;

        py::class_<NegativeHeightWettingModel, NegativeHeightWettingModel::Pointer, Process>
        (m, "NegativeHeightWettingModel")
        .def(py::init<ModelPart&, Parameters>())
        .def(py::init<ModelPart&, double>())
        ;

        py::class_<IdRenumberingProcess, IdRenumberingProcess::Pointer, Process>
        (m, "IdRenumberingProcess")
        .def(py::init<Model&>())
        .def(py::init<Model&, StringVectorType&>())
        .def("RenumberNodes", &IdRenumberingProcess::RenumberNodes)
        .def("RenumberElements", &IdRenumberingProcess::RenumberElements)
        .def("RenumberConditions", &IdRenumberingProcess::RenumberConditions)
        .def("RestoreNodes", &IdRenumberingProcess::RestoreNodes)
        .def("RestoreElements", &IdRenumberingProcess::RestoreElements)
        .def("RestoreConditions", &IdRenumberingProcess::RestoreConditions)
        ;

        py::class_<ComputeVelocityProcess, ComputeVelocityProcess::Pointer, Process>
        (m, "ComputeVelocityProcess")
        .def(py::init<ModelPart&, double>())
        ;

        py::class_<MoveShallowParticlesProcess<2>, MoveShallowParticlesProcess<2>::Pointer, Process>
        (m, "MoveShallowParticlesProcess2D")
        .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3>>&, Variable<double>&, Parameters>())
        ;

    }

}  // namespace Python.

} // Namespace Kratos
