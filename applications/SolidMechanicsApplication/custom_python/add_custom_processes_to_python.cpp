//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <pybind11/stl.h>

// External includes

// Project includes
#include "custom_python/add_custom_processes_to_python.h"

// Processes
#include "custom_processes/transfer_solving_model_part_entities_process.hpp"
#include "custom_processes/transfer_entities_between_model_parts_process.hpp"
#include "custom_processes/assign_flags_to_model_part_entities_process.hpp"
#include "custom_processes/assign_scalar_variable_to_entities_process.hpp"
#include "custom_processes/assign_vector_variable_to_conditions_process.hpp"
#include "custom_processes/assign_vector_field_to_entities_process.hpp"
#include "custom_processes/fix_scalar_dof_process.hpp"
#include "custom_processes/free_scalar_dof_process.hpp"
#include "custom_processes/add_dofs_process.hpp"
#include "custom_processes/assign_rotation_field_about_an_axis_to_nodes_process.hpp"
#include "custom_processes/assign_torque_field_about_an_axis_to_conditions_process.hpp"
#include "custom_processes/build_string_skin_process.hpp"


// Solver Processes
#include "custom_processes/solver_process.hpp"

namespace Kratos
{

namespace Python
{

typedef std::vector<Flags>  FlagsContainer;

void  AddCustomProcessesToPython(pybind11::module& m)
{

  namespace py = pybind11;

  //**********ASSIGN FLAGS TO MODEL PART ENTITIES*********//

  py::class_<AssignFlagsToModelPartEntitiesProcess, AssignFlagsToModelPartEntitiesProcess::Pointer, Process>(m,"AssignFlagsToEntitiesProcess")
      .def(py::init<ModelPart&, const std::string, const FlagsContainer&>())
      .def(py::init<ModelPart&, const std::string, const FlagsContainer&, const FlagsContainer& >())
      .def("Execute", &AssignFlagsToModelPartEntitiesProcess::Execute)
      ;

  //**********TRANSFER ENTITIES BETWEEN MODEL PARTS*********//

  py::class_<TransferEntitiesBetweenModelPartsProcess, TransferEntitiesBetweenModelPartsProcess::Pointer, Process>(m,"TransferEntitiesProcess")
      .def(py::init<ModelPart&, ModelPart&, const std::string>())
      .def(py::init<ModelPart&, ModelPart&, const std::string, const FlagsContainer&>())
      .def(py::init<ModelPart&, ModelPart&, const std::string, const FlagsContainer&, const FlagsContainer& >())
      .def("Execute", &TransferEntitiesBetweenModelPartsProcess::Execute)
      ;

  py::class_<TransferSolvingModelPartEntitiesProcess, TransferSolvingModelPartEntitiesProcess::Pointer, Process>(m,"TransferSolvingModelPartProcess")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init<ModelPart&, Parameters& >())
      ;


  //**********ASSIGN VALUES TO VARIABLES PROCESSES*********//

  py::class_<AssignScalarVariableToEntitiesProcess, AssignScalarVariableToEntitiesProcess::Pointer, Process>(m,"AssignScalarToEntitiesProcess")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init< ModelPart&, Parameters& >())
      .def("Execute", &AssignScalarVariableToEntitiesProcess::Execute)
      ;

  py::class_<AssignScalarFieldToEntitiesProcess, AssignScalarFieldToEntitiesProcess::Pointer, AssignScalarVariableToEntitiesProcess>(m,"AssignScalarFieldToEntitiesProcess")
      .def(py::init<ModelPart&, pybind11::object&, const std::string, const bool, Parameters>())
      .def(py::init< ModelPart&, pybind11::object&, const std::string, const bool, Parameters& >())
      .def("Execute", &AssignScalarFieldToEntitiesProcess::Execute)
      ;

  py::class_<AssignVectorFieldToEntitiesProcess, AssignVectorFieldToEntitiesProcess::Pointer, AssignScalarFieldToEntitiesProcess>(m,"AssignVectorFieldToEntitiesProcess")
      .def(py::init<ModelPart&, pybind11::object&,const std::string,const bool, Parameters>())
      .def(py::init< ModelPart&, pybind11::object&,const std::string,const bool, Parameters& >())
      .def("Execute", &AssignVectorFieldToEntitiesProcess::Execute)
      ;

  py::class_<AssignVectorVariableToConditionsProcess, AssignVectorVariableToConditionsProcess::Pointer, AssignScalarVariableToEntitiesProcess>(m,"AssignVectorToConditionsProcess")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init< ModelPart&, Parameters& >())
      .def(py::init<ModelPart&, const Variable<array_1d<double,3> >&, array_1d<double,3>&>())
      .def("Execute", &AssignVectorVariableToConditionsProcess::Execute)
      ;

  //**********FIX AND FREE DOFS PROCESSES*********//

  py::class_<FixScalarDofProcess, FixScalarDofProcess::Pointer, Process>(m,"FixScalarDofProcess")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init<ModelPart&, Parameters&>())
      .def(py::init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&>())
      .def(py::init<ModelPart&, const Variable<double>&>())
      .def(py::init<ModelPart&, const Variable<int>&>())
      .def(py::init<ModelPart&, const Variable<bool>&>())
      .def("Execute", &FixScalarDofProcess::Execute)

      ;


  py::class_<FreeScalarDofProcess, FreeScalarDofProcess::Pointer, Process>(m,"FreeScalarDofProcess")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init<ModelPart&, Parameters&>())
      .def(py::init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&>())
      .def(py::init<ModelPart&, const Variable<double>&>())
      .def(py::init<ModelPart&, const Variable<int>&>())
      .def(py::init<ModelPart&, const Variable<bool>&>())
      .def("Execute", &FreeScalarDofProcess::Execute)

      ;


  //**********ADD DOFS PROCESS*********//

  py::class_<AddDofsProcess, AddDofsProcess::Pointer, Process>(m,"AddDofsProcess")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init<ModelPart&, Parameters&>())
      .def(py::init<ModelPart&, const pybind11::list&, const pybind11::list&>())
      .def("Execute", &AddDofsProcess::Execute)

      ;


  //**********ASSIGN ROTATION ABOUT AND AXIS*********//

  py::class_<AssignRotationAboutAnAxisToNodesProcess, AssignRotationAboutAnAxisToNodesProcess::Pointer, Process>(m,"AssignRotationAboutAnAxisToNodesProcess")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init< ModelPart&, Parameters& >())
      .def("Execute", &AssignRotationAboutAnAxisToNodesProcess::Execute)

      ;


  py::class_<AssignRotationFieldAboutAnAxisToNodesProcess, AssignRotationFieldAboutAnAxisToNodesProcess::Pointer, Process>(m,"AssignRotationFieldAboutAnAxisToNodesProcess")
      .def(py::init<ModelPart&, pybind11::object&, const std::string,const bool, Parameters>())
      .def(py::init< ModelPart&, pybind11::object&, const std::string,const bool, Parameters& >())
      .def("Execute", &AssignRotationFieldAboutAnAxisToNodesProcess::Execute)

      ;

  //**********ASSIGN TORQUE ABOUT AN AXIS*********//

  py::class_<AssignTorqueAboutAnAxisToConditionsProcess, AssignTorqueAboutAnAxisToConditionsProcess::Pointer, Process>(m,"AssignTorqueAboutAnAxisToConditionsProcess")
      .def(py::init< ModelPart&, Parameters >())
      .def(py::init< ModelPart&, Parameters& >())
      .def("Execute", &AssignTorqueAboutAnAxisToConditionsProcess::Execute)

      ;


  py::class_<AssignTorqueFieldAboutAnAxisToConditionsProcess, AssignTorqueFieldAboutAnAxisToConditionsProcess::Pointer, Process>(m,"AssignTorqueFieldAboutAnAxisToConditionsProcess")
      .def(py::init< ModelPart&, pybind11::object&,const std::string,const bool, Parameters >())
      .def(py::init< ModelPart&, pybind11::object&,const std::string,const bool, Parameters& >())
      .def("Execute", &AssignTorqueFieldAboutAnAxisToConditionsProcess::Execute)

      ;


  //**********BUILD STRING SKIN PROCESS*********//

  py::class_<BuildStringSkinProcess, BuildStringSkinProcess::Pointer, Process>(m,"BuildStringSkinProcess")
      .def(py::init<ModelPart&, unsigned int, double>())
      .def("ExecuteInitialize", &BuildStringSkinProcess::ExecuteInitialize)
      .def("ExecuteFinalizeSolutionStep", &BuildStringSkinProcess::ExecuteFinalizeSolutionStep)
      .def("ExecuteBeforeOutputStep", &BuildStringSkinProcess::ExecuteBeforeOutputStep)
      .def("ExecuteAfterOutputStep", &BuildStringSkinProcess::ExecuteAfterOutputStep)
      ;



  //**********SOLVER PROCESS*********//

  py::class_<SolverProcess, SolverProcess::Pointer, Process>(m,"SolverProcess")
      .def(py::init<>())
      ;



}

}  // namespace Python.

} // Namespace Kratos
