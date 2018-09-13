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
#include "custom_processes/transfer_computing_model_part_entities_process.hpp"

// Solver Processes
#include "custom_processes/solver_process.hpp"

namespace Kratos
{

namespace Python
{

typedef std::vector<Flags>  FlagsContainer;

void  AddCustomProcessesToPython(pybind11::module& m)
{

  using namespace pybind11;

  //**********ASSIGN FLAGS TO MODEL PART ENTITIES*********//

  class_<AssignFlagsToModelPartEntitiesProcess, AssignFlagsToModelPartEntitiesProcess::Pointer, Process>(m,"AssignFlagsToEntitiesProcess")
      .def(init<ModelPart&, const std::string, const FlagsContainer&>())
      .def(init<ModelPart&, const std::string, const FlagsContainer&, const FlagsContainer& >())
      .def("Execute", &AssignFlagsToModelPartEntitiesProcess::Execute)
      ;

  //**********TRANSFER ENTITIES BETWEEN MODEL PARTS*********//

  class_<TransferEntitiesBetweenModelPartsProcess, TransferEntitiesBetweenModelPartsProcess::Pointer, Process>(m,"TransferEntitiesProcess")
      .def(init<ModelPart&, ModelPart&, const std::string>())
      .def(init<ModelPart&, ModelPart&, const std::string, const FlagsContainer&>())
      .def(init<ModelPart&, ModelPart&, const std::string, const FlagsContainer&, const FlagsContainer& >())
      .def("Execute", &TransferEntitiesBetweenModelPartsProcess::Execute)
      ;

  class_<TransferComputingModelPartEntitiesProcess, TransferComputingModelPartEntitiesProcess::Pointer, Process>(m,"TransferComputingModelPartProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init<ModelPart&, Parameters& >())
      ;


  //**********ASSIGN VALUES TO VARIABLES PROCESSES*********//

  class_<AssignScalarVariableToEntitiesProcess, AssignScalarVariableToEntitiesProcess::Pointer, Process>(m,"AssignScalarToEntitiesProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init< ModelPart&, Parameters& >())
      .def("Execute", &AssignScalarVariableToEntitiesProcess::Execute)
      ;

  class_<AssignScalarFieldToEntitiesProcess, AssignScalarFieldToEntitiesProcess::Pointer, AssignScalarVariableToEntitiesProcess>(m,"AssignScalarFieldToEntitiesProcess")
      .def(init<ModelPart&, pybind11::object&, const std::string, const bool, Parameters>())
      .def(init< ModelPart&, pybind11::object&, const std::string, const bool, Parameters& >())
      .def("Execute", &AssignScalarFieldToEntitiesProcess::Execute)
      ;

  class_<AssignVectorFieldToEntitiesProcess, AssignVectorFieldToEntitiesProcess::Pointer, AssignScalarFieldToEntitiesProcess>(m,"AssignVectorFieldToEntitiesProcess")
      .def(init<ModelPart&, pybind11::object&,const std::string,const bool, Parameters>())
      .def(init< ModelPart&, pybind11::object&,const std::string,const bool, Parameters& >())
      .def("Execute", &AssignVectorFieldToEntitiesProcess::Execute)
      ;

  class_<AssignVectorVariableToConditionsProcess, AssignVectorVariableToConditionsProcess::Pointer, AssignScalarVariableToEntitiesProcess>(m,"AssignVectorToConditionsProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init< ModelPart&, Parameters& >())
      .def(init<ModelPart&, const Variable<array_1d<double,3> >&, array_1d<double,3>&>())
      .def("Execute", &AssignVectorVariableToConditionsProcess::Execute)
      ;

  //**********FIX AND FREE DOFS PROCESSES*********//

  class_<FixScalarDofProcess, FixScalarDofProcess::Pointer, Process>(m,"FixScalarDofProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init<ModelPart&, Parameters&>())
      .def(init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&>())
      .def(init<ModelPart&, const Variable<double>&>())
      .def(init<ModelPart&, const Variable<int>&>())
      .def(init<ModelPart&, const Variable<bool>&>())
      .def("Execute", &FixScalarDofProcess::Execute)

      ;


  class_<FreeScalarDofProcess, FreeScalarDofProcess::Pointer, Process>(m,"FreeScalarDofProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init<ModelPart&, Parameters&>())
      .def(init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&>())
      .def(init<ModelPart&, const Variable<double>&>())
      .def(init<ModelPart&, const Variable<int>&>())
      .def(init<ModelPart&, const Variable<bool>&>())
      .def("Execute", &FreeScalarDofProcess::Execute)

      ;


  //**********ADD DOFS PROCESS*********//

  class_<AddDofsProcess, AddDofsProcess::Pointer, Process>(m,"AddDofsProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init<ModelPart&, Parameters&>())
      .def(init<ModelPart&, const pybind11::list&, const pybind11::list&>())
      .def("Execute", &AddDofsProcess::Execute)

      ;


  //**********ASSIGN ROTATION ABOUT AND AXIS*********//

  class_<AssignRotationAboutAnAxisToNodesProcess, AssignRotationAboutAnAxisToNodesProcess::Pointer, Process>(m,"AssignRotationAboutAnAxisToNodesProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init< ModelPart&, Parameters& >())
      .def("Execute", &AssignRotationAboutAnAxisToNodesProcess::Execute)

      ;


  class_<AssignRotationFieldAboutAnAxisToNodesProcess, AssignRotationFieldAboutAnAxisToNodesProcess::Pointer, Process>(m,"AssignRotationFieldAboutAnAxisToNodesProcess")
      .def(init<ModelPart&, pybind11::object&, const std::string,const bool, Parameters>())
      .def(init< ModelPart&, pybind11::object&, const std::string,const bool, Parameters& >())
      .def("Execute", &AssignRotationFieldAboutAnAxisToNodesProcess::Execute)

      ;

  //**********ASSIGN TORQUE ABOUT AN AXIS*********//

  class_<AssignTorqueAboutAnAxisToConditionsProcess, AssignTorqueAboutAnAxisToConditionsProcess::Pointer, Process>(m,"AssignTorqueAboutAnAxisToConditionsProcess")
      .def(init< ModelPart&, Parameters >())
      .def(init< ModelPart&, Parameters& >())
      .def("Execute", &AssignTorqueAboutAnAxisToConditionsProcess::Execute)

      ;


  class_<AssignTorqueFieldAboutAnAxisToConditionsProcess, AssignTorqueFieldAboutAnAxisToConditionsProcess::Pointer, Process>(m,"AssignTorqueFieldAboutAnAxisToConditionsProcess")
      .def(init< ModelPart&, pybind11::object&,const std::string,const bool, Parameters >())
      .def(init< ModelPart&, pybind11::object&,const std::string,const bool, Parameters& >())
      .def("Execute", &AssignTorqueFieldAboutAnAxisToConditionsProcess::Execute)

      ;


  //**********BUILD STRING SKIN PROCESS*********//

  class_<BuildStringSkinProcess, BuildStringSkinProcess::Pointer, Process>(m,"BuildStringSkinProcess")
      .def(init<ModelPart&, unsigned int, double>())
      .def("ExecuteInitialize", &BuildStringSkinProcess::ExecuteInitialize)
      .def("ExecuteFinalizeSolutionStep", &BuildStringSkinProcess::ExecuteFinalizeSolutionStep)
      .def("ExecuteBeforeOutputStep", &BuildStringSkinProcess::ExecuteBeforeOutputStep)
      .def("ExecuteAfterOutputStep", &BuildStringSkinProcess::ExecuteAfterOutputStep)
      ;



  //**********SOLVER PROCESS*********//

  class_<SolverProcess, SolverProcess::Pointer, Process>(m,"SolverProcess")
      .def(init<>())
      ;



}

}  // namespace Python.

} // Namespace Kratos
