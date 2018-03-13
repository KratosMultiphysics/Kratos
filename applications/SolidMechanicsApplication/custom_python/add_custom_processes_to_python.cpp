//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes 
#include <boost/python.hpp>

// External includes 

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"
#include "containers/flags.h"

//Application includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/transfer_entities_between_model_parts_process.h"
#include "custom_processes/transfer_nodes_to_model_part_process.h"
#include "custom_processes/assign_scalar_field_to_nodes_process.h"
#include "custom_processes/assign_scalar_field_to_conditions_process.h"
#include "custom_processes/assign_scalar_variable_to_nodes_process.h"
#include "custom_processes/assign_scalar_variable_to_conditions_process.h"
#include "custom_processes/assign_vector_variable_to_conditions_process.h"
#include "custom_processes/assign_vector_field_to_conditions_process.h"
#include "custom_processes/fix_scalar_dof_process.h"
#include "custom_processes/free_scalar_dof_process.h"
#include "custom_processes/add_dofs_process.h"
#include "custom_processes/assign_rotation_field_about_an_axis_to_nodes_process.h"
#include "custom_processes/assign_torque_field_about_an_axis_to_conditions_process.h"

namespace Kratos
{
	
  namespace Python
  {

    typedef std::vector<Flags>  FlagsContainer;
    
    
    void Push_Back_Flags( FlagsContainer& ThisFlagContainer,
			  Flags ThisFlag )
    {
      ThisFlagContainer.push_back( ThisFlag );
    }

    
    void  AddCustomProcessesToPython()
    {

      using namespace boost::python;
      typedef Process                                         ProcessBaseType;


      class_< FlagsContainer >( "FlagsContainer", init<>() )
	.def( "PushBack", Push_Back_Flags )
	;
    

      //**********TRANSFER NODES TO MODEL PART*********//

      class_<TransferNodesToModelPartProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "TransferNodesProcess", init<ModelPart&, ModelPart&, const FlagsContainer&>()
      	)
	.def(init<ModelPart&, ModelPart&, const FlagsContainer&, const FlagsContainer& >())
        .def("Execute", &TransferNodesToModelPartProcess::Execute)
      	;
      
      //**********TRANSFER ENTITIES BETWEEN MODEL PARTS*********//

      class_<TransferEntitiesBetweenModelPartsProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "TransferEntitiesProcess", init<ModelPart&, ModelPart&, const std::string>()
      	)	
	.def(init<ModelPart&, ModelPart&, const std::string, const FlagsContainer&>())
	.def(init<ModelPart&, ModelPart&, const std::string, const FlagsContainer&, const FlagsContainer& >())
        .def("Execute", &TransferEntitiesBetweenModelPartsProcess::Execute)
      	;
      
      
      //**********ASSIGN VALUES TO VARIABLES PROCESSES*********//

      class_<AssignScalarVariableToNodesProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AssignScalarToNodesProcess", init<ModelPart&, Parameters>()
      	)
        .def(init< ModelPart&, Parameters& >())
        .def(init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&, double>())
        .def(init<ModelPart&, const Variable<double>&, double>())
        .def(init<ModelPart&, const Variable<int>&, int>())
        .def(init<ModelPart&, const Variable<bool>&, bool>())
        .def("Execute", &AssignScalarVariableToNodesProcess::Execute)

      	;

      class_<AssignScalarFieldToNodesProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AssignScalarFieldToNodesProcess", init<ModelPart&, PyObject* ,const char* ,const bool, Parameters>()
      	)
        .def(init< ModelPart&, PyObject* ,const char* ,const bool, Parameters& >())
        .def("Execute", &AssignScalarFieldToNodesProcess::Execute)

      	;

      class_<AssignScalarVariableToConditionsProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AssignScalarToConditionsProcess", init<ModelPart&, Parameters>()
      	)
        .def(init< ModelPart&, Parameters& >())
        .def(init<ModelPart&, const Variable<double>&, double>())
        .def(init<ModelPart&, const Variable<int>&, int>())
        .def(init<ModelPart&, const Variable<bool>&, bool>())
        .def("Execute", &AssignScalarVariableToConditionsProcess::Execute)
      	;

      class_<AssignVectorVariableToConditionsProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AssignVectorToConditionsProcess", init<ModelPart&, Parameters>()
      	)
        .def(init< ModelPart&, Parameters& >())
        .def(init<ModelPart&, const Variable<array_1d<double,3> >&, array_1d<double,3>&>())
        .def("Execute", &AssignVectorVariableToConditionsProcess::Execute)
      	;

      class_<AssignScalarFieldToConditionsProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AssignScalarFieldToConditionsProcess", init<ModelPart&, PyObject* ,const char* ,const bool, Parameters>()
      	)
        .def(init< ModelPart&, PyObject* ,const char* ,const bool, Parameters& >())
        .def("Execute", &AssignScalarFieldToConditionsProcess::Execute)

      	;

      class_<AssignVectorFieldToConditionsProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AssignVectorFieldToConditionsProcess", init<ModelPart&, PyObject* ,const char* ,const bool, Parameters>()
      	)
        .def(init< ModelPart&, PyObject* ,const char* ,const bool, Parameters& >())
        .def("Execute", &AssignVectorFieldToConditionsProcess::Execute)

      	;
	
      //**********FIX AND FREE DOFS PROCESSES*********//

      class_<FixScalarDofProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "FixScalarDofProcess", init<ModelPart&, Parameters>()
      	)
        .def(init<ModelPart&, Parameters&>())
        .def(init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&>())
        .def(init<ModelPart&, const Variable<double>&>())
        .def(init<ModelPart&, const Variable<int>&>())
        .def(init<ModelPart&, const Variable<bool>&>())
        .def("Execute", &FixScalarDofProcess::Execute)

      	;


      class_<FreeScalarDofProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "FreeScalarDofProcess", init<ModelPart&, Parameters>()
      	)
        .def(init<ModelPart&, Parameters&>())
        .def(init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&>())
        .def(init<ModelPart&, const Variable<double>&>())
        .def(init<ModelPart&, const Variable<int>&>())
        .def(init<ModelPart&, const Variable<bool>&>())
        .def("Execute", &FreeScalarDofProcess::Execute)

      	;

 
      //**********ADD DOFS PROCESS*********//

      class_<AddDofsProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AddDofsProcess", init<ModelPart&, Parameters>()
      	)
        .def(init<ModelPart&, Parameters&>())
	.def(init<ModelPart&, const boost::python::list&, const boost::python::list&>())
        .def("Execute", &AddDofsProcess::Execute)

      	;


      //**********ASSIGN ROTATION ABOUT AND AXIS*********//

      class_<AssignRotationAboutAnAxisToNodesProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AssignRotationAboutAnAxisToNodesProcess", init<ModelPart&, Parameters>()
      	)
        .def(init< ModelPart&, Parameters& >())
        .def("Execute", &AssignRotationAboutAnAxisToNodesProcess::Execute)

      	;

      
      class_<AssignRotationFieldAboutAnAxisToNodesProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AssignRotationFieldAboutAnAxisToNodesProcess", init<ModelPart&, PyObject* ,const char* ,const bool, Parameters>()
      	)
	.def(init< ModelPart&, PyObject* ,const char* ,const bool, Parameters& >())
        .def("Execute", &AssignRotationFieldAboutAnAxisToNodesProcess::Execute)

      	;

      //**********ASSIGN TORQUE ABOUT AN AXIS*********//

      class_<AssignTorqueAboutAnAxisToConditionsProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AssignTorqueAboutAnAxisToConditionsProcess", init<ModelPart&, Parameters>()
      	)
        .def(init< ModelPart&, Parameters& >())
        .def("Execute", &AssignTorqueAboutAnAxisToConditionsProcess::Execute)

      	;

      
      class_<AssignTorqueFieldAboutAnAxisToConditionsProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AssignTorqueFieldAboutAnAxisToConditionsProcess", init<ModelPart&, PyObject* ,const char* ,const bool, Parameters>()
      	)
	.def(init< ModelPart&, PyObject* ,const char* ,const bool, Parameters& >())
        .def("Execute", &AssignTorqueFieldAboutAnAxisToConditionsProcess::Execute)

      	;
      
      
    }
 
  }  // namespace Python.

} // Namespace Kratos

