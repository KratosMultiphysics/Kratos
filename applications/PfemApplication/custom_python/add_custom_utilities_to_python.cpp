//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes 

// External includes 

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

//Utilities
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "custom_utilities/modeler_utilities.hpp"

namespace Kratos
{
	
namespace Python
{

// refining methods
    
void SetRefiningParameters(ModelerUtilities::RefiningParameters& rRefiningParameters, double AlphaParameter, double CriticalRadius, double CriticalSide)
{
  rRefiningParameters.SetParameters(AlphaParameter, CriticalRadius, CriticalSide);
}
    
    
void SetThresholdVariable(ModelerUtilities::RefiningParameters& rRefiningParameters, const Variable<double>& rVariable)
{
  rRefiningParameters.SetThresholdVariable(rVariable);
}

void SetErrorVariable(ModelerUtilities::RefiningParameters& rRefiningParameters, const Variable<double>& rVariable)
{
  rRefiningParameters.SetErrorVariable(rVariable);
}

Variable<double> GetThresholdVariable(ModelerUtilities::RefiningParameters& rRefiningParameters)
{
  return rRefiningParameters.GetThresholdVariable();

}
// remeshing methods

void SetReferenceElement(ModelerUtilities::MeshingParameters& rMeshingParameters, char* ElementName)
{
  rMeshingParameters.SetReferenceElement(KratosComponents<Element>::Get(ElementName));
}

void SetReferenceCondition(ModelerUtilities::MeshingParameters& rMeshingParameters, char* ConditionName)
{
  rMeshingParameters.SetReferenceCondition(KratosComponents<Condition>::Get(ConditionName));
}


// transfer methods
void TransferNodesToElementsOnThreshold( MeshDataTransferUtilities& rMeshDataTransfer, MeshDataTransferUtilities::TransferParameters& rTransferParameters, ModelerUtilities::RefiningParameters& rRefiningParameters, ModelPart& rModelPart)
{
  rMeshDataTransfer.TransferNodalValuesToElements( rTransferParameters, rRefiningParameters.GetThresholdVariable(), rRefiningParameters.GetReferenceThreshold(), rModelPart );
}

void SetDoubleVariable( MeshDataTransferUtilities::TransferParameters& rTransferParameters, const Variable<double>& rVariable)
{
  rTransferParameters.SetVariable( rVariable );
}
    
void SetArray1DVariable( MeshDataTransferUtilities::TransferParameters& rTransferParameters, const Variable<array_1d<double,3> >& rVariable )
{
  rTransferParameters.SetVariable( rVariable );
}
    
void SetVectorVariable( MeshDataTransferUtilities::TransferParameters& rTransferParameters, const Variable<Vector>& rVariable )
{
  rTransferParameters.SetVariable( rVariable );
}
    
void SetMatrixVariable( MeshDataTransferUtilities::TransferParameters& rTransferParameters, const Variable<Matrix>& rVariable )
{
  rTransferParameters.SetVariable( rVariable );
}


void  AddCustomUtilitiesToPython(pybind11::module& m)
{

  using namespace  pybind11;

 
  //***************DOMAIN SET**************//
  class_<ModelerUtilities>(m,"ModelerUtilities")
      .def(init<>())
      .def("SetModelPartNameToConditions",&ModelerUtilities::SetModelPartNameToConditions)
      .def("SetModelPartNameToNodes",&ModelerUtilities::SetModelPartNameToNodes)
      .def("CheckCriticalRadius",&ModelerUtilities::CheckCriticalRadius)
      .def("ComputeModelPartVolume",&ModelerUtilities::ComputeModelPartVolume)
	
      .def_readonly_static("REMESH",&ModelerUtilities::REMESH)
      .def_readonly_static("REFINE",&ModelerUtilities::REFINE)
      .def_readonly_static("RECONNECT",&ModelerUtilities::RECONNECT)
      .def_readonly_static("TRANSFER",&ModelerUtilities::TRANSFER)
      .def_readonly_static("CONSTRAINED",&ModelerUtilities::CONSTRAINED)
      .def_readonly_static("CONTACT_SEARCH",&ModelerUtilities::CONTACT_SEARCH)
      .def_readonly_static("MESH_SMOOTHING",&ModelerUtilities::MESH_SMOOTHING)
      .def_readonly_static("VARIABLES_SMOOTHING",&ModelerUtilities::VARIABLES_SMOOTHING)
	
      .def_readonly_static("REMOVE_NODES",&ModelerUtilities::REMOVE_NODES)
      .def_readonly_static("REMOVE_NODES_ON_DISTANCE",&ModelerUtilities::REMOVE_NODES_ON_DISTANCE)
      .def_readonly_static("REMOVE_NODES_ON_ERROR",&ModelerUtilities::REMOVE_NODES_ON_ERROR)
      .def_readonly_static("REMOVE_NODES_ON_THRESHOLD",&ModelerUtilities::REMOVE_NODES_ON_THRESHOLD)
      .def_readonly_static("REMOVE_BOUNDARY_NODES",&ModelerUtilities::REMOVE_BOUNDARY_NODES)
      .def_readonly_static("REMOVE_BOUNDARY_NODES_ON_DISTANCE",&ModelerUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE)
      .def_readonly_static("REMOVE_BOUNDARY_NODES_ON_ERROR",&ModelerUtilities::REMOVE_BOUNDARY_NODES_ON_ERROR)
      .def_readonly_static("REMOVE_BOUNDARY_NODES_ON_THRESHOLD",&ModelerUtilities::REMOVE_BOUNDARY_NODES_ON_THRESHOLD)
	
      .def_readonly_static("REFINE_ADD_NODES",&ModelerUtilities::REFINE_ADD_NODES)
      .def_readonly_static("REFINE_INSERT_NODES",&ModelerUtilities::REFINE_INSERT_NODES)
      .def_readonly_static("REFINE_ELEMENTS",&ModelerUtilities::REFINE_ELEMENTS)
      .def_readonly_static("REFINE_ELEMENTS_ON_DISTANCE",&ModelerUtilities::REFINE_ELEMENTS_ON_DISTANCE)
      .def_readonly_static("REFINE_ELEMENTS_ON_ERROR",&ModelerUtilities::REFINE_ELEMENTS_ON_ERROR)
      .def_readonly_static("REFINE_ELEMENTS_ON_THRESHOLD",&ModelerUtilities::REFINE_ELEMENTS_ON_THRESHOLD)
      .def_readonly_static("REFINE_BOUNDARY",&ModelerUtilities::REFINE_BOUNDARY)
      .def_readonly_static("REFINE_BOUNDARY_ON_DISTANCE",&ModelerUtilities::REFINE_BOUNDARY_ON_DISTANCE)
      .def_readonly_static("REFINE_BOUNDARY_ON_ERROR",&ModelerUtilities::REFINE_BOUNDARY_ON_ERROR)
      .def_readonly_static("REFINE_BOUNDARY_ON_THRESHOLD",&ModelerUtilities::REFINE_BOUNDARY_ON_THRESHOLD)

      .def_readonly_static("INITIALIZE_MESHER_INPUT",&ModelerUtilities::INITIALIZE_MESHER_INPUT)
      .def_readonly_static("FINALIZE_MESHER_INPUT",&ModelerUtilities::FINALIZE_MESHER_INPUT)
      .def_readonly_static("TRANSFER_KRATOS_NODES_TO_MESHER",&ModelerUtilities::TRANSFER_KRATOS_NODES_TO_MESHER)
      .def_readonly_static("TRANSFER_KRATOS_ELEMENTS_TO_MESHER",&ModelerUtilities::TRANSFER_KRATOS_ELEMENTS_TO_MESHER)
      .def_readonly_static("TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER",&ModelerUtilities::TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER)
      .def_readonly_static("TRANSFER_KRATOS_FACES_TO_MESHER",&ModelerUtilities::TRANSFER_KRATOS_FACES_TO_MESHER)

      .def_readonly_static("SELECT_TESSELLATION_ELEMENTS",&ModelerUtilities::SELECT_TESSELLATION_ELEMENTS)
      .def_readonly_static("KEEP_ISOLATED_NODES",&ModelerUtilities::KEEP_ISOLATED_NODES)
      ;
        
        
  //***************NORMALS**************//

  class_<BoundaryNormalsCalculationUtilities>(m,"BoundaryNormalsCalculation")
      .def(init<>())
      .def("CalculateWeightedBoundaryNormals", &BoundaryNormalsCalculationUtilities::CalculateWeightedBoundaryNormals)
      .def("CalculateUnitBoundaryNormals", &BoundaryNormalsCalculationUtilities::CalculateUnitBoundaryNormals)
      ;
      
  //***************TRANSFER UTILITIES**************//

  typedef  void (MeshDataTransferUtilities::*TransferElementalValuesToNodes)(const MeshDataTransferUtilities::TransferParameters&, ModelPart&);
  typedef  void (MeshDataTransferUtilities::*TransferNodalValuesToElements)(const MeshDataTransferUtilities::TransferParameters&, ModelPart&);
  typedef  void (MeshDataTransferUtilities::*TransferBoundaryData)(const MeshDataTransferUtilities::TransferParameters&, ModelPart&);

  TransferElementalValuesToNodes   TransferElementsToNodes = &MeshDataTransferUtilities::TransferElementalValuesToNodes;
  TransferNodalValuesToElements    TransferNodesToElements = &MeshDataTransferUtilities::TransferNodalValuesToElements;
  TransferBoundaryData             TransferDataToBoundary  = &MeshDataTransferUtilities::TransferBoundaryData;

  class_<MeshDataTransferUtilities>(m,"MeshDataTransferUtilities")
      .def(init<>())
      .def("TransferElementalValuesToNodes", TransferElementsToNodes)
      .def("TransferNodalValuesToElements", TransferNodesToElements)
      .def("TransferNodalValuesToElementsOnThreshold", TransferNodesToElementsOnThreshold)
      .def("TransferBoundaryData", TransferDataToBoundary)
      .def_readonly_static("NODE_TO_ELEMENT",&MeshDataTransferUtilities::NODE_TO_ELEMENT)
      .def_readonly_static("ELEMENT_TO_NODE",&MeshDataTransferUtilities::ELEMENT_TO_NODE)
      .def_readonly_static("ELEMENT_TO_ELEMENT",&MeshDataTransferUtilities::ELEMENT_TO_ELEMENT)
      .def_readonly_static("INITIALIZE_MASTER_CONDITION",&MeshDataTransferUtilities::INITIALIZE_MASTER_CONDITION)
      .def_readonly_static("MASTER_ELEMENT_TO_MASTER_CONDITION",&MeshDataTransferUtilities::MASTER_ELEMENT_TO_MASTER_CONDITION)
      ;



  // Remeshing modeler information parameters
  class_< MeshDataTransferUtilities::TransferParameters, typename MeshDataTransferUtilities::TransferParameters::Pointer>
      (m,"TransferParameters")
      .def(init<>())    
      .def("Set",&MeshDataTransferUtilities::TransferParameters::Set)
      .def("Reset",&MeshDataTransferUtilities::TransferParameters::Reset)
      .def("SetOptions",&MeshDataTransferUtilities::TransferParameters::SetOptions)
      .def("GetOptions",&MeshDataTransferUtilities::TransferParameters::GetOptions)
      .def("SetVariable",SetDoubleVariable)
      .def("SetVariable",SetArray1DVariable)
      .def("SetVariable",SetVectorVariable)
      .def("SetVariable",SetMatrixVariable)
      ;

  //***************MODELER UTILITIES**************//

  // Remeshing modeler information parameters
  class_< ModelerUtilities::MeshingInfoParameters, typename ModelerUtilities::MeshingInfoParameters::Pointer>
      (m,"MeshingInfoParameters")
      .def(init<>())    
      .def("Initialize",&ModelerUtilities::MeshingInfoParameters::Initialize)
      .def("CheckMechanicalSmoothing",&ModelerUtilities::MeshingInfoParameters::CheckMechanicalSmoothing)
      .def("SetNumberOfNodes",&ModelerUtilities::MeshingInfoParameters::SetNumberOfNodes)
      .def("SetNumberOfElements",&ModelerUtilities::MeshingInfoParameters::SetNumberOfElements)
      .def("SetNumberOfConditions",&ModelerUtilities::MeshingInfoParameters::SetNumberOfConditions)
      .def("GetNumberOfNodes",&ModelerUtilities::MeshingInfoParameters::GetNumberOfNodes)
      .def("GetNumberOfElements",&ModelerUtilities::MeshingInfoParameters::GetNumberOfElements)
      .def("GetNumberOfConditions",&ModelerUtilities::MeshingInfoParameters::GetNumberOfConditions)
      .def("SetNumberOfNewNodes",&ModelerUtilities::MeshingInfoParameters::SetNumberOfNewNodes)
      .def("SetNumberOfNewElements",&ModelerUtilities::MeshingInfoParameters::SetNumberOfNewElements)
      .def("SetNumberOfNewConditions",&ModelerUtilities::MeshingInfoParameters::SetNumberOfNewConditions)
      .def("SetInitialMeshVolume",&ModelerUtilities::MeshingInfoParameters::SetInitialMeshVolume)
      .def("GetInitialMeshVolume",&ModelerUtilities::MeshingInfoParameters::GetInitialMeshVolume)
      .def("GetInsertedNodes",&ModelerUtilities::MeshingInfoParameters::GetInsertedNodes)
      ;

  // Remeshing modeler refining parameters
  class_< ModelerUtilities::RefiningParameters, typename ModelerUtilities::RefiningParameters::Pointer>
      (m,"RefiningParameters")
      .def(init<>())
      .def("Initialize",&ModelerUtilities::RefiningParameters::Initialize)
      .def("SetRefiningOptions",&ModelerUtilities::RefiningParameters::SetRefiningOptions)
      .def("SetRemovingOptions",&ModelerUtilities::RefiningParameters::SetRemovingOptions)
      .def("SetAlphaParameter",&ModelerUtilities::RefiningParameters::SetAlphaParameter)
      .def("SetMeanVolume",&ModelerUtilities::RefiningParameters::SetMeanVolume)
      .def("GetMeanVolume",&ModelerUtilities::RefiningParameters::GetMeanVolume)
      .def("SetCriticalRadius",&ModelerUtilities::RefiningParameters::SetCriticalRadius)
      .def("SetInitialRadius",&ModelerUtilities::RefiningParameters::SetInitialRadius)
      .def("SetCriticalSide",&ModelerUtilities::RefiningParameters::SetCriticalSide)
      .def("SetRefiningParameters",SetRefiningParameters)
      .def("SetRefiningBox",&ModelerUtilities::RefiningParameters::SetRefiningBox)
      .def("SetReferenceThreshold",&ModelerUtilities::RefiningParameters::SetReferenceThreshold)
      .def("SetReferenceError",&ModelerUtilities::RefiningParameters::SetReferenceError)
      .def("SetThresholdVariable",SetThresholdVariable)
      .def("SetErrorVariable",SetErrorVariable)
      .def("GetThresholdVariable",GetThresholdVariable)
      .def("GetReferenceThreshold",&ModelerUtilities::RefiningParameters::GetReferenceThreshold)
      .def("GetRefiningOptions",&ModelerUtilities::RefiningParameters::GetRefiningOptions)
      .def("GetRemovingOptions",&ModelerUtilities::RefiningParameters::GetRemovingOptions)
      ;

  // Remeshing modeler remeshing parameters
  class_< ModelerUtilities::MeshingParameters, typename ModelerUtilities::MeshingParameters::Pointer>
      (m,"MeshingParameters")
      .def(init<>())
      .def("Initialize",&ModelerUtilities::MeshingParameters::Initialize)
      .def("Set",&ModelerUtilities::MeshingParameters::Set)
      .def("Reset",&ModelerUtilities::MeshingParameters::Reset)
      .def("SetOptions",&ModelerUtilities::MeshingParameters::SetOptions)
      .def("SetSubModelPartName",&ModelerUtilities::MeshingParameters::SetSubModelPartName)
      .def("SetExecutionOptions",&ModelerUtilities::MeshingParameters::SetExecutionOptions)
      .def("SetTessellationFlags",&ModelerUtilities::MeshingParameters::SetTessellationFlags)
      .def("SetTessellationInfo",&ModelerUtilities::MeshingParameters::SetTessellationInfo)
      .def("SetAlphaParameter",&ModelerUtilities::MeshingParameters::SetAlphaParameter)
      .def("SetOffsetFactor",&ModelerUtilities::MeshingParameters::SetOffsetFactor)
      .def("SetInfoParameters",&ModelerUtilities::MeshingParameters::SetInfoParameters)
      .def("SetRefiningParameters",&ModelerUtilities::MeshingParameters::SetRefiningParameters)
      .def("SetProperties",&ModelerUtilities::MeshingParameters::SetProperties)
      .def("SetMeshingBox",&ModelerUtilities::MeshingParameters::SetMeshingBox)
      .def("SetTransferParameters",&ModelerUtilities::MeshingParameters::SetTransferParameters)
      .def("SetTransferVariable",&ModelerUtilities::MeshingParameters::SetTransferVariable)
      .def("SetReferenceElement",SetReferenceElement)
      .def("SetReferenceCondition",SetReferenceCondition)
      .def("GetInfoParameters",&ModelerUtilities::MeshingParameters::GetInfoParameters)
      .def("GetTransferParameters",&ModelerUtilities::MeshingParameters::GetTransferParameters)
      .def("GetRefiningParameters",&ModelerUtilities::MeshingParameters::GetRefiningParameters)
      .def("GetSubModelPartName",&ModelerUtilities::MeshingParameters::GetSubModelPartName)
      .def("GetOptions",&ModelerUtilities::MeshingParameters::GetOptions)
      .def("InitializeMeshing",&ModelerUtilities::MeshingParameters::InitializeMeshing)
      .def("FinalizeMeshing",&ModelerUtilities::MeshingParameters::FinalizeMeshing)
      ;

}

}  // namespace Python.

} // Namespace Kratos

