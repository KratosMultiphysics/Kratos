//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

//Utilities
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"

namespace Kratos
{

namespace Python
{

// refining methods

void SetRefiningParameters(MesherUtilities::RefiningParameters& rRefiningParameters, double AlphaParameter, double CriticalRadius, double CriticalSide)
{
  rRefiningParameters.SetParameters(AlphaParameter, CriticalRadius, CriticalSide);
}


void SetThresholdVariable(MesherUtilities::RefiningParameters& rRefiningParameters, const Variable<double>& rVariable)
{
  rRefiningParameters.SetThresholdVariable(rVariable);
}

void SetErrorVariable(MesherUtilities::RefiningParameters& rRefiningParameters, const Variable<double>& rVariable)
{
  rRefiningParameters.SetErrorVariable(rVariable);
}

const Variable<double>& GetThresholdVariable(MesherUtilities::RefiningParameters& rRefiningParameters)
{
  return rRefiningParameters.GetThresholdVariable();
}
// remeshing methods

void SetReferenceElement(MesherUtilities::MeshingParameters& rMeshingParameters, char* ElementName)
{
  rMeshingParameters.SetReferenceElement(KratosComponents<Element>::Get(ElementName));
}

void SetReferenceCondition(MesherUtilities::MeshingParameters& rMeshingParameters, char* ConditionName)
{
  rMeshingParameters.SetReferenceCondition(KratosComponents<Condition>::Get(ConditionName));
}


// transfer methods
void TransferNodesToElementsOnThreshold( MeshDataTransferUtilities& rMeshDataTransfer, MeshDataTransferUtilities::TransferParameters& rTransferParameters, MesherUtilities::RefiningParameters& rRefiningParameters, ModelPart& rModelPart)
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

  namespace py = pybind11;


  //***************DOMAIN SET**************//
  py::class_<MesherUtilities>(m,"MesherUtilities")
      .def(py::init<>())
      .def("SetModelPartNameToElements",&MesherUtilities::SetModelPartNameToElements)
      .def("SetModelPartNameToConditions",&MesherUtilities::SetModelPartNameToConditions)
      .def("SetModelPartNameToNodes",&MesherUtilities::SetModelPartNameToNodes)
      .def("CheckCriticalRadius",&MesherUtilities::CheckCriticalRadius)
      .def("ComputeModelPartVolume",&MesherUtilities::ComputeModelPartVolume)

      .def_readonly_static("REMESH",&MesherUtilities::REMESH)
      .def_readonly_static("REFINE",&MesherUtilities::REFINE)
      .def_readonly_static("RECONNECT",&MesherUtilities::RECONNECT)
      .def_readonly_static("TRANSFER",&MesherUtilities::TRANSFER)
      .def_readonly_static("CONSTRAINED",&MesherUtilities::CONSTRAINED)
      .def_readonly_static("CONTACT_SEARCH",&MesherUtilities::CONTACT_SEARCH)
      .def_readonly_static("MESH_SMOOTHING",&MesherUtilities::MESH_SMOOTHING)
      .def_readonly_static("VARIABLES_SMOOTHING",&MesherUtilities::VARIABLES_SMOOTHING)

      .def_readonly_static("REMOVE_NODES",&MesherUtilities::REMOVE_NODES)
      .def_readonly_static("REMOVE_NODES_ON_DISTANCE",&MesherUtilities::REMOVE_NODES_ON_DISTANCE)
      .def_readonly_static("REMOVE_NODES_ON_ERROR",&MesherUtilities::REMOVE_NODES_ON_ERROR)
      .def_readonly_static("REMOVE_NODES_ON_THRESHOLD",&MesherUtilities::REMOVE_NODES_ON_THRESHOLD)
      .def_readonly_static("REMOVE_BOUNDARY_NODES",&MesherUtilities::REMOVE_BOUNDARY_NODES)
      .def_readonly_static("REMOVE_BOUNDARY_NODES_ON_DISTANCE",&MesherUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE)
      .def_readonly_static("REMOVE_BOUNDARY_NODES_ON_ERROR",&MesherUtilities::REMOVE_BOUNDARY_NODES_ON_ERROR)
      .def_readonly_static("REMOVE_BOUNDARY_NODES_ON_THRESHOLD",&MesherUtilities::REMOVE_BOUNDARY_NODES_ON_THRESHOLD)

      .def_readonly_static("REFINE_ADD_NODES",&MesherUtilities::REFINE_ADD_NODES)
      .def_readonly_static("REFINE_INSERT_NODES",&MesherUtilities::REFINE_INSERT_NODES)
      .def_readonly_static("REFINE_ELEMENTS",&MesherUtilities::REFINE_ELEMENTS)
      .def_readonly_static("REFINE_ELEMENTS_ON_DISTANCE",&MesherUtilities::REFINE_ELEMENTS_ON_DISTANCE)
      .def_readonly_static("REFINE_ELEMENTS_ON_ERROR",&MesherUtilities::REFINE_ELEMENTS_ON_ERROR)
      .def_readonly_static("REFINE_ELEMENTS_ON_THRESHOLD",&MesherUtilities::REFINE_ELEMENTS_ON_THRESHOLD)
      .def_readonly_static("REFINE_BOUNDARY",&MesherUtilities::REFINE_BOUNDARY)
      .def_readonly_static("REFINE_BOUNDARY_ON_DISTANCE",&MesherUtilities::REFINE_BOUNDARY_ON_DISTANCE)
      .def_readonly_static("REFINE_BOUNDARY_ON_ERROR",&MesherUtilities::REFINE_BOUNDARY_ON_ERROR)
      .def_readonly_static("REFINE_BOUNDARY_ON_THRESHOLD",&MesherUtilities::REFINE_BOUNDARY_ON_THRESHOLD)

      .def_readonly_static("INITIALIZE_MESHER_INPUT",&MesherUtilities::INITIALIZE_MESHER_INPUT)
      .def_readonly_static("FINALIZE_MESHER_INPUT",&MesherUtilities::FINALIZE_MESHER_INPUT)
      .def_readonly_static("TRANSFER_KRATOS_NODES_TO_MESHER",&MesherUtilities::TRANSFER_KRATOS_NODES_TO_MESHER)
      .def_readonly_static("TRANSFER_KRATOS_ELEMENTS_TO_MESHER",&MesherUtilities::TRANSFER_KRATOS_ELEMENTS_TO_MESHER)
      .def_readonly_static("TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER",&MesherUtilities::TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER)
      .def_readonly_static("TRANSFER_KRATOS_FACES_TO_MESHER",&MesherUtilities::TRANSFER_KRATOS_FACES_TO_MESHER)

      .def_readonly_static("SELECT_TESSELLATION_ELEMENTS",&MesherUtilities::SELECT_TESSELLATION_ELEMENTS)
      .def_readonly_static("KEEP_ISOLATED_NODES",&MesherUtilities::KEEP_ISOLATED_NODES)
      .def_readonly_static("REFINE_WALL_CORNER",&MesherUtilities::REFINE_WALL_CORNER)
      ;


  //***************NORMALS**************//

  py::class_<BoundaryNormalsCalculationUtilities>(m,"BoundaryNormalsCalculation")
      .def(py::init<>())
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

  py::class_<MeshDataTransferUtilities>(m,"MeshDataTransferUtilities")
      .def(py::init<>())
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



  // Remeshing mesher information parameters
  py::class_< MeshDataTransferUtilities::TransferParameters, typename MeshDataTransferUtilities::TransferParameters::Pointer>
      (m,"TransferParameters")
      .def(py::init<>())
      .def("Set",&MeshDataTransferUtilities::TransferParameters::Set)
      .def("Reset",&MeshDataTransferUtilities::TransferParameters::Reset)
      .def("SetOptions",&MeshDataTransferUtilities::TransferParameters::SetOptions)
      .def("GetOptions",&MeshDataTransferUtilities::TransferParameters::GetOptions)
      .def("SetVariable",SetDoubleVariable)
      .def("SetVariable",SetArray1DVariable)
      .def("SetVariable",SetVectorVariable)
      .def("SetVariable",SetMatrixVariable)
      ;

  //***************MESHER UTILITIES**************//

  // Remeshing mesher information parameters
  py::class_< MesherUtilities::MeshingInfoParameters, typename MesherUtilities::MeshingInfoParameters::Pointer>
      (m,"MeshingInfoParameters")
      .def(py::init<>())
      .def("Initialize",&MesherUtilities::MeshingInfoParameters::Initialize)
      .def("CheckMechanicalSmoothing",&MesherUtilities::MeshingInfoParameters::CheckMechanicalSmoothing)
      .def("SetNumberOfNodes",&MesherUtilities::MeshingInfoParameters::SetNumberOfNodes)
      .def("SetNumberOfElements",&MesherUtilities::MeshingInfoParameters::SetNumberOfElements)
      .def("SetNumberOfConditions",&MesherUtilities::MeshingInfoParameters::SetNumberOfConditions)
      .def("GetNumberOfNodes",&MesherUtilities::MeshingInfoParameters::GetNumberOfNodes)
      .def("GetNumberOfElements",&MesherUtilities::MeshingInfoParameters::GetNumberOfElements)
      .def("GetNumberOfConditions",&MesherUtilities::MeshingInfoParameters::GetNumberOfConditions)
      .def("SetNumberOfNewNodes",&MesherUtilities::MeshingInfoParameters::SetNumberOfNewNodes)
      .def("SetNumberOfNewElements",&MesherUtilities::MeshingInfoParameters::SetNumberOfNewElements)
      .def("SetNumberOfNewConditions",&MesherUtilities::MeshingInfoParameters::SetNumberOfNewConditions)
      .def("SetInitialMeshVolume",&MesherUtilities::MeshingInfoParameters::SetInitialMeshVolume)
      .def("GetInitialMeshVolume",&MesherUtilities::MeshingInfoParameters::GetInitialMeshVolume)
      .def("GetInsertedNodes",&MesherUtilities::MeshingInfoParameters::GetInsertedNodes)
      ;

  // Remeshing mesher refining parameters
  py::class_< MesherUtilities::RefiningParameters, typename MesherUtilities::RefiningParameters::Pointer>
      (m,"RefiningParameters")
      .def(py::init<>())
      .def("Initialize",&MesherUtilities::RefiningParameters::Initialize)
      .def("SetRefiningOptions",&MesherUtilities::RefiningParameters::SetRefiningOptions)
      .def("SetRemovingOptions",&MesherUtilities::RefiningParameters::SetRemovingOptions)
      .def("SetAlphaParameter",&MesherUtilities::RefiningParameters::SetAlphaParameter)
      .def("SetMeanVolume",&MesherUtilities::RefiningParameters::SetMeanVolume)
      .def("GetMeanVolume",&MesherUtilities::RefiningParameters::GetMeanVolume)
      .def("SetCriticalRadius",&MesherUtilities::RefiningParameters::SetCriticalRadius)
      .def("SetInitialRadius",&MesherUtilities::RefiningParameters::SetInitialRadius)
      .def("SetCriticalSide",&MesherUtilities::RefiningParameters::SetCriticalSide)
      .def("SetRefiningParameters",SetRefiningParameters)
      .def("SetRefiningBox",&MesherUtilities::RefiningParameters::SetRefiningBox)
      .def("SetReferenceThreshold",&MesherUtilities::RefiningParameters::SetReferenceThreshold)
      .def("SetReferenceError",&MesherUtilities::RefiningParameters::SetReferenceError)
      .def("SetThresholdVariable",SetThresholdVariable)
      .def("SetErrorVariable",SetErrorVariable)
      .def("GetThresholdVariable",GetThresholdVariable)
      .def("GetReferenceThreshold",&MesherUtilities::RefiningParameters::GetReferenceThreshold)
      .def("GetRefiningOptions",&MesherUtilities::RefiningParameters::GetRefiningOptions)
      .def("GetRemovingOptions",&MesherUtilities::RefiningParameters::GetRemovingOptions)
      ;

  // Remeshing mesher remeshing parameters
  py::class_< MesherUtilities::MeshingParameters, typename MesherUtilities::MeshingParameters::Pointer>
      (m,"MeshingParameters")
      .def(py::init<>())
      .def("Initialize",&MesherUtilities::MeshingParameters::Initialize)
      .def("Set",&MesherUtilities::MeshingParameters::Set)
      .def("Reset",&MesherUtilities::MeshingParameters::Reset)
      .def("SetOptions",&MesherUtilities::MeshingParameters::SetOptions)
      .def("SetSubModelPartName",&MesherUtilities::MeshingParameters::SetSubModelPartName)
      .def("SetExecutionOptions",&MesherUtilities::MeshingParameters::SetExecutionOptions)
      .def("SetTessellationFlags",&MesherUtilities::MeshingParameters::SetTessellationFlags)
      .def("SetTessellationInfo",&MesherUtilities::MeshingParameters::SetTessellationInfo)
      .def("SetAlphaParameter",&MesherUtilities::MeshingParameters::SetAlphaParameter)
      .def("SetOffsetFactor",&MesherUtilities::MeshingParameters::SetOffsetFactor)
      .def("SetInfoParameters",&MesherUtilities::MeshingParameters::SetInfoParameters)
      .def("SetRefiningParameters",&MesherUtilities::MeshingParameters::SetRefiningParameters)
      .def("SetProperties",&MesherUtilities::MeshingParameters::SetProperties)
      .def("SetMeshingBox",&MesherUtilities::MeshingParameters::SetMeshingBox)
      .def("SetTransferParameters",&MesherUtilities::MeshingParameters::SetTransferParameters)
      .def("SetTransferVariable",&MesherUtilities::MeshingParameters::SetTransferVariable)
      .def("SetReferenceElement",SetReferenceElement)
      .def("SetReferenceCondition",SetReferenceCondition)
      .def("GetInfoParameters",&MesherUtilities::MeshingParameters::GetInfoParameters)
      .def("GetTransferParameters",&MesherUtilities::MeshingParameters::GetTransferParameters)
      .def("GetRefiningParameters",&MesherUtilities::MeshingParameters::GetRefiningParameters)
      .def("GetSubModelPartName",&MesherUtilities::MeshingParameters::GetSubModelPartName)
      .def("GetOptions",&MesherUtilities::MeshingParameters::GetOptions)
      .def("InitializeMeshing",&MesherUtilities::MeshingParameters::InitializeMeshing)
      .def("FinalizeMeshing",&MesherUtilities::MeshingParameters::FinalizeMeshing)
      .def("SetUseBoundingBox",&MesherUtilities::MeshingParameters::SetUseBoundingBox)
      .def("SetBoundingBoxLowerPoint",&MesherUtilities::MeshingParameters::SetBoundingBoxLowerPoint)
      .def("SetBoundingBoxUpperPoint",&MesherUtilities::MeshingParameters::SetBoundingBoxUpperPoint)
      .def("SetBoundingBoxTimeInterval",&MesherUtilities::MeshingParameters::SetBoundingBoxTimeInterval)
      .def("SetUseRefiningBox",&MesherUtilities::MeshingParameters::SetUseRefiningBox)
      .def("SetRefiningBoxMinimumPoint",&MesherUtilities::MeshingParameters::SetRefiningBoxMinimumPoint)
      .def("SetRefiningBoxMaximumPoint",&MesherUtilities::MeshingParameters::SetRefiningBoxMaximumPoint)
      .def("SetRefiningBoxTimeInterval",&MesherUtilities::MeshingParameters::SetRefiningBoxTimeInterval)
      .def("SetRefiningBoxMeshSize",&MesherUtilities::MeshingParameters::SetRefiningBoxMeshSize)
      ;

}

}  // namespace Python.

} // Namespace Kratos
