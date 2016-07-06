//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes 

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

//Application includes
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
    void TransferNodesToElementsOnThreshold( MeshDataTransferUtilities& rMeshDataTransfer, MeshDataTransferUtilities::TransferParameters& rTransferParameters, const Variable<double>& rVariable, double Value, ModelPart& rModelPart, int MeshId)
    {
      rMeshDataTransfer.TransferNodalValuesToElements( rTransferParameters, rVariable, Value, rModelPart, MeshId );
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


    void  AddCustomUtilitiesToPython()
    {

      using namespace boost::python;

      
      typedef UblasSpace<double, CompressedMatrix, Vector>    SparseSpaceType;
      typedef UblasSpace<double, Matrix, Vector>               LocalSpaceType;
      typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

      typedef Scheme< SparseSpaceType, LocalSpaceType >            SchemeType;
      typedef SchemeType::Pointer                           SchemePointerType;

      typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >   BuilderAndSolverType;
      typedef BuilderAndSolverType::Pointer        BuilderAndSolverPointerType;

      //***************DOMAIN SET**************//
      class_< ModelerUtilities, boost::noncopyable > ("ModelerUtilities", init<>())
	.def("SetDomainLabels",&ModelerUtilities::SetDomainLabels)
	.def_readonly("REMESH",&ModelerUtilities::REMESH)
	.def_readonly("REFINE",&ModelerUtilities::REFINE)
	.def_readonly("RECONNECT",&ModelerUtilities::RECONNECT)
	.def_readonly("TRANSFER",&ModelerUtilities::TRANSFER)
	.def_readonly("CONSTRAINED",&ModelerUtilities::CONSTRAINED)
	.def_readonly("CONTACT_SEARCH",&ModelerUtilities::CONTACT_SEARCH)
	.def_readonly("MESH_SMOOTHING",&ModelerUtilities::MESH_SMOOTHING)
	.def_readonly("VARIABLES_SMOOTHING",&ModelerUtilities::VARIABLES_SMOOTHING)

	.def_readonly("REMOVE_NODES",&ModelerUtilities::REMOVE_NODES)
	.def_readonly("REMOVE_NODES_ON_DISTANCE",&ModelerUtilities::REMOVE_NODES_ON_DISTANCE)
	.def_readonly("REMOVE_NODES_ON_ERROR",&ModelerUtilities::REMOVE_NODES_ON_ERROR)
	.def_readonly("REMOVE_NODES_ON_THRESHOLD",&ModelerUtilities::REMOVE_NODES_ON_THRESHOLD)
	.def_readonly("REMOVE_BOUNDARY_NODES",&ModelerUtilities::REMOVE_BOUNDARY_NODES)
	.def_readonly("REMOVE_BOUNDARY_NODES_ON_DISTANCE",&ModelerUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE)
	.def_readonly("REMOVE_BOUNDARY_NODES_ON_ERROR",&ModelerUtilities::REMOVE_BOUNDARY_NODES_ON_ERROR)
	.def_readonly("REMOVE_BOUNDARY_NODES_ON_THRESHOLD",&ModelerUtilities::REMOVE_BOUNDARY_NODES_ON_THRESHOLD)

	.def_readonly("REFINE_ADD_NODES",&ModelerUtilities::REFINE_ADD_NODES)
	.def_readonly("REFINE_INSERT_NODES",&ModelerUtilities::REFINE_INSERT_NODES)
	.def_readonly("REFINE_ELEMENTS",&ModelerUtilities::REFINE_ELEMENTS)
	.def_readonly("REFINE_ELEMENTS_ON_DISTANCE",&ModelerUtilities::REFINE_ELEMENTS_ON_DISTANCE)
	.def_readonly("REFINE_ELEMENTS_ON_ERROR",&ModelerUtilities::REFINE_ELEMENTS_ON_ERROR)
	.def_readonly("REFINE_ELEMENTS_ON_THRESHOLD",&ModelerUtilities::REFINE_ELEMENTS_ON_THRESHOLD)
	.def_readonly("REFINE_BOUNDARY",&ModelerUtilities::REFINE_BOUNDARY)
	.def_readonly("REFINE_BOUNDARY_ON_DISTANCE",&ModelerUtilities::REFINE_BOUNDARY_ON_DISTANCE)
	.def_readonly("REFINE_BOUNDARY_ON_ERROR",&ModelerUtilities::REFINE_BOUNDARY_ON_ERROR)
	.def_readonly("REFINE_BOUNDARY_ON_THRESHOLD",&ModelerUtilities::REFINE_BOUNDARY_ON_THRESHOLD)

	.def_readonly("INITIALIZE_MESHER_INPUT",&ModelerUtilities::INITIALIZE_MESHER_INPUT)
	.def_readonly("FINALIZE_MESHER_INPUT",&ModelerUtilities::FINALIZE_MESHER_INPUT)
	.def_readonly("TRANSFER_KRATOS_NODES_TO_MESHER",&ModelerUtilities::TRANSFER_KRATOS_NODES_TO_MESHER)
	.def_readonly("TRANSFER_KRATOS_ELEMENTS_TO_MESHER",&ModelerUtilities::TRANSFER_KRATOS_ELEMENTS_TO_MESHER)
	.def_readonly("TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER",&ModelerUtilities::TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER)
	.def_readonly("TRANSFER_KRATOS_FACES_TO_MESHER",&ModelerUtilities::TRANSFER_KRATOS_FACES_TO_MESHER)

	.def_readonly("SELECT_TESSELLATION_ELEMENTS",&ModelerUtilities::SELECT_TESSELLATION_ELEMENTS)
	.def_readonly("KEEP_ISOLATED_NODES",&ModelerUtilities::KEEP_ISOLATED_NODES)
	;
        
        

      //***************NORMALS**************//

      // This is required to recognize the different overloads 
      typedef  void (BoundaryNormalsCalculationUtilities::*CalculateMeshBoundaryNormals)(ModelPart&, int, int);
      typedef  void (BoundaryNormalsCalculationUtilities::*CalculateMeshUnitBoundaryNormals)(ModelPart&, int, int); 

      CalculateMeshBoundaryNormals          CalculateMeshNormals     = &BoundaryNormalsCalculationUtilities::CalculateMeshBoundaryNormals;
      CalculateMeshUnitBoundaryNormals      CalculateMeshUnitNormals = &BoundaryNormalsCalculationUtilities::CalculateMeshUnitBoundaryNormals;
      
      class_<BoundaryNormalsCalculationUtilities > ("BoundaryNormalsCalculation", init<>())
	.def("CalculateBoundaryNormals", &BoundaryNormalsCalculationUtilities::CalculateBoundaryNormals)
	.def("CalculateBoundaryUnitNormals", &BoundaryNormalsCalculationUtilities::CalculateUnitBoundaryNormals)
	.def("CalculateMeshBoundaryNormals", CalculateMeshNormals)
	.def("CalculateMeshBoundaryUnitNormals", CalculateMeshUnitNormals)
	;
      

      //***************TRANSFER UTILITIES**************//

      typedef  void (MeshDataTransferUtilities::*TransferElementalValuesToNodes)(const MeshDataTransferUtilities::TransferParameters&, ModelPart&, ModelPart::IndexType);
      typedef  void (MeshDataTransferUtilities::*TransferNodalValuesToElements)(const MeshDataTransferUtilities::TransferParameters&, ModelPart&, ModelPart::IndexType);

      TransferElementalValuesToNodes   TransferElementsToNodes  = &MeshDataTransferUtilities::TransferElementalValuesToNodes;
      TransferNodalValuesToElements    TransferNodesToElements  = &MeshDataTransferUtilities::TransferNodalValuesToElements;


      class_<MeshDataTransferUtilities> ("MeshDataTransferUtilities", init<>())
	.def("TransferElementalValuesToNodes", TransferElementsToNodes)
	.def("TransferNodalValuesToElements", TransferNodesToElements)
	.def("TransferNodalValuesToElementsOnThreshold", TransferNodesToElementsOnThreshold)
	;



      // Remeshing modeler information parameters
      class_< MeshDataTransferUtilities::TransferParameters, MeshDataTransferUtilities::TransferParameters::Pointer, boost::noncopyable>
	("TransferParameters", init<>() )    
	.def("Set",&MeshDataTransferUtilities::TransferParameters::Set)
	.def("Reset",&MeshDataTransferUtilities::TransferParameters::Reset)
	.def("SetOptions",&MeshDataTransferUtilities::TransferParameters::SetOptions)
	.def("SetVariable",SetDoubleVariable)
	.def("SetVariable",SetArray1DVariable)
	.def("SetVariable",SetVectorVariable)
	.def("SetVariable",SetMatrixVariable)
	;

      //***************MODELER UTILITIES**************//

      // Remeshing modeler information parameters
      class_< ModelerUtilities::InfoParameters, ModelerUtilities::InfoParameters::Pointer, boost::noncopyable>
	("InfoParameters", init<>() )    
	.def("Initialize",&ModelerUtilities::InfoParameters::Initialize)
	.def("CheckMechanicalSmoothing",&ModelerUtilities::InfoParameters::CheckMechanicalSmoothing)
	.def("SetNumberOfNodes",&ModelerUtilities::InfoParameters::SetNumberOfNodes)
	.def("SetNumberOfElements",&ModelerUtilities::InfoParameters::SetNumberOfElements)
	.def("SetNumberOfConditions",&ModelerUtilities::InfoParameters::SetNumberOfConditions)
	.def("SetNumberOfNewNodes",&ModelerUtilities::InfoParameters::SetNumberOfNewNodes)
	.def("SetNumberOfNewElements",&ModelerUtilities::InfoParameters::SetNumberOfNewElements)
	.def("SetNumberOfNewConditions",&ModelerUtilities::InfoParameters::SetNumberOfNewConditions)

	;

      // Remeshing modeler refining parameters
      class_< ModelerUtilities::RefiningParameters, ModelerUtilities::RefiningParameters::Pointer, boost::noncopyable>
	("RefiningParameters", init<>() )
	.def("Initialize",&ModelerUtilities::RefiningParameters::Initialize)
	.def("SetRefiningOptions",&ModelerUtilities::RefiningParameters::SetRefiningOptions)
	.def("SetRemovingOptions",&ModelerUtilities::RefiningParameters::SetRemovingOptions)
	.def("SetExecutionOptions",&ModelerUtilities::RefiningParameters::SetExecutionOptions)
	.def("SetAlphaParameter",&ModelerUtilities::RefiningParameters::SetAlphaParameter)
	.def("SetCriticalRadius",&ModelerUtilities::RefiningParameters::SetCriticalRadius)
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
      class_< ModelerUtilities::MeshingParameters, ModelerUtilities::MeshingParameters::Pointer, boost::noncopyable>
	("MeshingParameters", init<>() )
	.def("Initialize",&ModelerUtilities::MeshingParameters::Initialize)
	.def("Set",&ModelerUtilities::MeshingParameters::Set)
	.def("Reset",&ModelerUtilities::MeshingParameters::Reset)
	.def("SetOptions",&ModelerUtilities::MeshingParameters::SetOptions)
	.def("SetExecutionOptions",&ModelerUtilities::MeshingParameters::SetExecutionOptions)
	.def("SetTessellationFlags",&ModelerUtilities::MeshingParameters::SetTessellationFlags)
	.def("SetTessellationInfo",&ModelerUtilities::MeshingParameters::SetTessellationInfo)
	.def("SetAlphaParameter",&ModelerUtilities::MeshingParameters::SetAlphaParameter)
	.def("SetOffsetFactor",&ModelerUtilities::MeshingParameters::SetOffsetFactor)
	.def("SetInfoParameters",&ModelerUtilities::MeshingParameters::SetInfoParameters)
	.def("SetRefiningParameters",&ModelerUtilities::MeshingParameters::SetRefiningParameters)
	.def("SetMeshingBox",&ModelerUtilities::MeshingParameters::SetMeshingBox)
	.def("SetTransferParameters",&ModelerUtilities::MeshingParameters::SetTransferParameters)
	.def("SetTransferVariable",&ModelerUtilities::MeshingParameters::SetTransferVariable)
	.def("SetReferenceElement",SetReferenceElement)
	.def("SetReferenceCondition",SetReferenceCondition)
	.def("GetInfoParameters",&ModelerUtilities::MeshingParameters::GetInfoParameters)
	.def("GetRefiningParameters",&ModelerUtilities::MeshingParameters::GetRefiningParameters)
	.def("GetOptions",&ModelerUtilities::MeshingParameters::GetOptions)
	;


    }

  }  // namespace Python.

} // Namespace Kratos

