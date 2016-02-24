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
	.def_readonly("CONSTRAINED",&ModelerUtilities::CONSTRAINED)
	.def_readonly("CONSTACT_SEARCH",&ModelerUtilities::CONTACT_SEARCH)
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
	;
        
      
      class_< MeshDataTransferUtilities, boost::noncopyable > ("MeshDataTransferUtilities", init<>())
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
	;

      // Remeshing modeler refining parameters
      class_< ModelerUtilities::RefiningParameters, ModelerUtilities::RefiningParameters::Pointer, boost::noncopyable>
	("RefiningParameters", init<>() )
	.def("Initialize",&ModelerUtilities::RefiningParameters::Initialize)
	.def("SetRefiningOptions",&ModelerUtilities::RefiningParameters::SetRefiningOptions)
	.def("SetRemovingOptions",&ModelerUtilities::RefiningParameters::SetRemovingOptions)
	.def("SetAlphaParameter",&ModelerUtilities::RefiningParameters::SetAlphaParameter)
	.def("SetCriticalRadius",&ModelerUtilities::RefiningParameters::SetCriticalRadius)
	.def("SetCriticalSide",&ModelerUtilities::RefiningParameters::SetCriticalSide)
	.def("SetRefiningParameters",SetRefiningParameters)
	.def("SetRefiningBox",&ModelerUtilities::RefiningParameters::SetRefiningBox)
	.def("SetReferenceThreshold",&ModelerUtilities::RefiningParameters::SetReferenceThreshold)
	.def("SetReferenceError",&ModelerUtilities::RefiningParameters::SetReferenceError)
	.def("SetThresholdVariable",SetThresholdVariable)
	.def("SetErrorVariable",SetErrorVariable)
	;

      // Remeshing modeler remeshing parameters
      class_< ModelerUtilities::MeshingParameters, ModelerUtilities::MeshingParameters::Pointer, boost::noncopyable>
	("MeshingParameters", init<>() )
	.def("Initialize",&ModelerUtilities::MeshingParameters::Initialize)
	.def("Set",&ModelerUtilities::MeshingParameters::Set)
	.def("Reset",&ModelerUtilities::MeshingParameters::Reset)
	.def("SetOptions",&ModelerUtilities::MeshingParameters::SetOptions)
	.def("SetAlphaParameter",&ModelerUtilities::MeshingParameters::SetAlphaParameter)
	.def("SetOffsetFactor",&ModelerUtilities::MeshingParameters::SetOffsetFactor)
	.def("SetInfoParameters",&ModelerUtilities::MeshingParameters::SetInfoParameters)
	.def("SetRefiningParameters",&ModelerUtilities::MeshingParameters::SetRefiningParameters)
	.def("SetMeshingBox",&ModelerUtilities::MeshingParameters::SetMeshingBox)
	.def("SetTransferParameters",&ModelerUtilities::MeshingParameters::SetTransferParameters)
	.def("SetTransferVariable",&ModelerUtilities::MeshingParameters::SetTransferVariable)
	.def("SetReferenceElement",SetReferenceElement)
	.def("SetReferenceCondition",SetReferenceCondition)
	;


    }

  }  // namespace Python.

} // Namespace Kratos

