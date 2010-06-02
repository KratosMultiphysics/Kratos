/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-23 14:27:01 $
//   Revision:            $Revision: 1.20 $
//
//


#if !defined(KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED


// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/define.h"
#include "custom_utilities/deactivation_utility.h"
#include "custom_utilities/variable_transfer_utility.h"

#ifdef _OPENMP
#include "custom_utilities/parallel_variable_transfer_utility.h"
#endif

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/contact_utility.h"
#include "custom_utilities/restart_utility.h"
#include "custom_utilities/node_snapping_utility.h"
#include "custom_elements/rigid_body_3D.h"
#include "custom_utilities/output_utility.h"
#include "custom_utilities/smoothing_utility.h"


//#include "custom_utilities/detect_elements_utility.h"
#include "custom_utilities/intra_fracture_triangle_utility.h"
#include "custom_utilities/inter_fracture_triangle_utility.h"
#include "custom_utilities/inter_fracture_tetrahedra_utility.h"


namespace Kratos
{

    namespace Python
    {

        using namespace boost::python;
	
	void AddNewRigidBody3D( ModelPart& structural_model_part,
				ModelPart& skin_model_part,
				Variable<double>& rSelectionVariable,
    				double selection_value,
    				Node<3>::Pointer CenterNode,
				Element::PropertiesType::Pointer pProperties, 
				double nodal_mass,
    				Matrix& Inertia
				)
	{
		Geometry<Node<3> >::Pointer skin_nodes_geometry( new Geometry<Node<3> > ); ;
	  
		//selecting the nodes in the model part having rSelectionVariable==selection_value
		for(ModelPart::NodesContainerType::iterator it = skin_model_part.NodesBegin(); it!=skin_model_part.NodesEnd(); it++)
		{
			if(it->FastGetSolutionStepValue(rSelectionVariable) == selection_value)
				skin_nodes_geometry->push_back( *(it.base()) );
		}
		
		//creating a geometry containing the center node
		Geometry<Node<3> >::Pointer center_node_geometry( new Geometry<Node<3> > ) ;
		center_node_geometry->push_back( Node<3>::Pointer( CenterNode ) );
		
		unsigned int last_id = 1;
		if(structural_model_part.Elements().size() != 0)
			last_id = (structural_model_part.ElementsEnd()-1)->Id() + 1;

		array_1d<double,3> zero = ZeroVector(3);
		
		Element::Pointer new_el = RigidBody3D::Pointer( new  RigidBody3D(last_id,
					        	center_node_geometry, 
							pProperties,  
							skin_nodes_geometry, 
							nodal_mass, 
							Inertia,zero, zero  ) );
		
		structural_model_part.Elements().push_back(
				    new_el
				);
	}
				
	void AddNewRigidBodyAndSpring3D( ModelPart& structural_model_part,
				ModelPart& skin_model_part,
				Variable<double>& rSelectionVariable,
    				double selection_value,
    				Node<3>::Pointer CenterNode,
				Element::PropertiesType::Pointer pProperties, 
				double nodal_mass,
    				Matrix& Inertia,
				array_1d<double,3>& translational_stiffness,
				array_1d<double,3>& rotational_stiffness
				)
	{
		Geometry<Node<3> >::Pointer skin_nodes_geometry( new Geometry<Node<3> > ); ;
	  
		//selecting the nodes in the model part having rSelectionVariable==selection_value
		for(ModelPart::NodesContainerType::iterator it = skin_model_part.NodesBegin(); it!=skin_model_part.NodesEnd(); it++)
		{
			if(it->FastGetSolutionStepValue(rSelectionVariable) == selection_value)
				skin_nodes_geometry->push_back( *(it.base()) );
		}
		
		//creating a geometry containing the center node
		Geometry<Node<3> >::Pointer center_node_geometry( new Geometry<Node<3> > ) ;
		center_node_geometry->push_back( Node<3>::Pointer( CenterNode ) );
		
		unsigned int last_id = 1;
		if(structural_model_part.Elements().size() != 0)
			last_id = (structural_model_part.ElementsEnd()-1)->Id() + 1;

		array_1d<double,3> zero = ZeroVector(3);
		
		Element::Pointer new_el = RigidBody3D::Pointer( new  RigidBody3D(last_id,
					        	center_node_geometry, 
							pProperties,  
							skin_nodes_geometry, 
							nodal_mass, 
							Inertia,
							translational_stiffness, 
							rotational_stiffness  ) );
		
		structural_model_part.Elements().push_back(
				    new_el
				);
	}        



        void  AddCustomUtilitiesToPython()
        {
            class_<DeactivationUtility, boost::noncopyable >
                    ("DeactivationUtility",
                     init<>() )
                    .def("Deactivate", &DeactivationUtility::Deactivate )
                    .def("Reactivate", &DeactivationUtility::Reactivate )
                    .def("ReactivateStressFree", &DeactivationUtility::ReactivateStressFree )
                    .def("ReactivateAll", &DeactivationUtility::ReactivateAll )
                    .def("Initialize", &DeactivationUtility::Initialize )
                    ;
            
            class_<VariableTransferUtility, boost::noncopyable >
                    ("VariableTransferUtility",
                     init<>() )
                    .def("TransferNodalVariables", &VariableTransferUtility::TransferNodalVariables )
                    .def("TransferConstitutiveLawVariables", &VariableTransferUtility::TransferConstitutiveLawVariables )
                    .def("TransferInSituStress", &VariableTransferUtility::TransferInSituStress )
                    .def("InitializeModelPart", &VariableTransferUtility::InitializeModelPart )
                    .def("TransferVariablesToNodes", &VariableTransferUtility::DoubleTransferVariablesToNodes)
                    ;


#ifdef _OPENMP
           class_<ParallelVariableTransferUtility, boost::noncopyable >
                    ("ParallelVariableTransferUtility",
                     init<>() )
                    .def("TransferNodalVariables", &ParallelVariableTransferUtility::TransferNodalVariables )
                    .def("TransferConstitutiveLawVariables", &ParallelVariableTransferUtility::TransferConstitutiveLawVariables )
                    .def("TransferInSituStress", &ParallelVariableTransferUtility::TransferInSituStress )
                    .def("InitializeModelPart", &ParallelVariableTransferUtility::InitializeModelPart )
                    ;
#endif
            
            class_<ContactUtility, boost::noncopyable >
                    ("ContactUtility",
                     init<int>() )
                    .def("SetUpContactConditions", &ContactUtility::SetUpContactConditions )
                    .def("Update", &ContactUtility::Update )
                    .def("IsConverged", &ContactUtility::IsConverged )
                    .def("Clean", &ContactUtility::Clean )
                    ;

            class_<RestartUtility, boost::noncopyable >
                    ("RestartUtility",
                     init< std::string const& >() )
                    .def("ChangeFileName", &RestartUtility::ChangeFileName )
                    .def("StoreNodalVariables", &RestartUtility::StoreNodalVariables )
                    .def("WriteNodalVariables", &RestartUtility::WriteNodalVariables )
                    .def("StoreInSituStress", &RestartUtility::StoreInSituStress )
                    .def("WriteConstitutiveLawVariables", &RestartUtility::WriteConstitutiveLawVariables )
                    .def("StoreConstitutiveLawVariables", &RestartUtility::StoreConstitutiveLawVariables )
                    .def("WriteInSituStress", &RestartUtility::WriteInSituStress )
                    ;

            class_<NodeSnappingUtility, boost::noncopyable >
                    ("NodeSnappingUtility",
                     init<>() )
                    .def("MoveNode",&NodeSnappingUtility::MoveNode)
                    .def("AdjustNodes",&NodeSnappingUtility::AdjustNodes)
                    .def("AdjustToCircle",&NodeSnappingUtility::AdjustToCircle)
                    .def("AdjustToCylinder",&NodeSnappingUtility::AdjustToCylinder)
                    .def("AdjustToClosedCylinder",&NodeSnappingUtility::AdjustToClosedCylinder)
                    .def("IdentifyInsideElements",&NodeSnappingUtility::IdentifyInsideElements)
                    .def("SetInsituStress", &NodeSnappingUtility::SetInsituStress)
                    .def("ExtractCapNodes", &NodeSnappingUtility::ExtractCapNodes)
                    .def("TestElements", &NodeSnappingUtility::TestElements)
                    ;
            
            class_<OutputUtility, boost::noncopyable >
                    ("OutputUtility",
                     init<>() )
                    .def("GetStrain",&OutputUtility::GetStrain)
                    .def("GetStress",&OutputUtility::GetStress)
                    .def("GetInternalVariables",&OutputUtility::GetInternalVariables)
                    ;

	    
	    def("AddNewRigidBody3D", AddNewRigidBody3D);
	    def("AddNewRigidBodyAndSpring3D", AddNewRigidBodyAndSpring3D);
             ; 
      
/*
            class_<Detect_Elements_And_Nodes, boost::noncopyable >
                    ("DetectElementsAndNodes", init<ModelPart&, int >() )
		    .def("DetectNode",              &Detect_Elements_And_Nodes::Detect_Node_To_Be_Splitted)
                    .def("DetectElements",          &Detect_Elements_And_Nodes::Detect_Elements_To_Be_Splitted)
                    .def("CalculateMapFailure",     &Detect_Elements_And_Nodes::Calculate_Map_Failure) 
                    .def("Finalize",                &Detect_Elements_And_Nodes::Finalize)
                    ;
*/
            class_<Smoothing_Utility, boost::noncopyable >
                    ("SmoothingUtility", init<ModelPart&, int >() )
                    .def("WeightedRecoveryGradients", &Smoothing_Utility::WeightedRecoveryGradients) // for matrices
                    .def("DoubleWeightedRecoveryGradients", &Smoothing_Utility::DoubleWeightedRecoveryGradients) // for doubles
                    .def("InterpolatedRecoveryGradients", &Smoothing_Utility::InterpolatedRecoveryGradients)
                    .def("SettingNodalValues", &Smoothing_Utility::SettingNodalValues) 
                    .def("RecomputeValuesForNewMesh", &Smoothing_Utility::Recompute_Values_For_New_Mesh)
                    .def("Finalize", &Smoothing_Utility::Finalize)
                    .def("SettingNodalValues", &Smoothing_Utility::SettingNodalValues)   
                    ;


            class_<Intra_Fracture_Triangle, boost::noncopyable >
                    ("IntraFractureTriangle", init<ModelPart&, int >() )
		    .def("DetectAndSplitElements",              &Intra_Fracture_Triangle::Detect_And_Split_Elements)
                    ;

            class_<Inter_Fracture_Triangle, boost::noncopyable >
                    ("InterFractureTriangle", init<ModelPart&, int >() )
		    .def("DetectAndSplitElements",              &Inter_Fracture_Triangle::Detect_And_Split_Elements)
                    ;

             class_<Inter_Fracture_Tetrahedra, boost::noncopyable >
                    ("InterFractureTetrahedra", init<ModelPart&, int >() )
		    .def("DetectAndSplitElements",              &Inter_Fracture_Tetrahedra::Detect_And_Split_Elements)
                    ;  


    
        }
    }  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED  defined 
