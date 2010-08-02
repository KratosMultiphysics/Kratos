/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-10-13 07:04:58 $
//   Revision:            $Revision: 1.5 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// Project includes
#include "includes/define.h"
#include "includes/mesh.h"
#include "includes/properties.h"
#include "includes/element.h"
#include "includes/condition.h"

#include "python/add_mesh_to_python.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"
#include "python/solution_step_variable_indexing_python.h"

namespace Kratos
{
	namespace Python
	{
		using namespace boost::python;


	  template<class TVariableType> inline
	  void NodeFix(Node<3>& rThisNode, TVariableType const& rThisVariable)
	  {
	    rThisNode.Fix(rThisVariable);
	  }

	  template<class TVariableType> inline
	  void NodeFree(Node<3>& rThisNode, TVariableType const& rThisVariable)
	  {
	    rThisNode.Free(rThisVariable);
	  }

	  template<class TVariableType> inline
	  bool NodeIsFixed(Node<3>& rThisNode, TVariableType const& rThisVariable)
	  {
	    return rThisNode.IsFixed(rThisVariable);
	  }

	  template<class TVariableType> inline
	  bool NodeHasDofFor(Node<3>& rThisNode, TVariableType const& rThisVariable)
	  {
	    return rThisNode.HasDofFor(rThisVariable);
	  }

	  template<class TVariableType> inline
	  void NodeAddDof(Node<3>& rThisNode, TVariableType const& rThisVariable)
	  {
		rThisNode.pAddDof(rThisVariable);
	  }

	  template<class TVariableType> inline
	  void NodeAddDofwithReaction(Node<3>& rThisNode, TVariableType const& rThisVariable,  TVariableType const& rReactionVariable)
	  {
// 	      KRATOS_WATCH(rReactionVariable);
		rThisNode.pAddDof(rThisVariable, rReactionVariable);
	  }
	  
	  
	  void PointSetX0(Node<3>& ThisPoint, double Value)
	  {
	    ThisPoint.X0() = Value;
	  }
	  
	  void PointSetY0(Node<3>& ThisPoint, double Value)
	  {
	    ThisPoint.Y0() = Value;
	  }
	  
	  void PointSetZ0(Node<3>& ThisPoint, double Value)
	  {
	    ThisPoint.Z0() = Value;
	  }
	  
	  double PointGetX0(Node<3>& ThisPoint)
	  {
	    return ThisPoint.X0();
	  }
	  
	  double PointGetY0(Node<3>& ThisPoint)
	  {
	    return ThisPoint.Y0();
	  }
	  
	  double PointGetZ0(Node<3>& ThisPoint)
	  {
	    return ThisPoint.Z0();
	  }





		void  AddNodeToPython()
		{
			typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;

			typedef MeshType::NodeType NodeType;


			class_<IndexedObject>("IndexedObject")
				.add_property("Id", &IndexedObject::GetId, &IndexedObject::SetId)
				.def(self_ns::str(self))
				;


			class_<Dof<double>, Dof<double>::Pointer, bases<IndexedObject> >("Dof", no_init)
			  ;
			//class_<Dof, Dof::Pointer>("Dof", init<int, const Dof::VariableType&,  optional<const Dof::VariableType&, const Dof::VariableType&, const Dof::VariableType&> >())
			//.def("GetVariable", &Dof::GetVariable, return_internal_reference<>())
			//.def("GetReaction", &Dof::GetReaction, return_internal_reference<>())
			//.def("GetTimeDerivative", &Dof::GetTimeDerivative, return_internal_reference<>())
			//.def("GetSecondTimeDerivative", &Dof::GetSecondTimeDerivative, return_internal_reference<>())
			//.def("NodeIndex", &Dof::NodeIndex)
			//.add_property("EquationId", &Dof::EquationId, &Dof::SetEquationId)
			//.def("Fix", &Dof::FixDof)
			//.def("Free", &Dof::FreeDof)
			//.def("IsFixed", &Dof::IsFixed)
			//.def("HasTimeDerivative", &Dof::HasTimeDerivative)
			//.def("HasSecondTimeDerivative", &Dof::HasSecondTimeDerivative)
			//.def(self_ns::str(self))
			//      ;

//			void (NodeType::*pointer_to_double_variable_fix)(Variable<double> const&) = &NodeType::Fix;
//			void (NodeType::*pointer_to_double_component_fix)(VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > const&) = &NodeType::Fix;


			class_<NodeType, NodeType::Pointer, bases<NodeType::BaseType, IndexedObject > >("Node", init<int, double, double, double>())
				.def(init<int, const Point<3>& >())
                                .def(VariableIndexingPython<NodeType, Variable<bool> >())
				.def(VariableIndexingPython<NodeType, Variable<bool> >())
				.def(VariableIndexingPython<NodeType, Variable<int> >())
				.def(VariableIndexingPython<NodeType, Variable<double> >())
				.def(VariableIndexingPython<NodeType, Variable<array_1d<double, 3> > >())
				.def(VariableIndexingPython<NodeType, Variable<vector<double> > >())
				.def(VariableIndexingPython<NodeType, Variable<matrix<double> > >())
				.def(VariableIndexingPython<NodeType, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >())
                                .def(SolutionStepVariableIndexingPython<NodeType, Variable<bool> >()) 
				.def(SolutionStepVariableIndexingPython<NodeType, Variable<int> >())
				.def(SolutionStepVariableIndexingPython<NodeType, Variable<double> >())
				.def(SolutionStepVariableIndexingPython<NodeType, Variable<array_1d<double, 3> > >())
				.def(SolutionStepVariableIndexingPython<NodeType, Variable<vector<double> > >())
				.def(SolutionStepVariableIndexingPython<NodeType, Variable<matrix<double> > >())
				.def(SolutionStepVariableIndexingPython<NodeType, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >())
				.def("GetBufferSize", &NodeType::GetBufferSize)
 				//.def("AddDof", &NodeType::pAddDof, NodeType_padd_dof_overloads())
 				.def("AddDof", NodeAddDof<Variable<double> >)
 				.def("AddDof", NodeAddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
 				.def("AddDof", NodeAddDofwithReaction<Variable<double> >)
 				.def("AddDof", NodeAddDofwithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
				.def("Fix", NodeFix<Variable<double> >)
 				.def("Fix", NodeFix<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
 				.def("Free", NodeFree<Variable<double> >)
 				.def("Free", NodeFree<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
 				.def("IsFixed", NodeIsFixed<Variable<double> >)
 				.def("IsFixed", NodeIsFixed<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
 				.def("HasDofFor", NodeHasDofFor<Variable<double> >)
 				.def("HasDofFor", NodeHasDofFor<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
// 				.def("IsFixed", &NodeType::IsFixed)
// 				.def("HasDofFor", &NodeType::HasDofFor)
				.def(self_ns::str(self))
                                .def("OverwriteSolutionStepData", &NodeType::OverwriteSolutionStepData)
				.add_property("X0", PointGetX0, PointSetX0)
				.add_property("Y0", PointGetY0, PointSetY0)
				.add_property("Z0", PointGetZ0, PointSetZ0)
				;

			
			PointerVectorSetPythonInterface<MeshType::NodesContainerType>::CreateInterface("NodesArray")
				;

		}

	}  // namespace Python.

} // Namespace Kratos

