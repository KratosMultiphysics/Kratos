//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
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

template<class TVariableType>
bool NodeSolutionStepsDataHas(Node<3>& rThisNode, const TVariableType& rThisVariable)
{
    return rThisNode.SolutionStepsDataHas(rThisVariable);
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


    class_<NodeType, NodeType::Pointer, bases<NodeType::BaseType, IndexedObject, Flags >, boost::noncopyable >("Node", init<int, double, double, double>())
    .def(init<int, const Point& >())
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
//    .def("SolutionStepsDataHas", &NodeType::SolutionStepsDataHas<bool>)
    .def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<bool> >)
    .def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<int> >)
    .def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<double> >)
    .def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<array_1d<double, 3> > >)
    .def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<vector<double> > >)
    .def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<matrix<double> > >)
    .def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
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

