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
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/mesh.h"
#include "includes/properties.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "python/containers_interface.h"
#include "python/add_mesh_to_python.h"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

typedef Node<3> NodeType;

template<class TVariableType> inline
void NodeFix(NodeType& rThisNode, TVariableType const& rThisVariable)
{
    rThisNode.Fix(rThisVariable);
}

template<class TVariableType> inline
void NodeFree(NodeType& rThisNode, TVariableType const& rThisVariable)
{
    rThisNode.Free(rThisVariable);
}

template<class TVariableType> inline
bool NodeIsFixed(NodeType& rThisNode, TVariableType const& rThisVariable)
{
    return rThisNode.IsFixed(rThisVariable);
}

template<class TVariableType> inline
bool NodeHasDofFor(NodeType& rThisNode, TVariableType const& rThisVariable)
{
    return rThisNode.HasDofFor(rThisVariable);
}

template<class TVariableType> inline
void NodeAddDof(NodeType& rThisNode, TVariableType const& rThisVariable)
{
    rThisNode.pAddDof(rThisVariable);
}

template<class TVariableType> inline
void NodeAddDofwithReaction(NodeType& rThisNode, TVariableType const& rThisVariable,  TVariableType const& rReactionVariable)
{
// 	      KRATOS_WATCH(rReactionVariable);
    rThisNode.pAddDof(rThisVariable, rReactionVariable);
}

template<class TVariableType>
bool NodeSolutionStepsDataHas(NodeType& rThisNode, const TVariableType& rThisVariable)
{
    return rThisNode.SolutionStepsDataHas(rThisVariable);
}

void PointSetX0(NodeType& ThisPoint, double Value)
{
    ThisPoint.X0() = Value;
}

void PointSetY0(NodeType& ThisPoint, double Value)
{
    ThisPoint.Y0() = Value;
}

void PointSetZ0(NodeType& ThisPoint, double Value)
{
    ThisPoint.Z0() = Value;
}

double PointGetX0(NodeType& ThisPoint)
{
    return ThisPoint.X0();
}

double PointGetY0(NodeType& ThisPoint)
{
    return ThisPoint.Y0();
}

double PointGetZ0(NodeType& ThisPoint)
{
    return ThisPoint.Z0();
}

template< class TBinderType, typename TContainerType, typename TVariableType > void IndexingUtility(TBinderType& binder)
    {
        //data value container
        binder.def("__contains__", [](const TContainerType& node, const TVariableType& rV){return node.Has(rV);} );
        binder.def("__setitem__", [](TContainerType& node, const TVariableType& rV, const typename TVariableType::Type rValue){node.SetValue(rV, rValue);} );
        binder.def("__getitem__", [](TContainerType& node, const TVariableType& rV){return node.GetValue(rV);} );
        binder.def("Has", [](const TContainerType& node, const TVariableType& rV){return node.Has(rV);} );
        binder.def("SetValue",  [](TContainerType& node, const TVariableType& rV, const typename TVariableType::Type& rValue){node.SetValue(rV, rValue);} );
        binder.def("GetValue", [](TContainerType& node, const TVariableType& rV){return node.GetValue(rV);} );

        //solution steps data value container
        binder.def("HasSolutionStepValue", [](const TContainerType& node, const TVariableType& rV){return node.Has(rV);} ); //WARNING Previously it was Has identical to the previous!!
        binder.def("SetSolutionStepValue",  [](TContainerType& node, const TVariableType& rV,const typename TVariableType::Type rValue){node.GetSolutionStepValue(rV) = rValue;} );
        binder.def("SetSolutionStepValue",  [](TContainerType& node,
                                               const TVariableType& rV,
                                               typename TContainerType::IndexType SolutionStepIndex,
                                               const typename TVariableType::Type rValue){node.GetSolutionStepValue(rV, SolutionStepIndex) = rValue;} );
        binder.def("GetSolutionStepValue", [](TContainerType& node, const TVariableType& rV){return node.GetSolutionStepValue(rV);} );
        binder.def("GetSolutionStepValue", [](TContainerType& node, const TVariableType& rV, typename TContainerType::IndexType SolutionStepIndex ){return node.GetSolutionStepValue(rV, SolutionStepIndex);} );
    }



void  AddNodeToPython(pybind11::module& m)
{
    typedef Mesh<NodeType, Properties, Element, Condition> MeshType;

    typedef MeshType::NodeType NodeType;


    class_<IndexedObject, IndexedObject::Pointer>(m,"IndexedObject")
    .def_property("Id", &IndexedObject::GetId, &IndexedObject::SetId)
    .def("__str__", PrintObject<IndexedObject>)
    ;

    class_<Dof<double>, Dof<double>::Pointer, IndexedObject >(m,"Dof")
    ;

    typedef  class_<NodeType, NodeType::Pointer, NodeType::BaseType, IndexedObject, Flags > NodeBinderType;
    NodeBinderType node_binder(m,"Node");
    node_binder.def(init<NodeType::IndexType, double, double, double>());
    node_binder.def(init<NodeType::IndexType, const Point& >());

    IndexingUtility<NodeBinderType,NodeType,Variable<bool> >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,Variable<int> >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,Variable<double> >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,Variable<array_1d<double, 3> > >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,Variable<array_1d<double, 4> > >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,Variable<array_1d<double, 6> > >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,Variable<array_1d<double, 9> > >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,Variable<Vector > >(node_binder);
    IndexingUtility<NodeBinderType,NodeType,Variable<Matrix > >(node_binder);

    node_binder.def("GetBufferSize", &NodeType::GetBufferSize);
    node_binder.def("AddDof", NodeAddDof<Variable<double> >);
    node_binder.def("AddDof", NodeAddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >);
    node_binder.def("AddDof", NodeAddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >);
    node_binder.def("AddDof", NodeAddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >);
    node_binder.def("AddDof", NodeAddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >);
    node_binder.def("AddDof", NodeAddDofwithReaction<Variable<double> >);
    node_binder.def("AddDof", NodeAddDofwithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >);
    node_binder.def("AddDof", NodeAddDofwithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >);
    node_binder.def("AddDof", NodeAddDofwithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >);
    node_binder.def("AddDof", NodeAddDofwithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >);
    node_binder.def("Fix", NodeFix<Variable<double> >);
    node_binder.def("Fix", NodeFix<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >);
    node_binder.def("Fix", NodeFix<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >);
    node_binder.def("Fix", NodeFix<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >);
    node_binder.def("Fix", NodeFix<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >);
    node_binder.def("Free", NodeFree<Variable<double> >);
    node_binder.def("Free", NodeFree<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >);
    node_binder.def("Free", NodeFree<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >);
    node_binder.def("Free", NodeFree<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >);
    node_binder.def("Free", NodeFree<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >);
    node_binder.def("IsFixed", NodeIsFixed<Variable<double> >);
    node_binder.def("IsFixed", NodeIsFixed<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >);
    node_binder.def("IsFixed", NodeIsFixed<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >);
    node_binder.def("IsFixed", NodeIsFixed<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >);
    node_binder.def("IsFixed", NodeIsFixed<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >);
    node_binder.def("HasDofFor", NodeHasDofFor<Variable<double> >);
    node_binder.def("HasDofFor", NodeHasDofFor<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >);
    node_binder.def("HasDofFor", NodeHasDofFor<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >);
    node_binder.def("HasDofFor", NodeHasDofFor<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >);
    node_binder.def("HasDofFor", NodeHasDofFor<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<bool> >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<int> >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<double> >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<array_1d<double, 3> > >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<array_1d<double, 4> > >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<array_1d<double, 6> > >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<array_1d<double, 9> > >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<Vector> >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<Variable<DenseMatrix<double> > >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >);
    node_binder.def("SolutionStepsDataHas", &NodeSolutionStepsDataHas<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >);
    node_binder.def("__str__", PrintObject<NodeType>);
    node_binder.def("OverwriteSolutionStepData", &NodeType::OverwriteSolutionStepData);
    node_binder.def_property("X0", PointGetX0, PointSetX0);
    node_binder.def_property("Y0", PointGetY0, PointSetY0);
    node_binder.def_property("Z0", PointGetZ0, PointSetZ0);

    PointerVectorSetPythonInterface<MeshType::NodesContainerType>().CreateInterface(m,"NodesArray");

}

}  // namespace Python.

} // Namespace Kratos

