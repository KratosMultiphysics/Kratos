// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/mesh.h"
#include "includes/properties.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "python/add_mesh_to_python.h"
#include "python/pointer_vector_set_python_interface.h"
//#include "python/variable_indexing_python.h"
//#include "python/solution_step_variable_indexing_python.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

template< class TContainerType, class TVariableType >
void SetValueHelperFunction(TContainerType& el, const TVariableType& rVar,const typename TVariableType::Type& Data)
{
    el.SetValue(rVar,Data);
}

template< class TContainerType, class TVariableType >
typename TVariableType::Type GetValueHelperFunction(TContainerType& el, const TVariableType& rVar)
{
    return el.GetValue(rVar);
}

typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
typedef MeshType::NodeType NodeType;
typedef MeshType::NodesContainerType NodesContainerType;
typedef Geometry<Node<3> >::PointsArrayType NodesArrayType;
typedef Geometry<Node<3> >::IntegrationPointsArrayType IntegrationPointsArrayType;

double GetAreaFromCondition( Condition& dummy )
{
    return( dummy.GetGeometry().Area() );
}

double GetAreaFromElement( Element& dummy )
{
    return( dummy.GetGeometry().Area() );
}

Properties::Pointer GetPropertiesFromElement( Element& pelem )
{
    return( pelem.pGetProperties() );
}
void SetPropertiesFromElement( Element& pelem, Properties::Pointer pProperties )
{
     pelem.pGetProperties() =pProperties ;
}

Properties::Pointer GetPropertiesFromCondition( Condition& pcond )
{
    return( pcond.pGetProperties() );
}
void SetPropertiesFromCondition( Condition& pcond, Properties::Pointer pProperties )
{
     pcond.pGetProperties() =pProperties ;
}

NodeType::Pointer GetNodeFromElement( Element& dummy, unsigned int index )
{
    return( dummy.GetGeometry().pGetPoint(index) );
}

boost::python::list GetNodesFromElement( Element& dummy )
{
    boost::python::list nodes_list;
    for( unsigned int i=0; i<dummy.GetGeometry().size(); i++ )
    {
        nodes_list.append( dummy.GetGeometry().pGetPoint(i) );
    }
    return( nodes_list );
}

NodeType::Pointer GetNodeFromCondition( Condition& dummy, unsigned int index )
{
    return( dummy.GetGeometry().pGetPoint(index) );
}

boost::python::list GetNodesFromCondition( Condition& dummy )
{
    boost::python::list nodes_list;
    for( unsigned int i=0; i<dummy.GetGeometry().size(); i++ )
    {
        nodes_list.append( dummy.GetGeometry().pGetPoint(i) );
    }
    return( nodes_list );
}

boost::python::list GetIntegrationPointsFromElement( Element& dummy )
{
    boost::python::list integration_points_list;
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    for( unsigned int i=0; i< integration_points.size(); i++ )
    {
        boost::python::list item;
        Point<3> pnt;
        dummy.GetGeometry().GlobalCoordinates(pnt, integration_points[i]);
        for( unsigned int j=0; j<3; j++ )
            item.append( pnt[j] );
        integration_points_list.append( item );
    }
    return( integration_points_list );
}

boost::python::list CalculateOnIntegrationPointsVector( ModelPart& rModelPart,
        Element& dummy, const Variable<Vector>& rVariable )
{
    std::vector<Vector> Output;
    dummy.CalculateOnIntegrationPoints( rVariable, Output,
                                        rModelPart.GetProcessInfo() );
    boost::python::list result;
    for( unsigned int j=0; j<Output.size(); j++ )
        result.append( Output[j] );
    return result;
}

template< class TObject >
boost::python::list GetValuesOnIntegrationPointsDouble( TObject& dummy,
        const Variable<double>& rVariable, const ProcessInfo& rCurrentProcessInfo )
{
    boost::python::list values_list;
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    std::vector<double> values( integration_points.size() );
    dummy.CalculateOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
    for( unsigned int i=0; i<values.size(); i++ )
    {
        boost::python::list integration_point_value;
        integration_point_value.append( values[i] );
        values_list.append( integration_point_value );
    }
    return( values_list );
}

template< class TObject >
void SetValuesOnIntegrationPointsDouble( TObject& dummy, const Variable<double>& rVariable, boost::python::list values_list,  const ProcessInfo& rCurrentProcessInfo )
{
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    std::vector<double> values( integration_points.size() );
    for( unsigned int i=0; i<integration_points.size(); i++ )
    {
        boost::python::extract<double> x( values_list[i] );
        if( x.check() )
        {
            values[i] = x();
        }
        else
            break;
    }
    dummy.SetValueOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
}


template< class TObject >
boost::python::list GetValuesOnIntegrationPointsArray1d( TObject& dummy,
        const Variable<array_1d<double,3> >& rVariable, const ProcessInfo& rCurrentProcessInfo )
{
    boost::python::list values_list;
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    std::vector<array_1d<double,3> > values( integration_points.size() );
    dummy.CalculateOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
    for( unsigned int i=0; i<values.size(); i++ )
    {
        boost::python::list integration_point_value;
        for( int j=0; j<3; j++ )
            integration_point_value.append( values[i][j] );
        values_list.append( integration_point_value );
    }
    return( values_list );
}

template< class TObject >
void SetValuesOnIntegrationPointsArray1d( TObject& dummy, const Variable< array_1d<double,3> >& rVariable, boost::python::list values_list,  const ProcessInfo& rCurrentProcessInfo )
{
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    std::vector< array_1d<double,3> > values( integration_points.size() );
    for( unsigned int i=0; i<integration_points.size(); i++ )
    {
        boost::python::extract< array_1d<double,3> > x( values_list[i] );
        if( x.check() )
        {
            values[i] = x();
        }
        else
            break;
    }
    dummy.SetValueOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
}

template< class TObject >
boost::python::list GetValuesOnIntegrationPointsVector( TObject& dummy,
        const Variable<Vector>& rVariable, const ProcessInfo& rCurrentProcessInfo )
{
    boost::python::list values_list;
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    std::vector<Vector> values( integration_points.size() );
    dummy.CalculateOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
    for( unsigned int i=0; i<values.size(); i++ )
    {
        boost::python::list integration_point_value;
        for( unsigned int j=0; j<values[i].size(); j++ )
            integration_point_value.append( values[i][j] );
        values_list.append( integration_point_value );
    }
    return( values_list );
}

template< class TObject >
void SetValuesOnIntegrationPointsVector( TObject& dummy,
        const Variable<Vector>& rVariable, boost::python::list values_list, unsigned int len_values_list_item, const ProcessInfo& rCurrentProcessInfo )
{
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    std::vector<Vector> values( integration_points.size() );
    for( unsigned int i=0; i<integration_points.size(); i++ )
    {
        Vector value_item = ZeroVector(len_values_list_item);
        for( unsigned int j=0; j<len_values_list_item; j++ )
        {
            boost::python::extract<double> x( values_list[i][j] );
            if( x.check() )
            {
                value_item[j] = x();
            }
            else
                break;
        }
        values[i] = value_item;
    }
    dummy.SetValueOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
}


template< class TObject >
boost::python::list GetValuesOnIntegrationPointsMatrix( TObject& dummy,
        const Variable<Matrix>& rVariable, const ProcessInfo& rCurrentProcessInfo )
{
    boost::python::list values_list;
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    std::vector<Matrix> values( integration_points.size() );
    dummy.CalculateOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
    for( unsigned int i=0; i<values.size(); i++ )
    {
        boost::python::list integration_point_value;
        for( unsigned int j=0; j<values[i].size1(); j++ )
            for( unsigned int k=0; k<values[i].size2(); k++ )
                integration_point_value.append( values[i](j,k) );
        values_list.append( integration_point_value );
    }
    return( values_list );
}

template< class TDataType >
TDataType ElementCalculateInterface(Element& dummy, Variable<TDataType>& rVariable, ProcessInfo& rCurrentProcessInfo)
{
    TDataType aux;
    dummy.Calculate(rVariable, aux, rCurrentProcessInfo);
    return aux;
}

void SetValuesOnIntegrationPointsConstitutiveLaw( Element& dummy, const Variable<ConstitutiveLaw::Pointer>& rVariable, boost::python::list values_list, const ProcessInfo& rCurrentProcessInfo )
{
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    std::vector<ConstitutiveLaw::Pointer> values( integration_points.size() );
    for( unsigned int i=0; i<integration_points.size(); i++ )
    {
        ConstitutiveLaw::Pointer value_item;
        boost::python::extract<ConstitutiveLaw::Pointer> x( values_list[i] );
        values[i] = x();
    }
    dummy.SetValueOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
}


void  AddMeshToPython()
{
//             typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
//             typedef MeshType::NodeType NodeType;

    //     class_<Dof, Dof::Pointer>("Dof", init<int, const Dof::VariableType&,  optional<const Dof::VariableType&, const Dof::VariableType&, const Dof::VariableType&> >())
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
    
    class_<GeometricalObject, GeometricalObject::Pointer, bases<GeometricalObject::BaseType, Flags > >("GeometricalObject", init<int>())
    ;

    class_<Element, Element::Pointer, bases<Element::BaseType, Flags > >("Element", init<int>())
    //.def(init<int, const Point<3>& >())
    .add_property("Properties", GetPropertiesFromElement, SetPropertiesFromElement)
    .def("__setitem__", SetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< Vector > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< Vector > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< Vector > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< vector<int> > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< vector<int> > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< vector<int> > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< vector<int> > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< Matrix > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< Matrix > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< Matrix > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< Matrix > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< int > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< int > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< int > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< int > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< double > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< double > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< double > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< double > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< bool > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< bool > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< bool > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< bool > >)


    .def("GetArea", GetAreaFromElement )
    .def("GetNode", GetNodeFromElement )
    .def("GetNodes", GetNodesFromElement )
    .def("GetIntegrationPoints", GetIntegrationPointsFromElement )
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPointsVector)
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsDouble<Element>)
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsArray1d<Element>)
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsVector<Element>)
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsMatrix<Element>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsVector<Element>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsConstitutiveLaw)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsDouble<Element>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Element>)
    .def("ResetConstitutiveLaw", &Element::ResetConstitutiveLaw)
    
    .def("Calculate", &ElementCalculateInterface<double>)
    .def("Calculate", &ElementCalculateInterface<array_1d<double,3> >)
    .def("Calculate", &ElementCalculateInterface<Vector >)
    .def("Calculate", &ElementCalculateInterface<Matrix >)
    
    
    
    

    //.def("__setitem__", SetValueHelperFunction< Element, Variable< VectorComponentAdaptor< array_1d<double, 3>  > > >)
    //.def("__getitem__", GetValueHelperFunction< Element, Variable< VectorComponentAdaptor< array_1d<double, 3>  > > >)
    //.def("SetValue", SetValueHelperFunction< Element, Variable< VectorComponentAdaptor< array_1d<double, 3>  > > >)
    //.def("GetValue", GetValueHelperFunction< Element, Variable< VectorComponentAdaptor< array_1d<double, 3>  > > >)

    /*                  .def(VariableIndexingPython<Element, Variable<int> >())
                        .def(VariableIndexingPython<Element, Variable<double> >())
                        .def(VariableIndexingPython<Element, Variable<array_1d<double, 3> > >())
                        .def(VariableIndexingPython<Element, Variable< Vector > >())
                        .def(VariableIndexingPython<Element, Variable< Matrix > >())
                        .def(VariableIndexingPython<Element, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >())
                        .def(SolutionStepVariableIndexingPython<Element, Variable<int> >())
                        .def(SolutionStepVariableIndexingPython<Element, Variable<double> >())
                        .def(SolutionStepVariableIndexingPython<Element, Variable<array_1d<double, 3> > >())
                        .def(SolutionStepVariableIndexingPython<Element, Variable<vector<double> > >())
                        .def(SolutionStepVariableIndexingPython<Element, Variable<matrix<double> > >())
                        .def(SolutionStepVariableIndexingPython<Element, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >()) */
//                     .def("Initialize", &Element::Initialize)
    //.def("CalculateLocalSystem", &Element::CalculateLocalSystem)
    .def(self_ns::str(self))
    ;

    PointerVectorSetPythonInterface<MeshType::ElementsContainerType>::CreateInterface("ElementsArray")
    ;

    class_<Condition, Condition::Pointer, bases<Condition::BaseType, Flags > >("Condition", init<int>())
    //.def(init<int, const Point<3>& >())
    .add_property("Properties", GetPropertiesFromCondition, SetPropertiesFromCondition)
    .def("__setitem__", SetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< Vector > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< Vector > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< Vector > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< vector<int> > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< vector<int> > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< vector<int> > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< vector<int> > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< Matrix > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< Matrix > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< Matrix > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< Matrix > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< int > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< int > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< int > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< int > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< double > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< double > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< double > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< double > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< bool > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< bool > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< bool > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< bool > >)

    .def("GetNode", GetNodeFromCondition )
    .def("GetNodes", GetNodesFromCondition )

    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsDouble<Condition>)
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsArray1d<Condition>)
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsVector<Condition>)
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsMatrix<Condition>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsVector<Condition>)
    //.def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsConstitutiveLaw)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsDouble<Condition>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Condition>)
	.def("GetArea",GetAreaFromCondition)



//				.def(VariableIndexingPython<Condition, Variable<int> >())
//				.def(VariableIndexingPython<Condition, Variable<double> >())
//				.def(VariableIndexingPython<Condition, Variable<array_1d<double, 3> > >())
//				.def(VariableIndexingPython<Condition, Variable< Vector > >())
//				.def(VariableIndexingPython<Condition, Variable< Matrix > >())
//				.def(VariableIndexingPython<Condition, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >())
    /*				.def(SolutionStepVariableIndexingPython<Condition, Variable<int> >())
    				.def(SolutionStepVariableIndexingPython<Condition, Variable<double> >())
    				.def(SolutionStepVariableIndexingPython<Condition, Variable<array_1d<double, 3> > >())
    				.def(SolutionStepVariableIndexingPython<Condition, Variable<vector<double> > >())
    				.def(SolutionStepVariableIndexingPython<Condition, Variable<matrix<double> > >())
    				.def(SolutionStepVariableIndexingPython<Condition, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >())
    */				//.def("Initialize", &Condition::Initialize)
    //.def("CalculateLocalSystem", &Condition::CalculateLocalSystem)
    .def("Info", &Condition::Info)
    .def(self_ns::str(self))
    ;

    PointerVectorSetPythonInterface<MeshType::ConditionsContainerType>::CreateInterface("ConditionsArray")
    ;

    class_<MeshType, MeshType::Pointer, bases<DataValueContainer, Flags>, boost::noncopyable >("Mesh")
    .add_property("Nodes", &MeshType::pNodes,&MeshType::SetNodes)
    .def("NodesArray", &MeshType::NodesArray, return_internal_reference<>())
    .add_property("Elements", &MeshType::pElements,&MeshType::SetElements)
    .def("ElementsArray", &MeshType::ElementsArray, return_internal_reference<>())
    .add_property("Conditions", &MeshType::pConditions,&MeshType::SetConditions)
    .def("ConditionsArray", &MeshType::ConditionsArray, return_internal_reference<>())
    .add_property("Properties", &MeshType::pProperties,&MeshType::SetProperties)
    .def("PropertiesArray", &MeshType::PropertiesArray, return_internal_reference<>())
	.def("HasNode", &MeshType::HasNode)
	.def("HasProperties", &MeshType::HasProperties)
	.def("HasElement", &MeshType::HasElement)
	.def("HasCondition", &MeshType::HasCondition)
    .def(self_ns::str(self))
    ;
}
}  // namespace Python.
} // Namespace Kratos
