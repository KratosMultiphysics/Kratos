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
bool HasHelperFunction(TContainerType& el, const TVariableType& rVar)
{
    return el.Has(rVar);
}

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
typedef Point::CoordinatesArrayType CoordinatesArrayType;

array_1d<double,3> GetNormalFromCondition( 
    Condition& dummy, 
    CoordinatesArrayType& LocalCoords
    )
{
    return( dummy.GetGeometry().UnitNormal(LocalCoords) );
}

array_1d<double,3> FastGetNormalFromCondition(Condition& dummy)
{
    CoordinatesArrayType LocalCoords;
    LocalCoords.clear();
    return( dummy.GetGeometry().UnitNormal(LocalCoords) );
}

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
     pelem.SetProperties(pProperties) ;
}

Properties::Pointer GetPropertiesFromCondition( Condition& pcond )
{
    return( pcond.pGetProperties() );
}
void SetPropertiesFromCondition( Condition& pcond, Properties::Pointer pProperties )
{
     pcond.SetProperties(pProperties) ;
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

void ConditionCalculateLocalSystemStandard( Condition& dummy, 
                                                Matrix& rLeftHandSideMatrix,
                                                Vector& rRightHandSideVector,
                                                ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
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
        Point point;
        dummy.GetGeometry().GlobalCoordinates(point, integration_points[i]);
        for( unsigned int j=0; j<3; j++ )
            item.append( point[j] );
        integration_points_list.append( item );
    }
    return( integration_points_list );
}

boost::python::list CalculateOnIntegrationPointsDouble(
        Element& dummy, const Variable<double>& rVariable, ProcessInfo& rProcessInfo )
{
    std::vector<double> Output;
    dummy.CalculateOnIntegrationPoints( rVariable, Output, rProcessInfo);
    boost::python::list result;
    for( unsigned int j=0; j<Output.size(); j++ )
    {
        result.append( Output[j] );
    }
    return result;
}

boost::python::list CalculateOnIntegrationPointsArray1d(
        Element& dummy, const Variable<array_1d<double, 3>>& rVariable, ProcessInfo& rProcessInfo )
{
    std::vector<array_1d<double, 3>> Output;
    dummy.CalculateOnIntegrationPoints( rVariable, Output, rProcessInfo);
    boost::python::list result;
    for( unsigned int j=0; j<Output.size(); j++ )
    {
        result.append( Output[j] );
    }
    return result;
}

boost::python::list CalculateOnIntegrationPointsVector(
        Element& dummy, const Variable<Vector>& rVariable, ProcessInfo& rProcessInfo )
{
    std::vector<Vector> Output;
    dummy.CalculateOnIntegrationPoints( rVariable, Output, rProcessInfo);
    boost::python::list result;
    for( unsigned int j=0; j<Output.size(); j++ )
    {
        result.append( Output[j] );
    }
    return result;
}

boost::python::list CalculateOnIntegrationPointsMatrix(
        Element& dummy, const Variable<Matrix>& rVariable, ProcessInfo& rProcessInfo )
{
    std::vector<Matrix> Output;
    dummy.CalculateOnIntegrationPoints( rVariable, Output,rProcessInfo );
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

void ElementCalculateLocalSystem1(Element& dummy,
        Matrix& rLeftHandSideMatrix,
        Vector& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
}

template<class TDataType>
void ElementCalculateSensitivityMatrix(Element& dummy,
        const Variable<TDataType>& rDesignVariable,
        Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateSensitivityMatrix(rDesignVariable,rOutput,rCurrentProcessInfo);
}

void ElementGetFirstDerivativesVector1(Element& dummy,
        Vector& rOutput)
{
    dummy.GetFirstDerivativesVector(rOutput,0);
}

void ElementGetFirstDerivativesVector2(Element& dummy,
        Vector& rOutput,
        int step)
{
    dummy.GetFirstDerivativesVector(rOutput,step);
}

void ElementGetSecondDerivativesVector1(Element& dummy,
        Vector& rOutput)
{
    dummy.GetSecondDerivativesVector(rOutput,0);
}

void ElementGetSecondDerivativesVector2(Element& dummy,
        Vector& rOutput,
        int step)
{
    dummy.GetSecondDerivativesVector(rOutput,step);
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
    .add_property("Properties", GetPropertiesFromElement, SetPropertiesFromElement)
    .def("__setitem__", SetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("Has", HasHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< Vector > >)
    .def("Has", HasHelperFunction< Element, Variable< Vector > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< Vector > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< Vector > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< vector<int> > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< vector<int> > >)
    .def("Has", HasHelperFunction< Element, Variable< vector<int> > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< vector<int> > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< vector<int> > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< Matrix > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< Matrix > >)
    .def("Has", HasHelperFunction< Element, Variable< Matrix > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< Matrix > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< Matrix > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< int > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< int > >)
    .def("Has", HasHelperFunction< Element, Variable< int > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< int > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< int > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< double > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< double > >)
    .def("Has", HasHelperFunction< Element, Variable< double > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< double > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< double > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< bool > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< bool > >)
    .def("Has", HasHelperFunction< Element, Variable< bool > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< bool > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< bool > >)


    .def("GetArea", GetAreaFromElement )
    .def("GetNode", GetNodeFromElement )
    .def("GetNodes", GetNodesFromElement )
    .def("GetIntegrationPoints", GetIntegrationPointsFromElement )
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPointsDouble)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPointsArray1d)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPointsVector)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPointsMatrix)
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
    .def("CalculateMassMatrix", &Element::CalculateMassMatrix)
    .def("CalculateDampingMatrix", &Element::CalculateDampingMatrix)
    .def("CalculateLocalSystem", &ElementCalculateLocalSystem1)
    .def("CalculateFirstDerivativesLHS", &Element::CalculateFirstDerivativesLHS)
    .def("CalculateSecondDerivativesLHS", &Element::CalculateSecondDerivativesLHS)
    .def("CalculateLocalVelocityContribution", &Element::CalculateLocalVelocityContribution)
    .def("GetFirstDerivativesVector", &ElementGetFirstDerivativesVector1)
    .def("GetFirstDerivativesVector", &ElementGetFirstDerivativesVector2)
    .def("GetSecondDerivativesVector", &ElementGetSecondDerivativesVector1)
    .def("GetSecondDerivativesVector", &ElementGetSecondDerivativesVector2)
    .def("CalculateSensitivityMatrix", &ElementCalculateSensitivityMatrix<double>)
    .def("CalculateSensitivityMatrix", &ElementCalculateSensitivityMatrix<array_1d<double,3> >)
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
    .def("Initialize", &Element::Initialize)
    //.def("CalculateLocalSystem", &Element::CalculateLocalSystem)
    .def(self_ns::str(self))
    ;

    PointerVectorSetPythonInterface<MeshType::ElementsContainerType>::CreateInterface("ElementsArray")
    ;

    class_<Condition, Condition::Pointer, bases<Condition::BaseType, Flags > >("Condition", init<int>())
    .add_property("Properties", GetPropertiesFromCondition, SetPropertiesFromCondition)
    .def("__setitem__", SetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("Has", HasHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< Vector > >)
    .def("Has", HasHelperFunction< Condition, Variable< Vector > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< Vector > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< Vector > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< vector<int> > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< vector<int> > >)
    .def("Has", HasHelperFunction< Condition, Variable< vector<int> > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< vector<int> > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< vector<int> > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< Matrix > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< Matrix > >)
    .def("Has", HasHelperFunction< Condition, Variable< Matrix > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< Matrix > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< Matrix > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< int > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< int > >)
    .def("Has", HasHelperFunction< Condition, Variable< int > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< int > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< int > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< double > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< double > >)
    .def("Has", HasHelperFunction< Condition, Variable< double > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< double > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< double > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< bool > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< bool > >)
    .def("Has", HasHelperFunction< Condition, Variable< bool > >)
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
    .def("GetNormal",GetNormalFromCondition)
    .def("GetNormal",FastGetNormalFromCondition)
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
    */
    .def("Initialize", &Condition::Initialize)
    .def("CalculateLocalSystem", &ConditionCalculateLocalSystemStandard)
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
