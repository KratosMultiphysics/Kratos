//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Project includes
#include "includes/define_python.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "python/add_mesh_to_python.h"
#include "python/containers_interface.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

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
typedef Geometry<Node<3> > GeometryType;
typedef GeometryType::PointsArrayType NodesArrayType;
typedef GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
typedef Point::CoordinatesArrayType CoordinatesArrayType;

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

template <class T>
const GeometryType& GetGeometryFromObject( T& rObject )
{
    return rObject.GetGeometry();
}

NodeType::Pointer GetNodeFromElement( Element& dummy, unsigned int index )
{
    return( dummy.GetGeometry().pGetPoint(index) );
}

py::list GetNodesFromElement( Element& dummy )
{
    pybind11::list nodes_list;
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
                                                const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
}

void ConditionInitialize(Condition& dummy,
                       const ProcessInfo& rCurrentProcessInfo)
{
    dummy.Initialize(rCurrentProcessInfo);
}

py::list GetNodesFromCondition( Condition& dummy )
{
    pybind11::list nodes_list;
    for( unsigned int i=0; i<dummy.GetGeometry().size(); i++ )
    {
        nodes_list.append( dummy.GetGeometry().pGetPoint(i) );
    }
    return( nodes_list );
}

py::list GetIntegrationPointsFromElement( Element& dummy )
{
    pybind11::list integration_points_list;
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    for( unsigned int i=0; i< integration_points.size(); i++ )
    {
        pybind11::list item;
        Point point;
        dummy.GetGeometry().GlobalCoordinates(point, integration_points[i]);
        for( unsigned int j=0; j<3; j++ )
            item.append( point[j] );
        integration_points_list.append( item );
    }
    return( integration_points_list );
}

///@}
///@name Calculate on Integration Points
///@{

template< class TObject, class TDataType >
pybind11::list CalculateOnIntegrationPoints(
    TObject& dummy, const Variable<TDataType>& rVariable, const ProcessInfo& rProcessInfo)
{
    std::vector<TDataType> Output;
    dummy.CalculateOnIntegrationPoints(rVariable, Output, rProcessInfo);
    pybind11::list result;
    for (unsigned int j = 0; j < Output.size(); j++)
    {
        result.append(Output[j]);
    }
    return result;
}

///@}
///@name Get Values on Integration Points
///@{

template< class TObject >
void GetValuesOnIntegrationPoints(
    TObject& dummy,
    const Variable<Vector>& rVariable,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "GetValuesOnIntegrationPoints is deprecated. Use CalculateOnIntegrationPoints instead!" << std::endl;
}

///@}
///@name Set Values on Integration Points
///@{

template< class TObject, class TDataType >
void SetValuesOnIntegrationPoints(
    TObject& dummy, const Variable<TDataType>& rVariable, std::vector<TDataType> values, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(values.size() != dummy.GetGeometry().IntegrationPointsNumber())
        << "Sizes do not match. Size of values vector is: " << values.size() << ". The number of integration points is: "
        << dummy.GetGeometry().IntegrationPointsNumber() << std::endl;

    dummy.SetValuesOnIntegrationPoints(rVariable, values, rCurrentProcessInfo);
}

template< class TObject >
void SetValuesOnIntegrationPointsArray1d(
    TObject& dummy,
    const Variable< array_1d<double, 3> >& rVariable,
    pybind11::list values_list,
    const ProcessInfo& rCurrentProcessInfo)
{
    std::vector< array_1d<double, 3> > values(values_list.size());
    for (unsigned int i = 0; i < values_list.size(); i++)
    {
        if (py::isinstance<array_1d<double, 3>>(values_list[i])) {
            values[i] = (values_list[i]).cast<array_1d<double, 3> >();
        }
        else if (py::isinstance<pybind11::list>(values_list[i]) ||
            py::isinstance<Vector>(values_list[i]))
        {
            Vector value = (values_list[i]).cast<Vector>();
            KRATOS_ERROR_IF(value.size() != 3)
                << " parsed vector is not of size 3. Size of vector: " << value.size() << std::endl;
            values[i] = value;
        }
        else {
            KRATOS_ERROR << "expecting a list of array_1d<double,3> " << std::endl;
        }
    }
    dummy.SetValuesOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
}

template< class TObject >
void SetValuesOnIntegrationPointsVector( TObject& dummy,
        const Variable<Vector>& rVariable, pybind11::list values_list, unsigned int len_values_list_item, const ProcessInfo& rCurrentProcessInfo )
{
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    std::vector<Vector> values( integration_points.size() );
    for( unsigned int i=0; i<integration_points.size(); i++ )
    {
        if(py::isinstance<Vector>(values_list[i]))
            values[i] = (values_list[i]).cast<Vector>();
        else
            KRATOS_ERROR << "expecting a list of vectors";
    }
    dummy.SetValuesOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
}

template< class TDataType >
TDataType ElementCalculateInterface(Element& dummy, Variable<TDataType>& rVariable, const ProcessInfo& rCurrentProcessInfo)
{
    TDataType aux;
    dummy.Calculate(rVariable, aux, rCurrentProcessInfo);
    return aux;
}

void SetValuesOnIntegrationPointsConstitutiveLaw( Element& dummy, const Variable<ConstitutiveLaw::Pointer>& rVariable, pybind11::list values_list, const ProcessInfo& rCurrentProcessInfo )
{
    IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints(
                dummy.GetIntegrationMethod() );
    std::vector<ConstitutiveLaw::Pointer> values( integration_points.size() );
    for( unsigned int i=0; i<integration_points.size(); i++ )
    {
        if(py::isinstance<ConstitutiveLaw::Pointer>(values_list[i]))
            values[i] = (values_list[i]).cast<ConstitutiveLaw::Pointer>();
        else
            KRATOS_ERROR << "expecting a list of ConstitutiveLaw::Pointer";
     }
    dummy.SetValuesOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
}

void ElementCalculateLocalSystem1(Element& dummy,
        Matrix& rLeftHandSideMatrix,
        Vector& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
}

void ElementCalculateMassMatrix(Element& dummy,
                                Matrix& rMassMatrix,
                                const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);
}

void ElementCalculateDampingMatrix(Element& dummy,
                                   Matrix& rDampingMatrix,
                                   const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
}

void ElementCalculateFirstDerivativesLHS(Element& dummy,
                                         Matrix& rLeftHandSideMatrix,
                                         const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateFirstDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
}

void ElementCalculateSecondDerivativesLHS(Element& dummy,
                                          Matrix& rLeftHandSideMatrix,
                                          const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateSecondDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
}

void ElementCalculateLocalVelocityContribution(Element& dummy,
                                               Matrix& rDampingMatrix,
                                               Vector& rRightHandSideVector,
                                               const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateLocalVelocityContribution(rDampingMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

void ElementInitialize(Element& dummy,
                       const ProcessInfo& rCurrentProcessInfo)
{
    dummy.Initialize(rCurrentProcessInfo);
}

template<class TDataType>
void ElementCalculateSensitivityMatrix(Element& dummy,
        const Variable<TDataType>& rDesignVariable,
        Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateSensitivityMatrix(rDesignVariable,rOutput,rCurrentProcessInfo);
}

void ElementGetFirstDerivativesVector1(const Element& dummy,
        Vector& rOutput)
{
    dummy.GetFirstDerivativesVector(rOutput,0);
}

void ElementGetFirstDerivativesVector2(const Element& dummy,
        Vector& rOutput,
        int step)
{
    dummy.GetFirstDerivativesVector(rOutput,step);
}

void ElementGetSecondDerivativesVector1(const Element& dummy,
        Vector& rOutput)
{
    dummy.GetSecondDerivativesVector(rOutput,0);
}

void ElementGetSecondDerivativesVector2(const Element& dummy,
        Vector& rOutput,
        int step)
{
    dummy.GetSecondDerivativesVector(rOutput,step);
}

template<class TDataType>
void ConditionCalculateSensitivityMatrix(Condition& dummy,
        const Variable<TDataType>& rDesignVariable,
        Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateSensitivityMatrix(rDesignVariable,rOutput,rCurrentProcessInfo);
}

void  AddMeshToPython(pybind11::module& m)
{
//             typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
//             typedef MeshType::NodeType NodeType;

    //     py::class_<Dof, Dof::Pointer>("Dof", init<int, const Dof::VariableType&,  optional<const Dof::VariableType&, const Dof::VariableType&, const Dof::VariableType&> >())
    //.def("GetVariable", &Dof::GetVariable, py::return_value_policy::reference_internal)
    //.def("GetReaction", &Dof::GetReaction, py::return_value_policy::reference_internal)
    //.def("GetTimeDerivative", &Dof::GetTimeDerivative, py::return_value_policy::reference_internal)
    //.def("GetSecondTimeDerivative", &Dof::GetSecondTimeDerivative, py::return_value_policy::reference_internal)
    //.def("NodeIndex", &Dof::NodeIndex)
    //.def_property("EquationId", &Dof::EquationId, &Dof::SetEquationId)
    //.def("Fix", &Dof::FixDof)
    //.def("Free", &Dof::FreeDof)
    //.def("IsFixed", &Dof::IsFixed)
    //.def("HasTimeDerivative", &Dof::HasTimeDerivative)
    //.def("HasSecondTimeDerivative", &Dof::HasSecondTimeDerivative)
    //.def(self_ns::str(self))
    //      ;

    py::class_<GeometricalObject, GeometricalObject::Pointer, IndexedObject, Flags>(m,"GeometricalObject")
    .def(py::init<Kratos::GeometricalObject::IndexType>())
    ;

    py::class_<Element, Element::Pointer, Element::BaseType>(m,"Element")
    .def(py::init<Kratos::Element::IndexType>())
    .def_property("Properties", GetPropertiesFromElement, SetPropertiesFromElement)
    .def("GetGeometry", GetGeometryFromObject<Element>, py::return_value_policy::reference_internal)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("Has", HasHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< array_1d<double, 4>  > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< array_1d<double, 4>  > >)
    .def("Has", HasHelperFunction< Element, Variable< array_1d<double, 4>  > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< array_1d<double, 4>  > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< array_1d<double, 4>  > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< array_1d<double, 6>  > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< array_1d<double, 6>  > >)
    .def("Has", HasHelperFunction< Element, Variable< array_1d<double, 6>  > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< array_1d<double, 6>  > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< array_1d<double, 6>  > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< array_1d<double, 9>  > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< array_1d<double, 9>  > >)
    .def("Has", HasHelperFunction< Element, Variable< array_1d<double, 9>  > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< array_1d<double, 9>  > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< array_1d<double, 9>  > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< Vector > >)
    .def("Has", HasHelperFunction< Element, Variable< Vector > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< Vector > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< Vector > >)

    .def("__setitem__", SetValueHelperFunction< Element, Variable< DenseVector<int> > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< DenseVector<int> > >)
    .def("Has", HasHelperFunction< Element, Variable< DenseVector<int> > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< DenseVector<int> > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< DenseVector<int> > >)

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

    .def("__setitem__", SetValueHelperFunction< Element, Variable< std::string > >)
    .def("__getitem__", GetValueHelperFunction< Element, Variable< std::string > >)
    .def("Has", HasHelperFunction< Element, Variable< std::string > >)
    .def("SetValue", SetValueHelperFunction< Element, Variable< std::string > >)
    .def("GetValue", GetValueHelperFunction< Element, Variable< std::string > >)

    .def("GetNode", GetNodeFromElement )
    .def("GetNodes", GetNodesFromElement )
    .def("GetIntegrationPoints", GetIntegrationPointsFromElement )
    // CalculateOnIntegrationPoints
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, int>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, double>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, array_1d<double, 3>>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, Vector>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, Matrix>)
    // GetValuesOnIntegrationPoints
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPoints<Element>)
    // SetValuesOnIntegrationPoints
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPoints<Element, int>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsVector<Element>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsConstitutiveLaw)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPoints<Element, double>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Element>)
    .def("ResetConstitutiveLaw", &Element::ResetConstitutiveLaw)
    .def("Calculate", &ElementCalculateInterface<double>)
    .def("Calculate", &ElementCalculateInterface<array_1d<double,3> >)
    .def("Calculate", &ElementCalculateInterface<Vector >)
    .def("Calculate", &ElementCalculateInterface<Matrix >)
    .def("CalculateMassMatrix", &ElementCalculateMassMatrix)
    .def("CalculateDampingMatrix", &ElementCalculateDampingMatrix)
    .def("CalculateLocalSystem", &ElementCalculateLocalSystem1)
    .def("CalculateFirstDerivativesLHS", &ElementCalculateFirstDerivativesLHS)
    .def("CalculateSecondDerivativesLHS", &ElementCalculateSecondDerivativesLHS)
    .def("CalculateLocalVelocityContribution", &ElementCalculateLocalVelocityContribution)
    .def("GetFirstDerivativesVector", &ElementGetFirstDerivativesVector1)
    .def("GetFirstDerivativesVector", &ElementGetFirstDerivativesVector2)
    .def("GetSecondDerivativesVector", &ElementGetSecondDerivativesVector1)
    .def("GetSecondDerivativesVector", &ElementGetSecondDerivativesVector2)
    .def("CalculateSensitivityMatrix", &ElementCalculateSensitivityMatrix<double>)
    .def("CalculateSensitivityMatrix", &ElementCalculateSensitivityMatrix<array_1d<double,3> >)

//     .def(VariableIndexingPython<Element, Variable<int> >())
//     .def(VariableIndexingPython<Element, Variable<double> >())
//     .def(VariableIndexingPython<Element, Variable<array_1d<double, 3> > >())
//     .def(VariableIndexingPython<Element, Variable< Vector > >())
//     .def(VariableIndexingPython<Element, Variable< Matrix > >())
//     .def(SolutionStepVariableIndexingPython<Element, Variable<int> >())
//     .def(SolutionStepVariableIndexingPython<Element, Variable<double> >())
//     .def(SolutionStepVariableIndexingPython<Element, Variable<array_1d<double, 3> > >())
//     .def(SolutionStepVariableIndexingPython<Element, Variable<vector<double> > >())
//     .def(SolutionStepVariableIndexingPython<Element, Variable<DenseMatrix<double> > >())
    .def("Initialize", &ElementInitialize)
    //.def("CalculateLocalSystem", &Element::CalculateLocalSystem)
    .def("__str__", PrintObject<Element>)
    ;

    PointerVectorSetPythonInterface<MeshType::ElementsContainerType>().CreateInterface(m,"ElementsArray")
    ;

    py::class_<Condition, Condition::Pointer, Condition::BaseType>(m,"Condition")
    .def(py::init<Kratos::Condition::IndexType>())
    .def_property("Properties", GetPropertiesFromCondition, SetPropertiesFromCondition)
    .def("GetGeometry", GetGeometryFromObject<Condition>, py::return_value_policy::reference_internal)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("Has", HasHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< array_1d<double, 4>  > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< array_1d<double, 4>  > >)
    .def("Has", HasHelperFunction< Condition, Variable< array_1d<double, 4>  > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< array_1d<double, 4>  > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< array_1d<double, 4>  > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< array_1d<double, 6>  > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< array_1d<double, 6>  > >)
    .def("Has", HasHelperFunction< Condition, Variable< array_1d<double, 6>  > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< array_1d<double, 6>  > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< array_1d<double, 6>  > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< array_1d<double, 9>  > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< array_1d<double, 9>  > >)
    .def("Has", HasHelperFunction< Condition, Variable< array_1d<double, 9>  > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< array_1d<double, 9>  > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< array_1d<double, 9>  > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< Vector > >)
    .def("Has", HasHelperFunction< Condition, Variable< Vector > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< Vector > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< Vector > >)

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< DenseVector<int> > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< DenseVector<int> > >)
    .def("Has", HasHelperFunction< Condition, Variable< DenseVector<int> > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< DenseVector<int> > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< DenseVector<int> > >)

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

    .def("__setitem__", SetValueHelperFunction< Condition, Variable< std::string > >)
    .def("__getitem__", GetValueHelperFunction< Condition, Variable< std::string > >)
    .def("Has", HasHelperFunction< Condition, Variable< std::string > >)
    .def("SetValue", SetValueHelperFunction< Condition, Variable< std::string > >)
    .def("GetValue", GetValueHelperFunction< Condition, Variable< std::string > >)

    .def("GetNode", GetNodeFromCondition )
    .def("GetNodes", GetNodesFromCondition )

    // CalculateOnIntegrationPoints
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, int>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, double>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, array_1d<double, 3>>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, Vector>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, Matrix>)
    // GetValuesOnIntegrationPoints
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPoints<Condition>)
    // SetValuesOnIntegrationPoints
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPoints<Condition, int>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPoints<Condition, double>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsVector<Condition>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Condition>)
    //.def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsConstitutiveLaw)
    .def("CalculateSensitivityMatrix", &ConditionCalculateSensitivityMatrix<double>)
    .def("CalculateSensitivityMatrix", &ConditionCalculateSensitivityMatrix<array_1d<double,3> >)

//     .def(VariableIndexingPython<Condition, Variable<int> >())
//     .def(VariableIndexingPython<Condition, Variable<double> >())
//     .def(VariableIndexingPython<Condition, Variable<array_1d<double, 3> > >())
//     .def(VariableIndexingPython<Condition, Variable< Vector > >())
//     .def(VariableIndexingPython<Condition, Variable< Matrix > >())
//     .def(SolutionStepVariableIndexingPython<Condition, Variable<int> >())
//     .def(SolutionStepVariableIndexingPython<Condition, Variable<double> >())
//     .def(SolutionStepVariableIndexingPython<Condition, Variable<array_1d<double, 3> > >())
//     .def(SolutionStepVariableIndexingPython<Condition, Variable<vector<double> > >())
//     .def(SolutionStepVariableIndexingPython<Condition, Variable<DenseMatrix<double> > >())


    .def("Initialize", &ConditionInitialize)
    .def("CalculateLocalSystem", &ConditionCalculateLocalSystemStandard)
    .def("Info", &Condition::Info)
    .def("__str__", PrintObject<Condition>)
    ;

    PointerVectorSetPythonInterface<MeshType::ConditionsContainerType>().CreateInterface(m,"ConditionsArray")
    ;

    py::class_<MeshType, MeshType::Pointer, DataValueContainer, Flags >(m,"Mesh")
    .def_property("Nodes", &MeshType::pNodes,&MeshType::SetNodes)
    .def("NodesArray", &MeshType::NodesArray, py::return_value_policy::reference_internal)
    .def("NumberOfNodes", &MeshType::NumberOfNodes)

    .def_property("Elements", &MeshType::pElements,&MeshType::SetElements)
    .def("ElementsArray", &MeshType::ElementsArray, py::return_value_policy::reference_internal)
    .def("NumberOfElements", &MeshType::NumberOfElements)

    .def_property("Conditions", &MeshType::pConditions,&MeshType::SetConditions)
    .def("ConditionsArray", &MeshType::ConditionsArray, py::return_value_policy::reference_internal)
    .def("NumberOfConditions", &MeshType::NumberOfConditions)

    .def_property("Properties", &MeshType::pProperties,&MeshType::SetProperties)
    .def("PropertiesArray", &MeshType::PropertiesArray, py::return_value_policy::reference_internal)
    .def("NumberOfProperties", &MeshType::NumberOfProperties)

    .def("HasNode", &MeshType::HasNode)
    .def("HasProperties", &MeshType::HasProperties)
    .def("HasElement", &MeshType::HasElement)
    .def("HasCondition", &MeshType::HasCondition)
    .def("__str__", PrintObject<MeshType>)
    ;
}
}  // namespace Python.
} // Namespace Kratos
