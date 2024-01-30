//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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

namespace Kratos::Python
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

typedef Mesh<Node, Properties, Element, Condition> MeshType;
typedef MeshType::NodeType NodeType;
typedef MeshType::NodesContainerType NodesContainerType;
typedef Geometry<Node > GeometryType;
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
    for (std::size_t j = 0; j < Output.size(); j++) {
        result.append(Output[j]);
    }
    return result;
}

template< class TObject>
pybind11::list CalculateOnIntegrationPointsBool(
    TObject& dummy, const Variable<bool>& rVariable, const ProcessInfo& rProcessInfo)
{
    std::vector<bool> Output;
    dummy.CalculateOnIntegrationPoints(rVariable, Output, rProcessInfo);
    pybind11::list result;
    for (std::size_t j = 0; j < Output.size(); j++) {
        result.append(static_cast<int>(Output[j]));
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
    TObject& dummy, const Variable<TDataType>& rVariable, const std::vector<TDataType>& values, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(values.size() != dummy.GetGeometry().IntegrationPointsNumber())
        << "Sizes do not match. Size of values vector is: " << values.size() << ". The number of integration points is: "
        << dummy.GetGeometry().IntegrationPointsNumber() << std::endl;

    dummy.SetValuesOnIntegrationPoints(rVariable, values, rCurrentProcessInfo);
}

template< class TObject, std::size_t TSize>
void SetValuesOnIntegrationPointsArray1d(
    TObject& dummy,
    const Variable<array_1d<double, TSize>>& rVariable,
    pybind11::list values_list,
    const ProcessInfo& rCurrentProcessInfo)
{
    std::vector<array_1d<double, TSize>> values(values_list.size());
    for (std::size_t i = 0; i < values_list.size(); i++) {
        if (py::isinstance<array_1d<double, TSize>>(values_list[i])) {
            values[i] = (values_list[i]).cast<array_1d<double, TSize> >();
        } else if (py::isinstance<pybind11::list>(values_list[i]) ||
            py::isinstance<Vector>(values_list[i])) {
            Vector value = (values_list[i]).cast<Vector>();
            KRATOS_ERROR_IF(value.size() != TSize)
                << " parsed vector is not of size " << TSize << ". Size of vector: " << value.size() << std::endl;
            values[i] = value;
        } else {
            KRATOS_ERROR << "expecting a list of array_1d<double, " << TSize << ">" << std::endl;
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
    for( std::size_t i=0; i<integration_points.size(); i++ ) {
        if(py::isinstance<Vector>(values_list[i]))
            values[i] = (values_list[i]).cast<Vector>();
        else
            KRATOS_ERROR << "expecting a list of vectors";
    }
    dummy.SetValuesOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
}

template<class TEntityType, class TDataType >
TDataType EntityCalculateInterface(TEntityType& dummy, Variable<TDataType>& rVariable, const ProcessInfo& rCurrentProcessInfo)
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
    for( unsigned int i=0; i<integration_points.size(); i++ ) {
        if(py::isinstance<ConstitutiveLaw::Pointer>(values_list[i]))
            values[i] = (values_list[i]).cast<ConstitutiveLaw::Pointer>();
        else
            KRATOS_ERROR << "expecting a list of ConstitutiveLaw::Pointer";
     }
    dummy.SetValuesOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
}

template<class TEntityType>
void EntityCalculateRightHandSide(
    TEntityType& dummy,
    Vector& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

template<class TEntityType>
void EntityCalculateLocalSystem(
    TEntityType& dummy,
    Matrix& rLeftHandSideMatrix,
    Vector& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
}

template<class TEntityType>
void EntityCalculateMassMatrix(
    TEntityType& dummy,
    Matrix& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);
}

template<class TEntityType>
void EntityCalculateDampingMatrix(
    TEntityType& dummy,
    Matrix& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
}

void ElementCalculateLumpedMassVector(Element& dummy,
                                      Vector& rMassVector,
                                      const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateLumpedMassVector(rMassVector, rCurrentProcessInfo);
}

template<class TEntityType>
void EntityCalculateFirstDerivativesLHS(
    TEntityType& dummy,
    Matrix& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateFirstDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template<class TEntityType>
void EntityCalculateSecondDerivativesLHS(
    TEntityType& dummy,
    Matrix& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateSecondDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template<class TEntityType>
void EntityCalculateLocalVelocityContribution(
    TEntityType& dummy,
    Matrix& rDampingMatrix,
    Vector& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateLocalVelocityContribution(rDampingMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template<class TEntityType>
void EntityInitialize(
    TEntityType& dummy,
    const ProcessInfo& rCurrentProcessInfo)
{
    dummy.Initialize(rCurrentProcessInfo);
}

template<class TEntityType, class TDataType>
void EntityCalculateSensitivityMatrix(
    TEntityType& dummy,
    const Variable<TDataType>& rDesignVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    dummy.CalculateSensitivityMatrix(rDesignVariable,rOutput,rCurrentProcessInfo);
}

template<class TEntityType>
void EntityGetFirstDerivativesVector1(
    const TEntityType& dummy,
    Vector& rOutput)
{
    dummy.GetFirstDerivativesVector(rOutput,0);
}

template<class TEntityType>
void EntityGetFirstDerivativesVector2(
    const TEntityType& dummy,
    Vector& rOutput,
    int step)
{
    dummy.GetFirstDerivativesVector(rOutput,step);
}

template<class TEntityType>
void EntityGetSecondDerivativesVector1(
    const TEntityType& dummy,
    Vector& rOutput)
{
    dummy.GetSecondDerivativesVector(rOutput,0);
}

template<class TEntityType>
void EntityGetSecondDerivativesVector2(
    const TEntityType& dummy,
    Vector& rOutput,
    int step)
{
    dummy.GetSecondDerivativesVector(rOutput,step);
}

void  AddMeshToPython(pybind11::module& m)
{
    py::class_<GeometricalObject, GeometricalObject::Pointer, IndexedObject, Flags>(m,"GeometricalObject")
    .def(py::init<Kratos::GeometricalObject::IndexType>())
    .def("IsActive", &GeometricalObject::IsActive)
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
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPointsBool<Element>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, int>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, double>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, array_1d<double, 3>>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, array_1d<double, 4>>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, array_1d<double, 6>>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, array_1d<double, 9>>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, Vector>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Element, Matrix>)
    // GetValuesOnIntegrationPoints
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPoints<Element>)
    // SetValuesOnIntegrationPoints
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPoints<Element, bool>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPoints<Element, int>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsVector<Element>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsConstitutiveLaw)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPoints<Element, double>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Element, 3>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Element, 4>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Element, 6>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Element, 9>)
    .def("ResetConstitutiveLaw", &Element::ResetConstitutiveLaw)
    .def("Calculate", &EntityCalculateInterface<Element, double>)
    .def("Calculate", &EntityCalculateInterface<Element, array_1d<double,3> >)
    .def("Calculate", &EntityCalculateInterface<Element, Vector >)
    .def("Calculate", &EntityCalculateInterface<Element, Matrix >)
    .def("CalculateLumpedMassVector", &ElementCalculateLumpedMassVector)
    .def("CalculateMassMatrix", &EntityCalculateMassMatrix<Element>)
    .def("CalculateDampingMatrix", &EntityCalculateDampingMatrix<Element>)
    .def("CalculateLocalSystem", &EntityCalculateLocalSystem<Element>)
    .def("CalculateRightHandSide", &EntityCalculateRightHandSide<Element>)
    .def("CalculateFirstDerivativesLHS", &EntityCalculateFirstDerivativesLHS<Element>)
    .def("CalculateSecondDerivativesLHS", &EntityCalculateSecondDerivativesLHS<Element>)
    .def("CalculateLocalVelocityContribution", &EntityCalculateLocalVelocityContribution<Element>)
    .def("GetFirstDerivativesVector", &EntityGetFirstDerivativesVector1<Element>)
    .def("GetFirstDerivativesVector", &EntityGetFirstDerivativesVector2<Element>)
    .def("GetSecondDerivativesVector", &EntityGetSecondDerivativesVector1<Element>)
    .def("GetSecondDerivativesVector", &EntityGetSecondDerivativesVector2<Element>)
    .def("CalculateSensitivityMatrix", &EntityCalculateSensitivityMatrix<Element, double>)
    .def("CalculateSensitivityMatrix", &EntityCalculateSensitivityMatrix<Element, array_1d<double,3>>)

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
    .def("Initialize", &EntityInitialize<Element>)
    .def("EquationIdVector", [](const Element& self, const ProcessInfo& rProcessInfo){
        Element::EquationIdVectorType ids;
        self.EquationIdVector(ids,rProcessInfo);
        return ids;
    })
    .def("GetDofList", [](const Element& self, const ProcessInfo& rProcessInfo){
        std::vector<Dof<double>*> dofs_list;
        self.GetDofList(dofs_list,rProcessInfo);
        return dofs_list;
    })
    .def("CalculateLocalSystem", &Element::CalculateLocalSystem)
    .def("GetSpecifications", &Element::GetSpecifications)
    .def("Info", &Element::Info)
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
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPointsBool<Condition>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, int>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, double>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, array_1d<double, 3>>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, array_1d<double, 4>>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, array_1d<double, 6>>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, array_1d<double, 9>>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, Vector>)
    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPoints<Condition, Matrix>)
    // GetValuesOnIntegrationPoints
    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPoints<Condition>)
    // SetValuesOnIntegrationPoints
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPoints<Condition, bool>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPoints<Condition, int>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPoints<Condition, double>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsVector<Condition>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Condition, 3>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Condition, 4>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Condition, 6>)
    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsArray1d<Condition, 9>)
    //.def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsConstitutiveLaw)

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
    .def("Calculate", &EntityCalculateInterface<Condition, double>)
    .def("Calculate", &EntityCalculateInterface<Condition, array_1d<double,3> >)
    .def("Calculate", &EntityCalculateInterface<Condition, Vector >)
    .def("Calculate", &EntityCalculateInterface<Condition, Matrix >)

    .def("Initialize", &EntityInitialize<Condition>)
    .def("EquationIdVector", [](const Condition& self, const ProcessInfo& rProcessInfo){
        Condition::EquationIdVectorType ids;
        self.EquationIdVector(ids,rProcessInfo);
        return ids;
    })
    .def("GetDofList", [](const Condition& self, const ProcessInfo& rProcessInfo){
        std::vector<Dof<double>*> dofs_list;
        self.GetDofList(dofs_list,rProcessInfo);
        return dofs_list;
    })
    .def("CalculateMassMatrix", &EntityCalculateMassMatrix<Condition>)
    .def("CalculateDampingMatrix", &EntityCalculateDampingMatrix<Condition>)
    .def("CalculateLocalSystem", &EntityCalculateLocalSystem<Condition>)
    .def("CalculateRightHandSide", &EntityCalculateRightHandSide<Condition>)
    .def("CalculateFirstDerivativesLHS", &EntityCalculateFirstDerivativesLHS<Condition>)
    .def("CalculateSecondDerivativesLHS", &EntityCalculateSecondDerivativesLHS<Condition>)
    .def("CalculateLocalVelocityContribution", &EntityCalculateLocalVelocityContribution<Condition>)
    .def("GetFirstDerivativesVector", &EntityGetFirstDerivativesVector1<Condition>)
    .def("GetFirstDerivativesVector", &EntityGetFirstDerivativesVector2<Condition>)
    .def("GetSecondDerivativesVector", &EntityGetSecondDerivativesVector1<Condition>)
    .def("GetSecondDerivativesVector", &EntityGetSecondDerivativesVector2<Condition>)
    .def("CalculateSensitivityMatrix", &EntityCalculateSensitivityMatrix<Condition, double>)
    .def("CalculateSensitivityMatrix", &EntityCalculateSensitivityMatrix<Condition, array_1d<double,3>>)
    .def("GetSpecifications", &Condition::GetSpecifications)
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
    .def("HasMasterSlaveConstraint ", &MeshType::HasMasterSlaveConstraint )
    .def("__str__", PrintObject<MeshType>)
    ;
}
}  // namespace Kratos::Python.
