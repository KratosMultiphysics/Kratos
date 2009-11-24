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
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2008-10-23 11:34:04 $
//   Revision:            $Revision: 1.8 $
//
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
        
        boost::python::list GetValuesOnIntegrationPointsDouble( Element& dummy, 
                const Variable<double>& rVariable, const ProcessInfo& rCurrentProcessInfo )
        {
            boost::python::list values_list;
            IntegrationPointsArrayType integration_points = dummy.GetGeometry().IntegrationPoints( 
                    dummy.GetIntegrationMethod() );
            Vector values( integration_points.size() );
            dummy.CalculateOnIntegrationPoints( rVariable, values, rCurrentProcessInfo );
            for( unsigned int i=0; i<values.size(); i++ )
            {
                boost::python::list integration_point_value;
                integration_point_value.append( values[i] );
                values_list.append( integration_point_value );
            }
            return( values_list );
        }
        
        boost::python::list GetValuesOnIntegrationPointsArray1d( Element& dummy,  
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
        
        boost::python::list GetValuesOnIntegrationPointsVector( Element& dummy, 
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
        
        void SetValuesOnIntegrationPointsVector( Element& dummy,
                const Variable<Vector>& rVariable, boost::python::list values_list, int len_values_list_item, const ProcessInfo& rCurrentProcessInfo )
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
        
        boost::python::list GetValuesOnIntegrationPointsMatrix( Element& dummy, 
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
            
            class_<Element, Element::Pointer, bases<Element::BaseType > >("Element", init<int>())
                    //.def(init<int, const Point<3>& >())
                    //.add_property("Properties", &Element::pGetProperties(), &NodeType::SetId)
                    .def("__setitem__", SetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
                    .def("__getitem__", GetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
                    .def("SetValue", SetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
                    .def("GetValue", GetValueHelperFunction< Element, Variable< array_1d<double, 3>  > >)
                    
                    .def("__setitem__", SetValueHelperFunction< Element, Variable< Vector > >)
                    .def("__getitem__", GetValueHelperFunction< Element, Variable< Vector > >)
                    .def("SetValue", SetValueHelperFunction< Element, Variable< Vector > >)
                    .def("GetValue", GetValueHelperFunction< Element, Variable< Vector > >)
                    
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
                    
                    .def("GetNode", GetNodeFromElement )
                    .def("GetNodes", GetNodesFromElement )
                    .def("GetIntegrationPoints", GetIntegrationPointsFromElement )
                    .def("CalculateOnIntegrationPoints", CalculateOnIntegrationPointsVector)
                    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsDouble)
                    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsArray1d)
                    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsVector)
                    .def("GetValuesOnIntegrationPoints", GetValuesOnIntegrationPointsMatrix)
                    .def("SetValuesOnIntegrationPoints", SetValuesOnIntegrationPointsVector)
                    .def("ResetConstitutiveLaw", &Element::ResetConstitutiveLaw)
                    
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
            
            class_<Condition, Condition::Pointer, bases<Condition::BaseType > >("Condition", init<int>())
                    //.def(init<int, const Point<3>& >())
                    //.add_property("Properties", &Condition::pGetProperties(), &NodeType::SetId)
                    .def("__setitem__", SetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
                    .def("__getitem__", GetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
                    .def("SetValue", SetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)
                    .def("GetValue", GetValueHelperFunction< Condition, Variable< array_1d<double, 3>  > >)

                    .def("__setitem__", SetValueHelperFunction< Condition, Variable< Vector > >)
                    .def("__getitem__", GetValueHelperFunction< Condition, Variable< Vector > >)
                    .def("SetValue", SetValueHelperFunction< Condition, Variable< Vector > >)
                    .def("GetValue", GetValueHelperFunction< Condition, Variable< Vector > >)

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

                    .def("GetNode", GetNodeFromCondition )
                    .def("GetNodes", GetNodesFromCondition )

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
                    .def(self_ns::str(self))
                    ;
            
            PointerVectorSetPythonInterface<MeshType::ConditionsContainerType>::CreateInterface("ConditionsArray")
                    ;
            
            class_<MeshType, MeshType::Pointer, boost::noncopyable >("Mesh")
                    .add_property("Nodes", &MeshType::pNodes,&MeshType::SetNodes)
                    .def("NodesArray", &MeshType::NodesArray, return_internal_reference<>())
                    .add_property("Elements", &MeshType::pElements,&MeshType::SetElements)
                    .def("ElementsArray", &MeshType::ElementsArray, return_internal_reference<>())
                    .add_property("Conditions", &MeshType::pConditions,&MeshType::SetConditions)
                    .def("ConditionsArray", &MeshType::ConditionsArray, return_internal_reference<>())
                    .add_property("Properties", &MeshType::pProperties,&MeshType::SetProperties)
                    .def("PropertiesArray", &MeshType::PropertiesArray, return_internal_reference<>())
                    .def(self_ns::str(self))
                    ;
        }
    }  // namespace Python.
} // Namespace Kratos
