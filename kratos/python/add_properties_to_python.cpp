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
//   Date:                $Date: 2007-04-24 10:30:22 $
//   Revision:            $Revision: 1.4 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// Project includes
#include "includes/define.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "includes/constitutive_law.h"
//#include "python/add_mesh_to_python.h"
#include "python/pointer_vector_set_python_interface.h"
//#include "python/variable_indexing_python.h"

namespace Kratos
{
    namespace Python
    {
        using namespace boost::python;
        
        typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
        typedef ConstitutiveLaw ConstitutiveLawBaseType;
        
        template< class TContainerType, class TVariableType > void SetValueHelperFunction1(
                TContainerType& el, 
                const TVariableType& rVar,
                const typename TVariableType::Type& Data)
        {
            el.SetValue(rVar,Data);
        }
        
        template< class TContainerType, class TVariableType > 
                typename TVariableType::Type GetValueHelperFunction1( TContainerType& el, 
                        const TVariableType& rVar )
        {
            return el.GetValue(rVar);
        }
        
        void  AddPropertiesToPython()
        {
            class_<Properties, Properties::Pointer, bases<Properties::BaseType > >("Properties", init<int>())
                    .def("__setitem__", SetValueHelperFunction1< Element, Variable< array_1d<double, 6> > >)
                    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< array_1d<double, 6> > >)
                    .def("SetValue", SetValueHelperFunction1< Properties, Variable< array_1d<double, 6> > >)
                    .def("GetValue", GetValueHelperFunction1< Properties, Variable< array_1d<double, 6> > >)

                    
                    .def("__setitem__", SetValueHelperFunction1< Element, Variable< array_1d<double, 3> > >)
                    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< array_1d<double, 3> > >)
                    .def("SetValue", SetValueHelperFunction1< Properties, Variable< array_1d<double, 3> > >)
                    .def("GetValue", GetValueHelperFunction1< Properties, Variable< array_1d<double, 3> > >)
                    
                    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< Vector > >)
                    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< Vector > >)
                    .def("SetValue", SetValueHelperFunction1< Properties, Variable< Vector > >)
                    .def("GetValue", GetValueHelperFunction1< Properties, Variable< Vector > >)
                    
                    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< Matrix > >)
                    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< Matrix > >)
                    .def("SetValue", SetValueHelperFunction1< Properties, Variable< Matrix > >)
                    .def("GetValue", GetValueHelperFunction1< Properties, Variable< Matrix > >)
                    
                    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< int > >)
                    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< int > >)
                    .def("SetValue", SetValueHelperFunction1< Properties, Variable< int > >)
                    .def("GetValue", GetValueHelperFunction1< Properties, Variable< int > >)
                    
                    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< double > >)
                    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< double > >)
                    .def("SetValue", SetValueHelperFunction1< Properties, Variable< double > >)
                    .def("GetValue", GetValueHelperFunction1< Properties, Variable< double > >)
                    
                    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)
                    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)
                    .def("SetValue", SetValueHelperFunction1< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)
                    .def("GetValue", GetValueHelperFunction1< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)
                    
                    .def(self_ns::str(self))
                    ;
            
            PointerVectorSetPythonInterface<MeshType::PropertiesContainerType>::CreateInterface("PropertiesArray")
                    ;
        }
    }  // namespace Python.
} // Namespace Kratos
