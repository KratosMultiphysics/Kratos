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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-04-24 10:30:22 $
//   Revision:            $Revision: 1.3 $
//
//


// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "python/add_constitutive_law_to_python.h"

// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "add_constitutive_law_to_python.h"
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
//#include "python/add_mesh_to_python.h"
//#include "python/pointer_vector_set_python_interface.h"
//#include "python/variable_indexing_python.h"

namespace Kratos
{
    namespace Python
    {
        using namespace boost::python;
        
        typedef ConstitutiveLaw ConstitutiveLawBaseType;
        typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
        
        void  AddConstitutiveLawToPython()
        {
            class_< ConstitutiveLaw, boost::noncopyable >
                    ("ConstitutiveLaw",
                     init<>() )
                    ;
            
            class_<Variable<ConstitutiveLawBaseType::Pointer> , bases<VariableData>, boost::noncopyable >("ConstitutiveLawVariable", no_init)
                    ;
        }
    }  // namespace Python.
}  // namespace Kratos.
