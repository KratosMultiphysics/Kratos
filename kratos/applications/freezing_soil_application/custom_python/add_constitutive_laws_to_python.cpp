/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//   Last Modified by:    $Author: hurga $
//   Date:                $Date: 2009-02-02 14:03:23 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED


// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "custom_python/add_constitutive_laws_to_python.h"
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "constitutive_laws/viscoplastic_3d.h"
#include "constitutive_laws/clay_and_sand_model.h"
#include "constitutive_laws/casm.h"
#include "constitutive_laws/modified_barcelona_basic_model.h"
#include "constitutive_laws/freezing_soil_elastoplastic_model.h"
//#include "constitutive_laws/freezing_soil_elastoviscoplastic_model.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "python/add_mesh_to_python.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"


namespace Kratos
{

    namespace Python
    {

        using namespace boost::python;

        typedef ConstitutiveLaw ConstitutiveLawBaseType;
        typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
        typedef Properties::Pointer PropertiesPointer;

        typedef std::vector<ConstitutiveLaw::Pointer> MaterialsContainer;
        typedef ConstitutiveLaw::Pointer  ConstitutiveLawPointer;

        void  AddConstitutiveLawsToPython()
        {
            class_< Viscoplastic3D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
            ( "Viscoplastic3D",
              init<>() )
            ;
            class_< ClayAndSandModel, bases< ConstitutiveLawBaseType >, boost::noncopyable >
            ( "ClayAndSandModel",
              init<>() )
            ; 
            class_< CASM, bases< ConstitutiveLawBaseType >, boost::noncopyable >
            ( "CASM",
              init<>() )
            ; 
            class_< ModifiedBarcelonaBasicModel, bases< ConstitutiveLawBaseType >, boost::noncopyable >
            ( "ModifiedBarcelonaBasicModel",
              init<>() )
            ; 
            class_< FreezingSoilElastoplasticModel, bases< ConstitutiveLawBaseType >, boost::noncopyable >
            ( "FreezingSoilElastoplasticModel",
              init<>() )
            ; 
 //           class_< FreezingSoilElastoviscoplasticModel, bases< ConstitutiveLawBaseType >, boost::noncopyable >
   //         ( "FreezingSoilElastoviscoplasticModel",
     //         init<>() )
       //     ; 
        }
    }  // namespace Python.
}  // namespace Kratos.
#endif // KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED defined
