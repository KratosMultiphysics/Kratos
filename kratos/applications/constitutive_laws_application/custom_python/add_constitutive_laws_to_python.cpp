//
//   Project Name:        Kratos
//   Last modified by:    $Author: nagel $
//   Date:                $Date: 2009-01-12 08:17:36 $
//   Revision:            $Revision: 1.10 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "custom_python/add_constitutive_laws_to_python.h"

#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/condition.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"

//constitutive laws
#include "includes/constitutive_law.h"
#include "constitutive_laws/umat.h"

namespace Kratos
{

    namespace Python
    {

        using namespace boost::python;

        typedef ConstitutiveLaw ConstitutiveLawBaseType;
        typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;

        void  AddConstitutiveLawsToPython()
        {
            class_< ConstitutiveLaw, boost::noncopyable >
            ( "ConstitutiveLaw", init<>() )
            ;

            class_< Umat, bases< ConstitutiveLawBaseType >, boost::noncopyable >
            ( "Umat", init<>() )
            ;
        }
    }  // namespace Python.
} // Namespace Kratos
