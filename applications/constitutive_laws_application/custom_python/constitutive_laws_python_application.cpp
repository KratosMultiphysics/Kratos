//
//   Project Name:        Kratos
//   Last modified by:    $Author: stasch $
//   Date:                $Date: 2007-09-26 13:57:49 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "constitutive_laws_application.h"
#include "custom_python/add_constitutive_laws_to_python.h"

namespace Kratos
{

    namespace Python
    {

        using namespace boost::python;

        BOOST_PYTHON_MODULE( KratosConstitutiveLawsApplication )
        {

            class_ < KratosConstitutiveLawsApplication,
            KratosConstitutiveLawsApplication::Pointer,
            bases<KratosApplication>, boost::noncopyable > ( "KratosConstitutiveLawsApplication" )
            ;

            AddConstitutiveLawsToPython();
        }


    }  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
