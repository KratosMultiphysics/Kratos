//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:      JMCarbonell  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "umat_application.h"
#include "umat_application_variables.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

namespace Kratos
{

  namespace Python
  {

    using namespace boost::python;



    BOOST_PYTHON_MODULE(KratosUmatApplication)
    {

      class_<KratosUmatApplication,
	     KratosUmatApplication::Pointer,
	     bases<KratosApplication>, boost::noncopyable >("KratosUmatApplication")
	;

      AddCustomConstitutiveLawsToPython();

      //registering variables in python

    }


  }  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
