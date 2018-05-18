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

      KRATOS_REGISTER_IN_PYTHON_VARIABLE( ALPHA )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( BETA )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( MF )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( CC )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( MM )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( RHOS )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( KSIS )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( RHOM )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( KSIM )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( PC0 )

      KRATOS_REGISTER_IN_PYTHON_VARIABLE( VOID_RATIO )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( PS )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( PM )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( PLASTIC_MULTIPLIER )
    }


  }  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
