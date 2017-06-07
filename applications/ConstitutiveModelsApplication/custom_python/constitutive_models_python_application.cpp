//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "constitutive_models_application.h"
#include "constitutive_models_application_variables.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosConstitutiveModelsApplication)
  {

    class_<KratosConstitutiveModelsApplication,
	   KratosConstitutiveModelsApplication::Pointer,
	   bases<KratosApplication>, boost::noncopyable >("KratosConstitutiveModelsApplication")
      ;

    AddCustomConstitutiveLawsToPython(); 
    AddCustomUtilitiesToPython();
    AddCustomProcessesToPython();

  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
