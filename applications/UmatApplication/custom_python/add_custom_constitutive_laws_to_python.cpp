//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:      JMCarbonell  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//


// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "includes/properties.h"

#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"


//models
#include "custom_models/hypoplastic_umat_small_strain_model.hpp"
#include "custom_models/von_mises_umat_small_strain_model.hpp"
#include "custom_models/von_mises_umat_large_strain_model.hpp"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

typedef ConstitutiveLaw      ConstitutiveLawBaseType;
typedef ConstitutiveModel  ConstitutiveModelBaseType;

void  AddCustomConstitutiveLawsToPython()
{

    // models
    class_< VonMisesSmallStrainUmatModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
    	( "VonMisesSmallStrainUmatModel",
    	  init<>() )
     	;
    class_< VonMisesLargeStrainUmatModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
    	( "VonMisesLargeStrainUmatModel",
    	  init<>() )
     	;
    class_< HypoplasticSmallStrainUmatModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
    	( "HypoplasticSmallStrainUmatModel",
    	  init<>() )
     	;
}
}  // namespace Python.
} // Namespace Kratos
