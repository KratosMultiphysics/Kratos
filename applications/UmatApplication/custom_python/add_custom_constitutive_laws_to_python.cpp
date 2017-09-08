//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:       LlMonforte  $
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

//constitutive laws
#include "custom_laws/small_strain_umat_3D_law.hpp"
#include "custom_laws/large_strain_umat_3D_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

typedef ConstitutiveLaw ConstitutiveLawBaseType;

void  AddCustomConstitutiveLawsToPython()
{
    class_< SmallStrainUmat3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "SmallStrainUmat3DLaw", init<>() )
    ;
    class_< LargeStrainUmat3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LargeStrainUmat3DLaw", init<>() )
    ;
}
}  // namespace Python.
} // Namespace Kratos
