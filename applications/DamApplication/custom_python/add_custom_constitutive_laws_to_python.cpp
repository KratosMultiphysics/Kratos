//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:                 $
//   Revision:            $Revision:          $
//

// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"

#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"
#include "python/add_mesh_to_python.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//constitutive laws
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

typedef ConstitutiveLaw                  ConstitutiveLawBaseType;

void  AddCustomConstitutiveLawsToPython()
{
    class_< ThermalLinearElastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >( "ThermalLinearElastic3DLaw",init<>() );

    class_< ThermalLinearElastic2DPlaneStrain, bases< ConstitutiveLawBaseType >, boost::noncopyable >( "ThermalLinearElastic2DPlaneStrain",init<>() );

    class_< ThermalLinearElastic2DPlaneStress, bases< ConstitutiveLawBaseType >, boost::noncopyable >( "ThermalLinearElastic2DPlaneStress",init<>() );
}

}  // namespace Python.
}  // namespace Kratos.
