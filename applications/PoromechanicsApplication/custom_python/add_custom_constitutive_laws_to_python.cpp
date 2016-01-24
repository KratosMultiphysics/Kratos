//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
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
#include "custom_constitutive/linear_elastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_2D_plane_strain_law.hpp"
#include "custom_constitutive/linear_elastic_2D_plane_stress_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

typedef ConstitutiveLaw                  ConstitutiveLawBaseType;


void  AddCustomConstitutiveLawsToPython()
{    
    class_< LinearElastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >( "LinearElastic3DLaw",init<>() );

    class_< LinearElastic2DPlaneStrainLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >( "LinearElastic2DPlaneStrainLaw",init<>() );

    class_< LinearElastic2DPlaneStressLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >( "LinearElastic2DPlaneStressLaw",init<>() );
}

}  // namespace Python.
}  // namespace Kratos.
