// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/elastic_isotropic_3d.h"

namespace Kratos
{
namespace Python
{

using namespace boost::python;


void  AddCustomConstitutiveLawsToPython()
{

    class_< LinearPlaneStress, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "LinearElasticPlaneStress2DLaw",
      init<>() )
    ;
    
    class_< LinearPlaneStrain, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "LinearElasticPlaneStrain2DLaw",
      init<>() )
    ;

    class_< ElasticIsotropic3D, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "LinearElasticIsotropic3DLaw",
      init<>() )
    ;
    
}

}  // namespace Python.
}  // namespace Kratos.
