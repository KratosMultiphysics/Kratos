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
#include "custom_constitutive/truss_constitutive_law.h"
#include "custom_constitutive/beam_constitutive_law.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/axisym_elastic_isotropic.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_stress_2d.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_strain_2d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"
#include "custom_constitutive/linear_elastic_orthotropic_2D_law.h"

#include "custom_constitutive/linear_j2_plasticity_3d.h"
#include "custom_constitutive/linear_j2_plasticity_plane_strain_2d.h"

namespace Kratos
{
namespace Python
{

using namespace boost::python;


void  AddCustomConstitutiveLawsToPython()
{

    class_< TrussConstitutiveLaw, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "TrussConstitutiveLaw", init<>() )
    ;

    class_< BeamConstitutiveLaw, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "BeamConstitutiveLaw", init<>() )
    ;

    class_< LinearPlaneStress, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "LinearElasticPlaneStress2DLaw", init<>() )
    ;
    
    class_< LinearPlaneStrain, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "LinearElasticPlaneStrain2DLaw", init<>() )
    ;

    class_< ElasticIsotropic3D, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "LinearElastic3DLaw", init<>() )
    ;

    class_< AxisymElasticIsotropic, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "LinearElasticAxisym2DLaw", init<>() )
    ;

    class_< HyperElasticIsotropicKirchhoff3D, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "KirchhoffSaintVenant3DLaw", init<>() )
    ;

    class_< HyperElasticIsotropicKirchhoffPlaneStress2D, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "KirchhoffSaintVenantPlaneStress2DLaw", init<>() )
    ;

    class_< HyperElasticIsotropicKirchhoffPlaneStrain2D, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "KirchhoffSaintVenantPlaneStrain2DLaw", init<>() )
    ;

    class_< HyperElasticIsotropicNeoHookean3D, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "HyperElastic3DLaw", init<>() )
    ;
    
    class_< HyperElasticIsotropicNeoHookeanPlaneStrain2D, bases< ConstitutiveLaw >, boost::noncopyable >
    ( "HyperElasticPlaneStrain2DLaw", init<>() )
    ;

    class_< LinearElasticOrthotropic2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >
    ("LinearElasticOrthotropic2DLaw", init<>())
    ;

    class_< LinearJ2PlasticityPlaneStrain2D, bases< ConstitutiveLaw >, boost::noncopyable >
    ("LinearJ2PlasticityPlaneStrain2DLaw", init<>())
    ;

    class_< LinearJ2Plasticity3D, bases< ConstitutiveLaw >, boost::noncopyable >
    ("LinearJ2Plasticity3DLaw", init<>())
    ;
}

}  // namespace Python.
}  // namespace Kratos.
