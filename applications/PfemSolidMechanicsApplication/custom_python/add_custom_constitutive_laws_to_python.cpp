//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <pybind11/stl.h>

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

// Constitutive laws
#include "custom_constitutive/borja_hencky_cam_clay_3D_law.hpp"
#include "custom_constitutive/borja_hencky_cam_clay_axisym_2D_law.hpp"
#include "custom_constitutive/borja_hencky_cam_clay_plane_strain_2D_law.hpp"

#include "custom_constitutive/hencky_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_J2_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_tresca_axisym_2D_law.hpp"
#include "custom_constitutive/new_hencky_tresca_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_tresca_plane_strain_2D_law.hpp"
#include "custom_constitutive/new_hencky_tresca_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_tresca_3D_law.hpp"

#include "custom_constitutive/hencky_U_P_J2_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_U_P_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_U_P_Tresca_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_U_P_Tresca_plane_strain_2D_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{

  typedef typename FlowRule::Pointer                        FlowRulePointer; 
  typedef typename YieldCriterion::Pointer            YieldCriterionPointer; 
  typedef typename HardeningLaw::Pointer                HardeningLawPointer;
  
  // Constitutive Laws for soil plasticity
  class_<BorjaHenckyCamClayPlastic3DLaw, typename BorjaHenckyCamClayPlastic3DLaw::Pointer, ConstitutiveLaw >
      (m, "BorjaHenckyCamClayPlastic3DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;
  
  class_<BorjaHenckyCamClayPlasticAxisym2DLaw, typename BorjaHenckyCamClayPlasticAxisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "BorjaHenckyCamClayPlasticAxisym2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;
  class_<BorjaHenckyCamClayPlasticPlaneStrain2DLaw, typename BorjaHenckyCamClayPlasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "BorjaHenckyCamClayPlasticPlaneStrain2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;
  class_<HenckyJ2PlasticPlaneStrain2DLaw, typename HenckyJ2PlasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "HenckyJ2PlasticPlaneStrain2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;
       
  class_<HenckyJ2PlasticAxisym2DLaw, typename HenckyJ2PlasticAxisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "HenckyJ2PlasticAxisym2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;

  class_<HenckyPlasticUPJ2Axisym2DLaw, typename HenckyPlasticUPJ2Axisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "HenckyPlasticUPJ2Axisym2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;

  class_<HenckyPlasticUPJ2PlaneStrain2DLaw, typename HenckyPlasticUPJ2PlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "HenckyPlasticUPJ2PlaneStrain2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;

  class_<HenckyPlasticUPTrescaAxisym2DLaw, typename HenckyPlasticUPTrescaAxisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "HenckyPlasticUPTrescaAxisym2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;

  class_<HenckyPlasticUPTrescaPlaneStrain2DLaw, typename HenckyPlasticUPTrescaPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "HenckyPlasticUPTrescaPlaneStrain2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;

  class_<HenckyTrescaPlasticAxisym2DLaw, typename HenckyTrescaPlasticAxisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "HenckyTrescaPlasticAxisym2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;

 class_<NewHenckyTrescaPlasticAxisym2DLaw, typename NewHenckyTrescaPlasticAxisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "NewHenckyTrescaPlasticAxisym2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;
  class_<HenckyTresca3DLaw, typename HenckyTresca3DLaw::Pointer, ConstitutiveLaw >
      (m, "HenckyTresca3DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;
  class_<HenckyTrescaPlasticPlaneStrain2DLaw, typename HenckyTrescaPlasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "HenckyTrescaPlasticPlaneStrain2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;

 class_<NewHenckyTrescaPlasticPlaneStrain2DLaw, typename NewHenckyTrescaPlasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "NewHenckyTrescaPlasticPlaneStrain2DLaw")
      .def( init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      ;
}
    
}  // namespace Python.
}  // namespace Kratos. 
