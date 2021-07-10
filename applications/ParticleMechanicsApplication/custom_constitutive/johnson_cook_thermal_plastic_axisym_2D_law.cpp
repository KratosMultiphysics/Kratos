//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/johnson_cook_thermal_plastic_axisym_2D_law.hpp"

namespace Kratos
{
    JohnsonCookThermalPlastic2DAxisymLaw::JohnsonCookThermalPlastic2DAxisymLaw()
  : JohnsonCookThermalPlastic2DPlaneStrainLaw()
  { }

    JohnsonCookThermalPlastic2DAxisymLaw::JohnsonCookThermalPlastic2DAxisymLaw(const JohnsonCookThermalPlastic2DAxisymLaw& rOther)
  : JohnsonCookThermalPlastic2DPlaneStrainLaw(rOther)
  { }

  ConstitutiveLaw::Pointer JohnsonCookThermalPlastic2DAxisymLaw::Clone() const
  {
    return Kratos::make_shared<JohnsonCookThermalPlastic2DAxisymLaw>(*this);
  }

  JohnsonCookThermalPlastic2DAxisymLaw::~JohnsonCookThermalPlastic2DAxisymLaw()
  { }
} // Namespace Kratos
