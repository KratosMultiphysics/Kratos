//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/johnson_cook_thermal_plastic_plane_strain_2D_law.hpp"

namespace Kratos
{
    JohnsonCookThermalPlastic2DPlaneStrainLaw::JohnsonCookThermalPlastic2DPlaneStrainLaw()
  : JohnsonCookThermalPlastic3DLaw()
  { }

    JohnsonCookThermalPlastic2DPlaneStrainLaw::JohnsonCookThermalPlastic2DPlaneStrainLaw(const JohnsonCookThermalPlastic2DPlaneStrainLaw& rOther)
  : JohnsonCookThermalPlastic3DLaw(rOther)
  { }

  ConstitutiveLaw::Pointer JohnsonCookThermalPlastic2DPlaneStrainLaw::Clone() const
  {
    return Kratos::make_shared<JohnsonCookThermalPlastic2DPlaneStrainLaw>(*this);
  }

  JohnsonCookThermalPlastic2DPlaneStrainLaw::~JohnsonCookThermalPlastic2DPlaneStrainLaw()
  { }
} // Namespace Kratos
