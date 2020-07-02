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
#include "custom_constitutive/rht_concrete_plane_strain_2D_law.hpp"

namespace Kratos
{
    RHTConcretePlaneStrain2DLaw::RHTConcretePlaneStrain2DLaw()
  : RHTConcrete3DLaw()
  { }

    RHTConcretePlaneStrain2DLaw::RHTConcretePlaneStrain2DLaw(const RHTConcretePlaneStrain2DLaw& rOther)
  : RHTConcrete3DLaw(rOther)
  { }

  ConstitutiveLaw::Pointer RHTConcretePlaneStrain2DLaw::Clone() const
  {
    return Kratos::make_shared<RHTConcretePlaneStrain2DLaw>(*this);
  }

  RHTConcretePlaneStrain2DLaw::~RHTConcretePlaneStrain2DLaw()
  { }
} // Namespace Kratos
