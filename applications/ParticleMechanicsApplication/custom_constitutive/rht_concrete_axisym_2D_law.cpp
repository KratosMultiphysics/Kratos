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
#include "custom_constitutive/rht_concrete_axisym_2D_law.hpp"

namespace Kratos
{
    RHTConcreteAxisym2DLaw::RHTConcreteAxisym2DLaw()
  : RHTConcretePlaneStrain2DLaw()
  { }

    RHTConcreteAxisym2DLaw::RHTConcreteAxisym2DLaw(const RHTConcreteAxisym2DLaw& rOther)
  : RHTConcretePlaneStrain2DLaw(rOther)
  { }

  ConstitutiveLaw::Pointer RHTConcreteAxisym2DLaw::Clone() const
  {
    return Kratos::make_shared<RHTConcreteAxisym2DLaw>(*this);
  }

  RHTConcreteAxisym2DLaw::~RHTConcreteAxisym2DLaw()
  { }
} // Namespace Kratos
