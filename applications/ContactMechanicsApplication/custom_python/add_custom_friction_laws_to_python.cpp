//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                August 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_python/add_custom_friction_laws_to_python.h"

//friction laws
#include "custom_friction/friction_law.hpp"
#include "custom_friction/coulomb_adhesion_friction_law.hpp"
#include "custom_friction/hardening_coulomb_friction_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

typedef typename FrictionLaw::Pointer         FrictionLawPointer;
typedef std::vector<FrictionLaw::Pointer>   FrictionLawContainer;

void Push_Back_Friction_Laws( FrictionLawContainer& ThisFrictionLawContainer,
                              FrictionLawPointer ThisFrictionLaw )
{
  ThisFrictionLawContainer.push_back( ThisFrictionLaw );
}

void  AddCustomFrictionLawsToPython(pybind11::module& m)
{

  class_<FrictionLawContainer>(m,"FrictionLawContainer")
      .def( init<>() )
      .def( "PushBack", Push_Back_Friction_Laws )
      ;

  class_<Variable<FrictionLawPointer>, VariableData>(m,"FrictionLawVariable")
      ;

  //Friction laws
  class_< FrictionLaw, typename FrictionLaw::Pointer>(m,"FrictionLaw")
      .def( init<>() )
      .def("Clone",&FrictionLaw::Clone)
      ;

  class_< CoulombAdhesionFrictionLaw, typename CoulombAdhesionFrictionLaw::Pointer, FrictionLaw>(m,"CoulombAdhesionFrictionLaw")
      .def( init<>() )
      ;

  class_< HardeningCoulombFrictionLaw, typename HardeningCoulombFrictionLaw::Pointer, CoulombAdhesionFrictionLaw>(m,"HardeningCoulombFrictionLaw")
      .def( init<>() )
      ;

}

}  // namespace Python.

}  // namespace Kratos.
