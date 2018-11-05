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

typedef typename FrictionLaw::Pointer         FrictionLawPointerType;
//typedef std::vector<FrictionLawPointerType>     FrictionLawContainer;


void  AddCustomFrictionLawsToPython(pybind11::module& m)
{

  //Friction laws
  class_<FrictionLaw, FrictionLawPointerType>(m,"FrictionLaw")
      .def(init<>())
      .def("Clone", &FrictionLaw::Clone)
      .def("__repr__", &FrictionLaw::Info)
      DECLARE_HAS_THIS_TYPE_PROPERTIES_PYTHON_AS_POINTER(FrictionLaw)
      DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON_AS_POINTER(FrictionLaw)
      DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON_AS_POINTER(FrictionLaw)
      ;

  //to define it as a variable
  class_<Variable<FrictionLawPointerType>, VariableData>(m,"FrictionLawVariable")
      .def( "__repr__", &Variable<FrictionLawPointerType>::Info )
      ;

  class_<CoulombAdhesionFrictionLaw, typename CoulombAdhesionFrictionLaw::Pointer, FrictionLaw>(m,"CoulombAdhesionFrictionLaw")
      .def(init<>())
      ;

  class_<HardeningCoulombFrictionLaw, typename HardeningCoulombFrictionLaw::Pointer, CoulombAdhesionFrictionLaw>(m,"HardeningCoulombFrictionLaw")
      .def(init<>())
      ;

}

}  // namespace Python.

}  // namespace Kratos.
