//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes
#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
// External includes


// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"


//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//hardening laws

//yield criteria

//flow rules

//constitutive laws


namespace Kratos
{

  namespace Python
  {

    using namespace pybind11;

    typedef ConstitutiveLaw                  ConstitutiveLawBaseType;

    void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
    {

    }

  }  // namespace Python.
}  // namespace Kratos.
