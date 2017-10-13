/*
==============================================================================
KratosFreezingSoilApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author:  $
//   Date:                $Date: $
//   Revision:            $Revision: 1.3 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "freezing_soil_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_constitutive_laws_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosFreezingSoilApplication)
{

    class_<KratosFreezingSoilApplication,
    KratosFreezingSoilApplication::Pointer,
    bases<KratosApplication>, boost::noncopyable >("KratosFreezingSoilApplication")
    ;

    AddCustomStrategiesToPython();
    AddConstitutiveLawsToPython();
//     AddCustomUtilitiesToPython();
//     AddCustomProcessesToPython();

    //registering variables in python

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE_NULL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE_EINS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE_NULL_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE_EINS_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE_NULL_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE_EINS_ACCELERATION )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FACE_WATER_FLUX )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(FACE_LOAD_PRESSURE )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(KRATOS_WATCH_FLAG )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ASSIGN_PRESTRESS_FLAG )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PLASTIC_FLAG )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PRESTRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EPSILON )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(LINEAR_STRAIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EFFECTIVE_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TOTAL_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SUCTION )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ICE_MASS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(WATER_MASS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ICE_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ICE_SATURATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ICE_VOLUME_FRACTION )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(WATER_FLOW)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ICE_FLOW)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(HEAT_FLOW)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRAIN_TENSOR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SCALE_U )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SCALE_O ) 

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MECH_DISSIPATION )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ICE_DENSITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(WATER_DENSITY )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PLASTICITY_INDICATOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INSITU_STRESS_SCALE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PRECONSOLIDATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EQUIVALENT_VOLUMETRIC_STRAIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EQUIVALENT_DEVIATORIC_STRAIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EQUIVALENT_VOLUMETRIC_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EQUIVALENT_DEVIATORIC_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(LOG_EQUIVALENT_VOLUMETRIC_STRESS )
     
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELEMENT_PARAMETERS )

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
