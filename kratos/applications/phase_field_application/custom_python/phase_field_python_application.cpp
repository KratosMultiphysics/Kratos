/*
==============================================================================
KratosMultiphaseApplication 
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
hgbk2008@gmail.com
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
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Apr 30, 2013 $
//   Revision:            $Revision: 1.1 $
//
//
//Change log:

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "custom_utilities/isotropic_tensor_utility_tester.h"
#include "phase_field_application.h"

namespace Kratos
{


namespace Python
{

    using namespace boost::python;

    BOOST_PYTHON_MODULE(KratosPhaseFieldApplication)
    {
        class_<KratosPhaseFieldApplication,  KratosPhaseFieldApplication::Pointer, bases<KratosApplication>, boost::noncopyable>
        ("KratosPhaseFieldApplication");

        KRATOS_REGISTER_IN_PYTHON_VARIABLE(PHASE_FIELD)
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(PHASE_FIELD_DUAL_VARIABLE)
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(PHASE_FIELD_GRADIENT)
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(LENGTH_SCALE)
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(REFERENCE_ENERGY_DENSITY)
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(INTEGRATION_POINT_GLOBAL_COORDINATES)
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(PHASE_FIELD_ORDER)

        class_<IsotropicTensorUtilityTester, boost::noncopyable>
        ("IsotropicTensorUtilityTester", init<>())
        .def("Test1", &IsotropicTensorUtilityTester::Test1)
        .def("Test1_1", &IsotropicTensorUtilityTester::Test1_1)
        .def("Test1_2", &IsotropicTensorUtilityTester::Test1_2)
        .def("Test2", &IsotropicTensorUtilityTester::Test2)
        .def("Test2_1", &IsotropicTensorUtilityTester::Test2_1)
        .def("Test2_2", &IsotropicTensorUtilityTester::Test2_2)
        ;

    }
  
}  // namespace Python.

  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
