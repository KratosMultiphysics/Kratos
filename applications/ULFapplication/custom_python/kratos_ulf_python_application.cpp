/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).
Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-23 12:50:01 $
//   Revision:            $Revision: 1.6 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "ULF_application.h"
#include "custom_python/add_custom_io_to_python.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_processes_to_python.h"


namespace Kratos
{

namespace Python
{

namespace py = pybind11;


PYBIND11_MODULE(KratosULFApplication,m)
{

    py::class_<KratosULFApplication,
           KratosULFApplication::Pointer,
           KratosApplication>(m,"KratosULFApplication")
           .def(py::init<>())
           ;
    AddCustomUtilitiesToPython(m);
    //AddCustomIOToPython();
    AddCustomStrategiesToPython(m);
    AddProcessesToPython(m);

    //registering variables in python

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,VAUX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PRESSURE_FORCE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DISP_FRAC);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TAUONE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TAUTWO);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_LENGTH);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MEAN_CURVATURE_2D);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TRIPLE_POINT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,CONTACT_ANGLE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,CONTACT_ANGLE_STATIC );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,SURFACE_TENSION_COEF );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,MEAN_CURVATURE_3D );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,GAUSSIAN_CURVATURE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,PRINCIPAL_CURVATURE_1 );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,PRINCIPAL_CURVATURE_2 );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SUBSCALE_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,VISCOUS_STRESSX );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,VISCOUS_STRESSY );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,VISCOUS_STRESSZ );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,PRINCIPAL_DIRECTION_1 );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,PRINCIPAL_DIRECTION_2 );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,NORMAL_GEOMETRIC );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,ADHESION_FORCE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,NORMAL_EQUILIBRIUM );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,NORMAL_CONTACT_LINE_EQUILIBRIUM );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,NORMAL_TRIPLE_POINT );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,NORMAL_CONTACT_LINE );
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,SOLID_FRACTION_GRADIENT );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m,SOLID_FRACTION_GRADIENT_PROJECTED );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SUBSCALE_PRESSURE);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DISSIPATIVE_FORCE_COEFF_JM);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DISSIPATIVE_FORCE_COEFF_SM);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DISSIPATIVE_FORCE_COEFF_BM);
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SOLID_LIQIUD_SURFTENS_COEFF);
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SOLID_AIR_SURFTENS_COEFF);


}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
