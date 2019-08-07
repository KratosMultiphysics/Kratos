//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:      JMCarbonell  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

#include "umat_application.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;



PYBIND11_MODULE(KratosUmatApplication,m)
{

  py::class_<KratosUmatApplication,
         KratosUmatApplication::Pointer,
         KratosApplication>(m,"KratosUmatApplication")
      .def(py::init<>())
      ;

  AddCustomConstitutiveLawsToPython(m);

  //registering variables in python
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ALPHA )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BETA )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MF )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CC )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MM )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, KSIS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RHOM )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PC0 )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VOID_RATIO )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PLASTIC_MULTIPLIER )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, P_ATM  )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, E0 )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, LAMBDA_SANISAND )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, EPSI )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, MC )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, ME )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, M_SANISAND )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, G0_SANISAND )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NU_SANISAND )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, H0 )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, CH )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NB )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, A_SANISAND )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, ND )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, PTMULT )

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
