//
//   Project Name:        KratosAdjointFluidApplication $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:          February 2015 $
//   Revision:            $Revision:                0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)
// External includes
#include "pybind11/pybind11.h"

// Application includes
#include "includes/define_python.h"
#include "adjoint_fluid_application.h"
#include "custom_python/add_custom_schemes_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(KratosAdjointFluidApplication,m)
{
    class_<KratosAdjointFluidApplication,
           KratosAdjointFluidApplication::Pointer,
           KratosApplication >(m,"KratosAdjointFluidApplication")
           .def(init<>())
           ;

  AddCustomSchemesToPython(m);

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NUMERICAL_DIFFUSION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VMS_STEADY_TERM_PRIMAL_GRADIENT_MATRIX)

}

} // namespace Python

} // namespace Kratos

#endif // KRATOS_PYTHON defined
