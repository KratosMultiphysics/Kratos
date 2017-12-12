/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Thomas Oberbichler
*/

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "eigen_solvers_application.h"
#include "custom_python/add_custom_solvers_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

BOOST_PYTHON_MODULE(KratosEigenSolversApplication)
{

	class_<KratosEigenSolversApplication,
		   KratosEigenSolversApplication::Pointer,
		   bases<KratosApplication>, boost::noncopyable>("KratosEigenSolversApplication");

	AddCustomSolversToPython();
}

} // namespace Python

} // namespace Kratos

#endif // defined(KRATOS_PYTHON)
