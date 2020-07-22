/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Thomas Oberbichler
*/

#if !defined(KRATOS_ADD_SOLVERS_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_SOLVERS_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes

//MKL implements some of FEAST's function. This define excludes the definitions of MKL's FEAST
//functions in mkl/include/mkl_solvers_ee.h. Otherwise, some functions are defined multiple times
#if defined USE_EIGEN_FEAST
#define _MKL_SOLVERS_EE_H
#endif

namespace Kratos {
namespace Python {

void AddCustomSolversToPython(pybind11::module& m);

} // namespace Python
} // namespace Kratos

#endif // defined(KRATOS_ADD_SOLVERS_TO_PYTHON_H_INCLUDED)
