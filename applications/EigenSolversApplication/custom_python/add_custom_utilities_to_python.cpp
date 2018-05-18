// ==============================================================================
//  KratosEigenSolversApplication
//
//  License:         BSD License
//
//  Main authors:    Long Chen, https://github.com/longchentum
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// // ------------------------------------------------------------------------------
// // External includes
// // ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
//#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/singular_value_decomposition.h"
// ==============================================================================

namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    // ================================================================
    // For perfoming singular value decomposition
    // ================================================================
    class_<SingularValueDecomposition >(m, "SingularValueDecomposition")
        .def(init<>())    
        .def("Solve", &SingularValueDecomposition::Solve)
        ;
}


}  // namespace Python.

} // Namespace Kratos

