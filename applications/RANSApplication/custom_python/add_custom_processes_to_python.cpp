// System includes

#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif                                    // KRATOS_USE_AMATRIX

// External includes
#include "pybind11/pybind11.h"

// Project includes

// Application includes
#include "custom_python/add_custom_auxiliary_processes_to_python.h"
#include "custom_python/add_custom_solving_processes_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomProcessesToPython(pybind11::module& m)
{
    AddCustomSolvingProcessesToPython(m);
    AddCustomAuxiliaryProcessesToPython(m);
}

} // namespace Python.
} // Namespace Kratos
