// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "custom_python/add_custom_schemes_to_python.h"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

void AddCustomSchemesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;

}

} // namespace Python

} // namespace Kratos
