// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "custom_python/add_custom_schemes_to_python.h"
#include "custom_schemes/adjoint_bossak_scheme.h"
#include "custom_schemes/adjoint_steady_scheme.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

void AddCustomSchemesToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;

    class_<AdjointBossakScheme<SparseSpaceType, LocalSpaceType>, bases<SchemeType>, boost::noncopyable>(
        "AdjointBossakScheme", init<Parameters&, ObjectiveFunction::Pointer>())
        ;

    class_<AdjointSteadyScheme<SparseSpaceType, LocalSpaceType>, bases<SchemeType>, boost::noncopyable>(
        "AdjointSteadyScheme", init<Parameters&, ObjectiveFunction::Pointer>())
        ;
}

} // namespace Python

} // namespace Kratos
