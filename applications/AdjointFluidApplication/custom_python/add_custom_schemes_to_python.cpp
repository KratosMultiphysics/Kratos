// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "custom_python/add_custom_schemes_to_python.h"
#include "custom_schemes/adjoint_bossak_scheme.h"
#include "custom_schemes/adjoint_steady_velocity_pressure_scheme.h"

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

    class_<
        AdjointBossakScheme<SparseSpaceType, LocalSpaceType>,
        typename AdjointBossakScheme<SparseSpaceType, LocalSpaceType>::Pointer,
        SchemeType>(m,"AdjointBossakScheme")
        .def(init<Parameters&, ResponseFunction::Pointer>())
        ;

    class_<
        AdjointSteadyVelocityPressureScheme<SparseSpaceType, LocalSpaceType>,
        typename AdjointSteadyVelocityPressureScheme<SparseSpaceType, LocalSpaceType>::Pointer,
        SchemeType>(m,"AdjointSteadyVelocityPressureScheme")
        .def(init<Parameters&, ResponseFunction::Pointer>())
        ;
}

} // namespace Python

} // namespace Kratos
