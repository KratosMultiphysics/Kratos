//
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Manuel Messmer
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "custom_python/add_custom_solvers_to_python.h"
#include "custom_solvers/additive_schwarz_preconditioner.h"
#include "factories/standard_preconditioner_factory.h"

#include "iga_application_variables.h"


namespace Kratos {
namespace Python {

void AddCustomSolversToPython(
    pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;
    typedef AdditiveSchwarzPreconditioner<SpaceType,  LocalSpaceType> AdditiveSchwarzPreconditionerType;

    py::class_<AdditiveSchwarzPreconditionerType, AdditiveSchwarzPreconditionerType::Pointer, PreconditionerType>(m,"AdditiveSchwarzPreconditioner")
        .def(py::init<>() )
        .def("__str__", PrintObject<AdditiveSchwarzPreconditionerType>)
    ;

    static auto AdditiveSchwarzPreconditionerFactory = StandardPreconditionerFactory<SpaceType,LocalSpaceType,AdditiveSchwarzPreconditionerType>();
    KRATOS_REGISTER_PRECONDITIONER("additive_schwarz", AdditiveSchwarzPreconditionerFactory);
}

} // namespace Python
} // Namespace Kratos
