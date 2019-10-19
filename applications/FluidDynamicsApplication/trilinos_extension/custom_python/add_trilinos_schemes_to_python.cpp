//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "add_trilinos_utilities_to_python.h"

// Trilinos includes
#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_MpiComm.h"

// KratosCore dependencies
#include "containers/variable.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "spaces/ublas_space.h"

// TrilinosApplication dependencies
#include "trilinos_space.h"

// FluidDynamicsApplication dependencies
#include "custom_strategies/strategies/gear_scheme.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bdf_scheme_turbulent.h"
#include "custom_strategies/strategies/residualbased_simple_steady_scheme.h"

namespace Kratos {
namespace Python {

void AddTrilinosSchemesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using TrilinosSparseSpace = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using UblasLocalSpace = UblasSpace<double, Matrix, Vector>;

    using TrilinosBaseScheme = Scheme< TrilinosSparseSpace, UblasLocalSpace >;

    using TrilinosGearScheme = GearScheme<TrilinosSparseSpace, UblasLocalSpace>;
    py::class_< TrilinosGearScheme, typename TrilinosGearScheme::Pointer, TrilinosBaseScheme >( m,"TrilinosGearScheme")
    .def(py::init<Process::Pointer >() )
    .def(py::init<>()) // constructor without a turbulence model
    .def(py::init<const Variable<int>&>()) // constructor for periodic conditions
    ;

    using TrilinosVelocityBossakSchemeTurbulent = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TrilinosSparseSpace, UblasLocalSpace>;
    py::class_ < TrilinosVelocityBossakSchemeTurbulent, typename TrilinosVelocityBossakSchemeTurbulent::Pointer,TrilinosBaseScheme >
    (m,"TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent")
    .def(py::init<double, double, unsigned int, Process::Pointer >())
    .def(py::init<double, double, unsigned int, double, Process::Pointer >())
    .def(py::init<double,double,unsigned int >())
    .def(py::init<double,unsigned int, const Variable<int>&>())
    ;

    using TrilinosVelocityBDFSchemeTurbulent = ResidualBasedPredictorCorrectorBDFSchemeTurbulent<TrilinosSparseSpace, UblasLocalSpace>;
    py::class_ < TrilinosVelocityBDFSchemeTurbulent, typename TrilinosVelocityBDFSchemeTurbulent::Pointer, TrilinosBaseScheme >
    (m,"TrilinosResidualBasedPredictorCorrectorBDFScheme")
    .def(py::init<unsigned int, Kratos::Flags& >() )
    ;

    using TrilinosResidualBasedSimpleSteadyScheme = ResidualBasedSimpleSteadyScheme<TrilinosSparseSpace, UblasLocalSpace>;
    py::class_ < TrilinosResidualBasedSimpleSteadyScheme, typename TrilinosResidualBasedSimpleSteadyScheme::Pointer, TrilinosBaseScheme >
    (m,"TrilinosResidualBasedSimpleSteadyScheme")
    .def(py::init<double, double, unsigned int, Process::Pointer >())
    .def(py::init<double,double,unsigned int >())
    ;

}

}
}
