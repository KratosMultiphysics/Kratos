//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

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

// PoromechanicsApplication dependencies
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_quasistatic_damped_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_dynamic_U_Pw_scheme.hpp"

namespace Kratos {
namespace Python {

void AddTrilinosSchemesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using TrilinosSparseSpace = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using UblasLocalSpace = UblasSpace<double, Matrix, Vector>;

    using TrilinosBaseScheme = Scheme< TrilinosSparseSpace, UblasLocalSpace >;

    using TrilinosNewmarkQuasistaticUPwScheme = NewmarkQuasistaticUPwScheme<TrilinosSparseSpace, UblasLocalSpace>;
    py::class_< TrilinosNewmarkQuasistaticUPwScheme, typename TrilinosNewmarkQuasistaticUPwScheme::Pointer, TrilinosBaseScheme >
    (m, "TrilinosNewmarkQuasistaticUPwScheme")
    .def( py::init< double, double, double >() );

    using TrilinosNewmarkQuasistaticDampedUPwScheme = NewmarkQuasistaticDampedUPwScheme<TrilinosSparseSpace, UblasLocalSpace>;
    py::class_< TrilinosNewmarkQuasistaticDampedUPwScheme, typename TrilinosNewmarkQuasistaticDampedUPwScheme::Pointer, TrilinosBaseScheme >
    (m, "TrilinosNewmarkQuasistaticDampedUPwScheme")
    .def( py::init< double, double, double >() );

    using TrilinosNewmarkDynamicUPwScheme = NewmarkDynamicUPwScheme<TrilinosSparseSpace, UblasLocalSpace>;
    py::class_< TrilinosNewmarkDynamicUPwScheme, typename TrilinosNewmarkDynamicUPwScheme::Pointer, TrilinosBaseScheme >
    (m, "TrilinosNewmarkDynamicUPwScheme")
    .def( py::init< double, double, double >() );
}

}  // namespace Python.
} // Namespace Kratos
