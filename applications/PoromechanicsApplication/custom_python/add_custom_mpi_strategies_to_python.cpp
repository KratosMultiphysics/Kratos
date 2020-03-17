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

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"

#include "custom_python/add_custom_mpi_strategies_to_python.h"

// Trilinos includes
#include "Epetra_FEVector.h"

// Project includes
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

//linear solvers

//strategies

//builders and solvers

//schemes
#include "custom_strategies/schemes/trilinos_newmark_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/trilinos_newmark_quasistatic_damped_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/trilinos_newmark_dynamic_U_Pw_scheme.hpp"


namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void  AddCustomMPIStrategiesToPython(pybind11::module& m)
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
    typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;

    typedef TrilinosNewmarkQuasistaticUPwScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosNewmarkQuasistaticUPwSchemeType;
    typedef TrilinosNewmarkQuasistaticDampedUPwScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosNewmarkQuasistaticDampedUPwSchemeType;
    typedef TrilinosNewmarkDynamicUPwScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosNewmarkDynamicUPwSchemeType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Schemes
    py::class_< TrilinosNewmarkQuasistaticUPwSchemeType, typename TrilinosNewmarkQuasistaticUPwSchemeType::Pointer, TrilinosBaseSchemeType >
    (m, "TrilinosNewmarkQuasistaticUPwScheme")
    .def( py::init< double, double, double >() );
    py::class_< TrilinosNewmarkQuasistaticDampedUPwSchemeType, typename TrilinosNewmarkQuasistaticDampedUPwSchemeType::Pointer, TrilinosBaseSchemeType >
    (m, "TrilinosNewmarkQuasistaticDampedUPwScheme")
    .def( py::init<  double, double, double >());
    py::class_< TrilinosNewmarkDynamicUPwSchemeType, typename TrilinosNewmarkDynamicUPwSchemeType::Pointer, TrilinosBaseSchemeType >
    (m, "TrilinosNewmarkDynamicUPwScheme")
    .def( py::init<  double, double, double >());

}

}  // namespace Python.
} // Namespace Kratos

#endif // KRATOS_PYTHON defined