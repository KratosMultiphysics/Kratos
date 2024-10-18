//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if defined(KRATOS_PYTHON)

// System includes

// External includes

// Project includes
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"
#include "custom_python/add_trilinos_convergence_criterias_to_python.h"

// Application includes
#include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
#include "custom_strategies/convergencecriterias/trilinos_residual_criteria.h"
#include "custom_strategies/convergencecriterias/trilinos_mixed_generic_criteria.h"

namespace Kratos::Python
{

namespace py = pybind11;

void  AddConvergenceCriterias(pybind11::module& m)
{
    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

    //********************************************************************
    //********************************************************************
    // Convergence criteria base class
    using TrilinosConvergenceCriteria = ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType >;

    using TrilinosConvergenceCriteriaPointer = typename ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType > ::Pointer;

    py::class_< TrilinosConvergenceCriteria, TrilinosConvergenceCriteriaPointer > (m,"TrilinosConvergenceCriteria")
    .def(py::init<>())
    .def("SetActualizeRHSFlag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::SetActualizeRHSFlag)
    .def("GetActualizeRHSflag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::GetActualizeRHSflag)
    .def("PreCriteria", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::PreCriteria)
    .def("PostCriteria", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::PostCriteria)
    .def("Initialize", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Initialize)
    .def("InitializeSolutionStep", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::InitializeSolutionStep)
    .def("FinalizeSolutionStep", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::FinalizeSolutionStep)
    .def("Check", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Check)
    .def("SetEchoLevel", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::SetEchoLevel)
    ;

    py::class_< TrilinosDisplacementCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename TrilinosDisplacementCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria>(m,"TrilinosDisplacementCriteria")
            .def(py::init< TrilinosSparseSpaceType::DataType, TrilinosSparseSpaceType::DataType >())
            .def(py::init< Parameters >())
            ;

    py::class_< TrilinosResidualCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename TrilinosResidualCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria >
            (m,"TrilinosResidualCriteria")
            .def(py::init< TrilinosSparseSpaceType::DataType, TrilinosSparseSpaceType::DataType >())
            .def(py::init< Parameters >())
            ;

    py::class_<And_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename And_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria>
            (m,"TrilinosAndCriteria")
            .def(py::init<TrilinosConvergenceCriteriaPointer, TrilinosConvergenceCriteriaPointer > ());

    py::class_<Or_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename Or_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria>
            (m,"TrilinosOrCriteria")
            .def(py::init<TrilinosConvergenceCriteriaPointer, TrilinosConvergenceCriteriaPointer > ());

    using ConvergenceVariableListType = typename TrilinosMixedGenericCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType>::ConvergenceVariableListType;
    py::class_<
        TrilinosMixedGenericCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename TrilinosMixedGenericCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosConvergenceCriteria>(m, "TrilinosMixedGenericCriteria")
        .def(py::init< const ConvergenceVariableListType& >());
}

} // namespace Kratos::Python.

#endif // KRATOS_PYTHON defined
