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
//                   Vicente Mataix Ferrandiz
//

#if defined(KRATOS_PYTHON)

// System includes

// External includes

// Project includes
#include "trilinos_space.h"
#ifdef HAVE_TPETRA
#include "trilinos_space_experimental.h"
#endif
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

namespace {
    template<class TSparseSpace, class TLocalSpace>
    void RegisterConvergenceCriterias(pybind11::module& m, const std::string& Prefix)
    {
        using ConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TLocalSpace>;
        using ConvergenceCriteriaPointerType = typename ConvergenceCriteria<TSparseSpace, TLocalSpace>::Pointer;

        py::class_<ConvergenceCriteriaType, ConvergenceCriteriaPointerType>(m, (Prefix + "ConvergenceCriteria").c_str())
            .def(py::init<>())
            .def("SetActualizeRHSFlag",      &ConvergenceCriteriaType::SetActualizeRHSFlag)
            .def("GetActualizeRHSflag",      &ConvergenceCriteriaType::GetActualizeRHSflag)
            .def("PreCriteria",              &ConvergenceCriteriaType::PreCriteria)
            .def("PostCriteria",             &ConvergenceCriteriaType::PostCriteria)
            .def("Initialize",               &ConvergenceCriteriaType::Initialize)
            .def("InitializeSolutionStep",   &ConvergenceCriteriaType::InitializeSolutionStep)
            .def("FinalizeSolutionStep",     &ConvergenceCriteriaType::FinalizeSolutionStep)
            .def("Check",                    &ConvergenceCriteriaType::Check)
            .def("SetEchoLevel",             &ConvergenceCriteriaType::SetEchoLevel)
            ;

        py::class_<
            TrilinosDisplacementCriteria<TSparseSpace, TLocalSpace>,
            typename TrilinosDisplacementCriteria<TSparseSpace, TLocalSpace>::Pointer,
            ConvergenceCriteriaType>(m, (Prefix + "DisplacementCriteria").c_str())
            .def(py::init<typename TSparseSpace::DataType, typename TSparseSpace::DataType>())
            .def(py::init<Parameters>())
            ;

        py::class_<
            TrilinosResidualCriteria<TSparseSpace, TLocalSpace>,
            typename TrilinosResidualCriteria<TSparseSpace, TLocalSpace>::Pointer,
            ConvergenceCriteriaType>(m, (Prefix + "ResidualCriteria").c_str())
            .def(py::init<typename TSparseSpace::DataType, typename TSparseSpace::DataType>())
            .def(py::init<Parameters>())
            ;

        py::class_<
            And_Criteria<TSparseSpace, TLocalSpace>,
            typename And_Criteria<TSparseSpace, TLocalSpace>::Pointer,
            ConvergenceCriteriaType>(m, (Prefix + "AndCriteria").c_str())
            .def(py::init<ConvergenceCriteriaPointerType, ConvergenceCriteriaPointerType>())
            ;

        py::class_<
            Or_Criteria<TSparseSpace, TLocalSpace>,
            typename Or_Criteria<TSparseSpace, TLocalSpace>::Pointer,
            ConvergenceCriteriaType>(m, (Prefix + "OrCriteria").c_str())
            .def(py::init<ConvergenceCriteriaPointerType, ConvergenceCriteriaPointerType>())
            ;

        using ConvergenceVariableListType = typename TrilinosMixedGenericCriteria<TSparseSpace, TLocalSpace>::ConvergenceVariableListType;
        py::class_<
            TrilinosMixedGenericCriteria<TSparseSpace, TLocalSpace>,
            typename TrilinosMixedGenericCriteria<TSparseSpace, TLocalSpace>::Pointer,
            ConvergenceCriteriaType>(m, (Prefix + "MixedGenericCriteria").c_str())
            .def(py::init<const ConvergenceVariableListType&>())
            ;
    }
}

void  AddConvergenceCriterias(pybind11::module& m)
{
    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

    RegisterConvergenceCriterias<TrilinosSparseSpaceType, TrilinosLocalSpaceType>(m, "Trilinos");

#ifdef HAVE_TPETRA
    using TrilinosExperimentalSparseSpaceType = TrilinosSpaceExperimental<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>;
    RegisterConvergenceCriterias<TrilinosExperimentalSparseSpaceType, TrilinosLocalSpaceType>(m, "TrilinosExperimental");
#endif
}

} // namespace Kratos::Python.

#endif // KRATOS_PYTHON defined
