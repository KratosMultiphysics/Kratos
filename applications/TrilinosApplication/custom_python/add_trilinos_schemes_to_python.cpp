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
#include "includes/define_python.h"
#include "custom_python/add_trilinos_schemes_to_python.h"
#include "trilinos_space.h"
#ifdef HAVE_TPETRA
#include "trilinos_space_experimental.h"
#endif
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

/* Schemes */
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_bdf_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_custom_scheme.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "solving_strategies/schemes/residual_based_adjoint_steady_scheme.h"
#include "solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"

/* Response function */
#include "response_functions/adjoint_response_function.h"

namespace Kratos::Python
{
namespace py = pybind11;

using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

namespace {
    template<class TSparseSpace, class TLocalSpace>
    void RegisterSchemes(pybind11::module& m, const std::string& Prefix)
    {
        using BaseSchemeType = Scheme<TSparseSpace, TLocalSpace>;

        py::class_<BaseSchemeType, typename BaseSchemeType::Pointer>(m, (Prefix + "Scheme").c_str())
            .def(py::init<>())
            .def("Initialize",                  &BaseSchemeType::Initialize)
            .def("SchemeIsInitialized",          &BaseSchemeType::SchemeIsInitialized)
            .def("ElementsAreInitialized",       &BaseSchemeType::ElementsAreInitialized)
            .def("ConditionsAreInitialized",     &BaseSchemeType::ConditionsAreInitialized)
            .def("InitializeElements",           &BaseSchemeType::InitializeElements)
            .def("InitializeConditions",         &BaseSchemeType::InitializeConditions)
            .def("InitializeSolutionStep",       &BaseSchemeType::InitializeSolutionStep)
            .def("FinalizeSolutionStep",         &BaseSchemeType::FinalizeSolutionStep)
            .def("InitializeNonLinIteration",    &BaseSchemeType::InitializeNonLinIteration)
            .def("FinalizeNonLinIteration",      &BaseSchemeType::FinalizeNonLinIteration)
            .def("Predict",                      &BaseSchemeType::Predict)
            .def("Update",                       &BaseSchemeType::Update)
            .def("CalculateOutputData",          &BaseSchemeType::CalculateOutputData)
            .def("Clean",                        &BaseSchemeType::Clean)
            .def("Clear",                        &BaseSchemeType::Clear)
            .def("MoveMesh", [](BaseSchemeType& self, ModelPart::NodesContainerType& rNodes) {
                for (auto it = rNodes.begin(); it != rNodes.end(); ++it) {
                    const array_1d<double, 3>& disp = it->FastGetSolutionStepValue(DISPLACEMENT);
                    it->X() = it->X0() + disp[0];
                    it->Y() = it->Y0() + disp[1];
                    it->Z() = it->Z0() + disp[2];
                }
            })
            .def("Check", [](const BaseSchemeType& self, const ModelPart& rModelPart) { return self.Check(rModelPart); })
            ;

        py::class_<
            ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TLocalSpace>,
            typename ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TLocalSpace>::Pointer,
            BaseSchemeType>(m, (Prefix + "ResidualBasedIncrementalUpdateStaticScheme").c_str())
            .def(py::init<>())
            ;

        py::class_<
            ResidualBasedIncrementalUpdateStaticSchemeSlip<TSparseSpace, TLocalSpace>,
            typename ResidualBasedIncrementalUpdateStaticSchemeSlip<TSparseSpace, TLocalSpace>::Pointer,
            ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TLocalSpace>>(m, (Prefix + "ResidualBasedIncrementalUpdateStaticSchemeSlip").c_str())
            .def(py::init<unsigned int, unsigned int>())
            ;

        py::class_<
            ResidualBasedBossakDisplacementScheme<TSparseSpace, TLocalSpace>,
            typename ResidualBasedBossakDisplacementScheme<TSparseSpace, TLocalSpace>::Pointer,
            BaseSchemeType>(m, (Prefix + "ResidualBasedBossakDisplacementScheme").c_str())
            .def(py::init<double>())
            ;

        py::class_<
            ResidualBasedBDFDisplacementScheme<TSparseSpace, TLocalSpace>,
            typename ResidualBasedBDFDisplacementScheme<TSparseSpace, TLocalSpace>::Pointer,
            BaseSchemeType>(m, (Prefix + "ResidualBasedBDFDisplacementScheme").c_str())
            .def(py::init<>())
            .def(py::init<const std::size_t>())
            ;

        py::class_<
            ResidualBasedBDFCustomScheme<TSparseSpace, TLocalSpace>,
            typename ResidualBasedBDFCustomScheme<TSparseSpace, TLocalSpace>::Pointer,
            BaseSchemeType>(m, (Prefix + "ResidualBasedBDFCustomScheme").c_str())
            .def(py::init<>())
            .def(py::init<const std::size_t>())
            .def(py::init<const std::size_t, Parameters>())
            ;

        using AdjointStaticSchemeType = ResidualBasedAdjointStaticScheme<TSparseSpace, TLocalSpace>;
        py::class_<AdjointStaticSchemeType, typename AdjointStaticSchemeType::Pointer, BaseSchemeType>
            (m, (Prefix + "ResidualBasedAdjointStaticScheme").c_str())
            .def(py::init<AdjointResponseFunction::Pointer>())
            .def("SetResponseFunction", &AdjointStaticSchemeType::SetResponseFunction, py::arg("new_response_function"))
            ;

        using AdjointSteadySchemeType = ResidualBasedAdjointSteadyScheme<TSparseSpace, TLocalSpace>;
        py::class_<AdjointSteadySchemeType, typename AdjointSteadySchemeType::Pointer, AdjointStaticSchemeType>
            (m, (Prefix + "ResidualBasedAdjointSteadyScheme").c_str())
            .def(py::init<AdjointResponseFunction::Pointer>())
            ;

        using AdjointBossakSchemeType = ResidualBasedAdjointBossakScheme<TSparseSpace, TLocalSpace>;
        py::class_<AdjointBossakSchemeType, typename AdjointBossakSchemeType::Pointer, BaseSchemeType>
            (m, (Prefix + "ResidualBasedAdjointBossakScheme").c_str())
            .def(py::init<Kratos::Parameters, AdjointResponseFunction::Pointer>())
            ;
    }
}

void  AddSchemes(pybind11::module& m)
{
    RegisterSchemes<TrilinosSparseSpaceType, TrilinosLocalSpaceType>(m, "Trilinos");

#ifdef HAVE_TPETRA
    using TrilinosExperimentalSparseSpaceType = TrilinosSpaceExperimental<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>;
    RegisterSchemes<TrilinosExperimentalSparseSpaceType, TrilinosLocalSpaceType>(m, "TrilinosExperimental");
#endif
}

} // namespace Kratos::Python.

#endif // KRATOS_PYTHON defined
