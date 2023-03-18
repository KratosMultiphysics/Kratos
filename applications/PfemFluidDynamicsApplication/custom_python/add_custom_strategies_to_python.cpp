//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

// strategies
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "custom_strategies/strategies/v_p_strategy.h"
#include "custom_strategies/strategies/two_step_v_p_strategy.h"
#include "custom_strategies/strategies/two_step_v_p_thermal_strategy.h"
#include "custom_strategies/strategies/three_step_v_p_strategy.h"
#include "custom_strategies/strategies/gauss_seidel_linear_strategy.h"
#include "custom_strategies/strategies/nodal_two_step_v_p_strategy.h"
#include "custom_strategies/strategies/nodal_two_step_v_p_strategy_for_FSI.h"
#include "custom_strategies/strategies/two_step_v_p_DEM_coupling_strategy.h"

// builder_and_solvers

// convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// linear solvers
#include "linear_solvers/linear_solver.h"

// schemes

namespace Kratos
{

    namespace Python
    {
        namespace py = pybind11;

        void AddCustomStrategiesToPython(pybind11::module &m)
        {
            typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
            typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

            // base types
            typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
            typedef SolvingStrategy<SparseSpaceType, LocalSpaceType> BaseSolvingStrategyType;
            typedef ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> ImplicitBaseSolvingStrategyType;
            typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
            typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;
            // typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;

            // custom strategy types
            typedef ThreeStepVPStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> ThreeStepVPStrategyType;
            typedef VPStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> VPStrategyType;
            typedef TwoStepVPStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> TwoStepVPStrategyType;
            typedef TwoStepVPThermalStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> TwoStepVPThermalStrategyType;
            typedef TwoStepVPDEMcouplingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> TwoStepVPDEMcouplingStrategyType;
            typedef NodalTwoStepVPStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> NodalTwoStepVPStrategyType;
            typedef NodalTwoStepVPStrategyForFSI<SparseSpaceType, LocalSpaceType, LinearSolverType> NodalTwoStepVPStrategyForFSIType;
            typedef GaussSeidelLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> GaussSeidelLinearStrategyType;

            // Clustom Settings Types
            typedef TwoStepVPSolverSettings<SparseSpaceType, LocalSpaceType, LinearSolverType> TwoStepVPSolverSettingsType;
            //********************************************************************
            //*************************SHCHEME CLASSES****************************
            //********************************************************************

            py::class_<VPStrategyType, VPStrategyType::Pointer, BaseSolvingStrategyType>(m, "VPStrategy")
                .def(py::init<ModelPart &, TwoStepVPSolverSettingsType &>())
                .def("CalculateAccelerations", &VPStrategyType::CalculateAccelerations)
                .def("CalculateDisplacements", &VPStrategyType::CalculateDisplacementsAndPorosity)
                // .def("InitializeStressStrain",&VPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::InitializeStressStrain)
                ;

            py::class_<TwoStepVPStrategyType, TwoStepVPStrategyType::Pointer, VPStrategyType>(m, "TwoStepVPStrategy")
                .def(py::init<ModelPart &, LinearSolverType::Pointer, LinearSolverType::Pointer, bool, double, double, int, unsigned int, unsigned int>());
            ;

            py::class_<TwoStepVPThermalStrategyType, TwoStepVPThermalStrategyType::Pointer, TwoStepVPStrategyType>(m, "TwoStepVPThermalStrategy")
                .def(py::init<ModelPart &, LinearSolverType::Pointer, LinearSolverType::Pointer, bool, double, double, int, unsigned int, unsigned int>());
            ;

            py::class_<ThreeStepVPStrategyType, ThreeStepVPStrategyType::Pointer, VPStrategyType>(m, "ThreeStepVPStrategy")
                .def(py::init<ModelPart &, LinearSolverType::Pointer, LinearSolverType::Pointer, bool, double, double, int, unsigned int, unsigned int>());
            ;

            py::class_<TwoStepVPDEMcouplingStrategyType, TwoStepVPDEMcouplingStrategyType::Pointer, TwoStepVPStrategyType>(m, "TwoStepVPDEMcouplingStrategy")
                .def(py::init<ModelPart &, LinearSolverType::Pointer, LinearSolverType::Pointer, bool, double, double, int, unsigned int, unsigned int>());

            py::class_<NodalTwoStepVPStrategyType, NodalTwoStepVPStrategyType::Pointer, ImplicitBaseSolvingStrategyType>(m, "NodalTwoStepVPStrategy")
                .def(py::init<ModelPart &, LinearSolverType::Pointer, LinearSolverType::Pointer, bool, double, double, int, unsigned int, unsigned int>());

            py::class_<NodalTwoStepVPStrategyForFSIType, NodalTwoStepVPStrategyForFSIType::Pointer, NodalTwoStepVPStrategyType>(m, "NodalTwoStepVPStrategyForFSI")
                .def(py::init<ModelPart &, LinearSolverType::Pointer, LinearSolverType::Pointer, bool, double, double, int, unsigned int, unsigned int>());

            py::class_<GaussSeidelLinearStrategyType, GaussSeidelLinearStrategyType::Pointer, ImplicitBaseSolvingStrategyType>(m, "GaussSeidelLinearStrategy")
                .def(py::init<ModelPart &, BaseSchemeType::Pointer, LinearSolverType::Pointer, BuilderAndSolverType::Pointer, bool, bool>())
                .def("GetResidualNorm", &GaussSeidelLinearStrategyType::GetResidualNorm)
                .def("SetBuilderAndSolver", &GaussSeidelLinearStrategyType::SetBuilderAndSolver);
        }

    } // namespace Python.

} // Namespace Kratos
