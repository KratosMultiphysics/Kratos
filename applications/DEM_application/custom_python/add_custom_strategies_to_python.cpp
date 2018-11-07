//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/explicit_solver_strategy.h"
#include "custom_strategies/strategies/explicit_solver_continuum.h"
#include "custom_strategies/strategies/iterative_solver_strategy.h"
#include "custom_strategies/strategies/velocity_verlet_solver_strategy.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/dem_integration_scheme.h"
#include "custom_strategies/schemes/forward_euler_scheme.h"
#include "custom_strategies/schemes/symplectic_euler_scheme.h"
#include "custom_strategies/schemes/taylor_scheme.h"
#include "custom_strategies/schemes/velocity_verlet_scheme.h"
#include "custom_strategies/schemes/runge_kutta_scheme.h"
#include "custom_strategies/schemes/quaternion_integration_scheme.h"

//builder_and_solvers
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

//Search
#include "custom_utilities/omp_dem_search.h"
#include "custom_utilities/dem_fem_search.h"

// Create and Destroy
#include "custom_utilities/create_and_destroy.h"

namespace Kratos {
namespace Python {

namespace py = pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef OMP_DEMSearch                                                         OmpDemSearchType;
    typedef DEMSearch<OmpDemSearchType >                                          DemSearchType;

    py::class_<DEMIntegrationScheme, DEMIntegrationScheme::Pointer>(m, "DEMIntegrationScheme")
        .def(py::init<>())
        .def("SetTranslationalIntegrationSchemeInProperties", &DEMIntegrationScheme::SetTranslationalIntegrationSchemeInProperties)
        .def("SetRotationalIntegrationSchemeInProperties", &DEMIntegrationScheme::SetRotationalIntegrationSchemeInProperties)
        ;

    py::class_<Variable<DEMIntegrationScheme::Pointer>, Variable<DEMIntegrationScheme::Pointer>::Pointer >(m, "DEMIntegrationSchemePointerVariable")
        .def("__str__", &Variable<DEMIntegrationScheme::Pointer>::Info)
        ;

    py::class_<Variable<DEMIntegrationScheme*>, Variable<DEMIntegrationScheme*>::Pointer>(m, "DEMIntegrationSchemeRawPointerVariable")
        .def("__str__", &Variable<DEMIntegrationScheme*>::Info)
        ;

    py::class_<ForwardEulerScheme, ForwardEulerScheme::Pointer, DEMIntegrationScheme>(m, "ForwardEulerScheme")
        .def(py::init<>())
        ;

    py::class_<SymplecticEulerScheme, SymplecticEulerScheme::Pointer, DEMIntegrationScheme>(m, "SymplecticEulerScheme")
        .def(py::init<>())
        ;

    py::class_<TaylorScheme, TaylorScheme::Pointer, DEMIntegrationScheme>(m, "TaylorScheme")
        .def(py::init<>())
        ;

    py::class_<VelocityVerletScheme, VelocityVerletScheme::Pointer, DEMIntegrationScheme>(m, "VelocityVerletScheme")
        .def(py::init<>())
        ;

    py::class_<RungeKuttaScheme, RungeKuttaScheme::Pointer, DEMIntegrationScheme>(m, "RungeKuttaScheme")
        .def(py::init<>())
        ;

    py::class_<QuaternionIntegrationScheme, QuaternionIntegrationScheme::Pointer, DEMIntegrationScheme>(m, "QuaternionIntegrationScheme")
        .def(py::init<>())
        ;

    py::class_<DemSearchType, DemSearchType::Pointer, SpatialSearch>(m, "OMP_DEMSearch")
        .def(py::init<>())
        .def(py::init<const double, const double, const double, const double, const double, const double>(), py::arg("min_x"), py::arg("min_y"), py::arg("min_z"), py::arg("max_x"), py::arg("max_y"), py::arg("max_z"))
        ;

    py::class_<ExplicitSolverSettings, ExplicitSolverSettings::Pointer>(m, "ExplicitSolverSettings")
        .def(py::init<>())
        .def_readwrite("r_model_part",&ExplicitSolverSettings::r_model_part)
        .def_readwrite("contact_model_part",&ExplicitSolverSettings::contact_model_part)
        .def_readwrite("fem_model_part",&ExplicitSolverSettings::fem_model_part)
        .def_readwrite("cluster_model_part",&ExplicitSolverSettings::cluster_model_part)
        .def_readwrite("inlet_model_part",&ExplicitSolverSettings::inlet_model_part)
        ;

    py::class_<ExplicitSolverStrategy, ExplicitSolverStrategy::Pointer>(m, "ExplicitSolverStrategy")
        .def(py::init< ExplicitSolverSettings&, double, int, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters, const bool>())
        .def("Solve", &ExplicitSolverStrategy::Solve)
        .def("Initialize", &ExplicitSolverStrategy::Initialize)
        .def("SetSearchRadiiOnAllParticles", &ExplicitSolverStrategy::SetSearchRadiiOnAllParticles)
        .def("SetNormalRadiiOnAllParticles", &ExplicitSolverStrategy::SetNormalRadiiOnAllParticles)
        .def("SetSearchRadiiWithFemOnAllParticles", &ExplicitSolverStrategy::SetSearchRadiiWithFemOnAllParticles)
        .def("RebuildListOfDiscontinuumSphericParticles", &ExplicitSolverStrategy::RebuildListOfDiscontinuumSphericParticles)
        .def("InitialTimeStepCalculation", &ExplicitSolverStrategy::InitialTimeStepCalculation)
        .def("PrepareElementsForPrinting", &ExplicitSolverStrategy::PrepareElementsForPrinting)
        .def("ResetPrescribedMotionFlagsRespectingImposedDofs", &ExplicitSolverStrategy::ResetPrescribedMotionFlagsRespectingImposedDofs)
        ;

    py::class_<ContinuumExplicitSolverStrategy, ContinuumExplicitSolverStrategy::Pointer, ExplicitSolverStrategy>(m, "ContinuumExplicitSolverStrategy")
        .def(py::init< ExplicitSolverSettings&, double, int, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters>())
        .def("PrepareContactElementsForPrinting", &ContinuumExplicitSolverStrategy::PrepareContactElementsForPrinting)
        ;

    py::class_<IterativeSolverStrategy, IterativeSolverStrategy::Pointer, ExplicitSolverStrategy>(m, "IterativeSolverStrategy")
        .def(py::init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters, const bool>())
        ;

    py::class_<VelocityVerletSolverStrategy<ExplicitSolverStrategy>, VelocityVerletSolverStrategy<ExplicitSolverStrategy>::Pointer, ExplicitSolverStrategy>(m, "VelocityVerletSolverStrategy")
        .def(py::init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters>())
        ;

    py::class_<VelocityVerletSolverStrategy<ContinuumExplicitSolverStrategy>, VelocityVerletSolverStrategy<ContinuumExplicitSolverStrategy>::Pointer, ExplicitSolverStrategy>(m, "ContinuumVelocityVerletSolverStrategy")
        .def(py::init<ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters>())
        ;

}

}  // namespace Python.

} // Namespace Kratos
