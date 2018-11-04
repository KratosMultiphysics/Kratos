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

using namespace pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef OMP_DEMSearch                                                         OmpDemSearchType;
    typedef DEMSearch<OmpDemSearchType >                                          DemSearchType;

    class_<DEMIntegrationScheme, DEMIntegrationScheme::Pointer>(m, "DEMIntegrationScheme")
        .def(init<>())
        .def("SetTranslationalIntegrationSchemeInProperties", &DEMIntegrationScheme::SetTranslationalIntegrationSchemeInProperties)
        .def("SetRotationalIntegrationSchemeInProperties", &DEMIntegrationScheme::SetRotationalIntegrationSchemeInProperties)
        ;

    class_<Variable<DEMIntegrationScheme::Pointer>, Variable<DEMIntegrationScheme::Pointer>::Pointer >(m, "DEMIntegrationSchemePointerVariable")
        .def("__str__", &Variable<DEMIntegrationScheme::Pointer>::Info)
        ;

    class_<Variable<DEMIntegrationScheme*>, Variable<DEMIntegrationScheme*>::Pointer>(m, "DEMIntegrationSchemeRawPointerVariable")
        .def("__str__", &Variable<DEMIntegrationScheme*>::Info)
        ;

    class_<ForwardEulerScheme, ForwardEulerScheme::Pointer, DEMIntegrationScheme>(m, "ForwardEulerScheme")
        .def(init<>())
        ;

    class_<SymplecticEulerScheme, SymplecticEulerScheme::Pointer, DEMIntegrationScheme>(m, "SymplecticEulerScheme")
        .def(init<>())
        ;

    class_<TaylorScheme, TaylorScheme::Pointer, DEMIntegrationScheme>(m, "TaylorScheme")
        .def(init<>())
        ;

    class_<VelocityVerletScheme, VelocityVerletScheme::Pointer, DEMIntegrationScheme>(m, "VelocityVerletScheme")
        .def(init<>())
        ;

    class_<RungeKuttaScheme, RungeKuttaScheme::Pointer, DEMIntegrationScheme>(m, "RungeKuttaScheme")
        .def(init<>())
        ;

    class_<QuaternionIntegrationScheme, QuaternionIntegrationScheme::Pointer, DEMIntegrationScheme>(m, "QuaternionIntegrationScheme")
        .def(init<>())
        ;

    class_<DemSearchType, DemSearchType::Pointer, SpatialSearch>(m, "OMP_DEMSearch")
        .def(init<>())
        .def(init<const double, const double, const double, const double, const double, const double>(), arg("min_x"), arg("min_y"), arg("min_z"), arg("max_x"), arg("max_y"), arg("max_z"))
        ;

    class_<ExplicitSolverSettings, ExplicitSolverSettings::Pointer>(m, "ExplicitSolverSettings")
        .def(init<>())
        .def_readwrite("r_model_part",&ExplicitSolverSettings::r_model_part)
        .def_readwrite("contact_model_part",&ExplicitSolverSettings::contact_model_part)
        .def_readwrite("fem_model_part",&ExplicitSolverSettings::fem_model_part)
        .def_readwrite("cluster_model_part",&ExplicitSolverSettings::cluster_model_part)
        .def_readwrite("inlet_model_part",&ExplicitSolverSettings::inlet_model_part)
        ;

    class_<ExplicitSolverStrategy, ExplicitSolverStrategy::Pointer>(m, "ExplicitSolverStrategy")
        .def(init< ExplicitSolverSettings&, double, int, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters, const bool>())
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

    class_<ContinuumExplicitSolverStrategy, ContinuumExplicitSolverStrategy::Pointer, ExplicitSolverStrategy>(m, "ContinuumExplicitSolverStrategy")
        .def(init< ExplicitSolverSettings&, double, int, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters>())
        .def("PrepareContactElementsForPrinting", &ContinuumExplicitSolverStrategy::PrepareContactElementsForPrinting)
        ;

    class_<IterativeSolverStrategy, IterativeSolverStrategy::Pointer, ExplicitSolverStrategy>(m, "IterativeSolverStrategy")
        .def(init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters, const bool>())
        ;

    class_<VelocityVerletSolverStrategy<ExplicitSolverStrategy>, VelocityVerletSolverStrategy<ExplicitSolverStrategy>::Pointer, ExplicitSolverStrategy>(m, "VelocityVerletSolverStrategy")
        .def(init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters>())
        ;

    class_<VelocityVerletSolverStrategy<ContinuumExplicitSolverStrategy>, VelocityVerletSolverStrategy<ContinuumExplicitSolverStrategy>::Pointer, ExplicitSolverStrategy>(m, "ContinuumVelocityVerletSolverStrategy")
        .def(init<ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters>())
        ;

}

}  // namespace Python.

} // Namespace Kratos
