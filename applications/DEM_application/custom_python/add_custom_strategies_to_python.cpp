//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>
#include "spaces/ublas_space.h"

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

namespace Kratos
{
    namespace Python
    {
        using namespace boost::python;

        void  AddCustomStrategiesToPython()
        {
          //typedef UblasSpace<double, CompressedMatrix, Vector>                          SparseSpaceType;
          //typedef UblasSpace<double, Matrix, Vector>                                    LocalSpaceType;
          //typedef LinearSolver<SparseSpaceType, LocalSpaceType >                        LinearSolverType;
          typedef OMP_DEMSearch                                                         OmpDemSearchType;
          typedef DEMSearch<OmpDemSearchType >                                          DemSearchType;

          class_<DEMIntegrationScheme, boost::noncopyable >
                    ("DEMIntegrationScheme", init< >())
            .def("SetTranslationalIntegrationSchemeInProperties", &DEMIntegrationScheme::SetTranslationalIntegrationSchemeInProperties)
            .def("SetRotationalIntegrationSchemeInProperties", &DEMIntegrationScheme::SetRotationalIntegrationSchemeInProperties)
                  ;
          
          class_<Variable<DEMIntegrationScheme::Pointer>, boost::noncopyable >("DEMIntegrationSchemePointerVariable", no_init)
                    .def(self_ns::str(self)
          );
          
          class_<Variable<DEMIntegrationScheme*>, boost::noncopyable >("DEMIntegrationSchemeRawPointerVariable", no_init)
                    .def(self_ns::str(self)
          );

          class_< ForwardEulerScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>("ForwardEulerScheme", init<>());

          class_< SymplecticEulerScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>("SymplecticEulerScheme", init<>());
          
          class_< TaylorScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>("TaylorScheme", init<>());

          class_< VelocityVerletScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>("VelocityVerletScheme", init<>());

          class_< RungeKuttaScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>("RungeKuttaScheme", init<>());

          class_< QuaternionIntegrationScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>("QuaternionIntegrationScheme", init<>());

          class_<DemSearchType, bases<SpatialSearch>, boost::noncopyable>("OMP_DEMSearch")
                  .def(init<>())
                  .def(init<const double, const double, const double, const double, const double, const double>(args("min_x", "min_y", "min_z", "max_x", "max_y", "max_z")))
                  ;

          class_< ExplicitSolverSettings, boost::noncopyable >("ExplicitSolverSettings", init<>() )
          .def_readwrite("r_model_part",&ExplicitSolverSettings::r_model_part)
          .def_readwrite("contact_model_part",&ExplicitSolverSettings::contact_model_part)
          .def_readwrite("fem_model_part",&ExplicitSolverSettings::fem_model_part)
          .def_readwrite("cluster_model_part",&ExplicitSolverSettings::cluster_model_part)
          .def_readwrite("inlet_model_part",&ExplicitSolverSettings::inlet_model_part)
          ;

          class_< ExplicitSolverStrategy,  boost::noncopyable>
          ("ExplicitSolverStrategy", init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, const bool>())
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

          class_< ContinuumExplicitSolverStrategy, bases< ExplicitSolverStrategy >,  boost::noncopyable>
          (
          "ContinuumExplicitSolverStrategy", init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer,  SpatialSearch::Pointer>())
                  .def("PrepareContactElementsForPrinting", &ContinuumExplicitSolverStrategy::PrepareContactElementsForPrinting)
          ;

          class_< IterativeSolverStrategy, bases< ExplicitSolverStrategy >,  boost::noncopyable>
          (
          "IterativeSolverStrategy", init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer, const bool>())

          ;

          class_< VelocityVerletSolverStrategy<ExplicitSolverStrategy>, bases< ExplicitSolverStrategy >,  boost::noncopyable>
          (
          "VelocityVerletSolverStrategy", init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer>())

          ;

          class_< VelocityVerletSolverStrategy<ContinuumExplicitSolverStrategy>, bases< ContinuumExplicitSolverStrategy >,  boost::noncopyable>
          (
          "ContinuumVelocityVerletSolverStrategy", init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, SpatialSearch::Pointer>())

          ;


        }

    }  // namespace Python.

} // Namespace Kratos
