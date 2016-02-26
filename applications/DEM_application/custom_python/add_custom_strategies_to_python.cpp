//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/explicit_solver_strategy.h"
#include "custom_strategies/strategies/explicit_solver_continuum.h"
#include "custom_strategies/strategies/iterative_solver_strategy.h"
#include "custom_strategies/strategies/verlet_solver_strategy.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/dem_integration_scheme.h"
#include "custom_strategies/schemes/forward_euler_scheme.h"
#include "custom_strategies/schemes/symplectic_euler_scheme.h"
#include "custom_strategies/schemes/newmark_beta_scheme.h"
#include "custom_strategies/schemes/verlet_velocity_scheme.h"
#include "custom_strategies/schemes/constant_average_acceleration_scheme.h"
#include "custom_strategies/schemes/taylor_scheme.h"

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
          typedef UblasSpace<double, CompressedMatrix, Vector>                          SparseSpaceType;
          typedef UblasSpace<double, Matrix, Vector>                                    LocalSpaceType;

          typedef LinearSolver<SparseSpaceType, LocalSpaceType >                        LinearSolverType;
          typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >  BaseSolvingStrategyType;
          //typedef Scheme< SparseSpaceType, LocalSpaceType >                             BaseSchemeType;
          
          typedef OMP_DEMSearch                                                         OmpDemSearchType;
          typedef DEMSearch<OmpDemSearchType >                                          DemSearchType;

          class_<DEMIntegrationScheme, boost::noncopyable >
                    ("DEMIntegrationScheme", init< >())
            .def("AddSpheresVariables",&DEMIntegrationScheme::AddSpheresVariables)
            .def("AddClustersVariables",&DEMIntegrationScheme::AddClustersVariables)
                  ;

          class_< ForwardEulerScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>
          (
                    "ForwardEulerScheme", init<>()
                  )
                  ;
                  
          class_< SymplecticEulerScheme, bases<DEMIntegrationScheme>,  boost::noncopyable> 
          (
                    "SymplecticEulerScheme", init<>()
                  )
                  ;

          class_< NewmarkBetaScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>
          (
                    "NewmarkBetaScheme", init<const double, const double>()
                  )
                  ;

          class_< VerletVelocityScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>
          (
                    "VerletVelocityScheme", init<>()
                  )
                  ;
                  
          class_< TaylorScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>
          (
                    "TaylorScheme", init<>()
                  )
                  ;

          class_< ConstAverageAccelerationScheme, bases<DEMIntegrationScheme>,  boost::noncopyable>
          (
                    "ConstAverageAccelerationScheme", init<>()
                  )
                  ;
                    
          class_<DemSearchType, bases<SpatialSearch>, boost::noncopyable>
                    ("OMP_DEMSearch", init<>())
                    ;
        
          typedef ExplicitSolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitSolverStrategyType;
          typedef ContinuumExplicitSolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType > ContinuumExplicitSolverStrategyType;
          typedef IterativeSolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType > IterativeSolverStrategy;
          
          class_< ExplicitSolverSettings, boost::noncopyable >
          (
            "ExplicitSolverSettings", init<>() )
          .def_readwrite("r_model_part",&ExplicitSolverSettings::r_model_part)
          .def_readwrite("contact_model_part",&ExplicitSolverSettings::contact_model_part)
          .def_readwrite("fem_model_part",&ExplicitSolverSettings::fem_model_part)
          .def_readwrite("cluster_model_part",&ExplicitSolverSettings::cluster_model_part)
          .def_readwrite("inlet_model_part",&ExplicitSolverSettings::inlet_model_part)
          ;
          
          class_< ExplicitSolverStrategyType, bases< BaseSolvingStrategyType >,  boost::noncopyable>
          (
          "ExplicitSolverStrategy", init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, DEMIntegrationScheme::Pointer, SpatialSearch::Pointer>())
                  .def("Initialize", &ExplicitSolverStrategyType::Initialize)
                  .def("InitialTimeStepCalculation", &ExplicitSolverStrategyType::InitialTimeStepCalculation)
                  .def("PrepareElementsForPrinting", &ExplicitSolverStrategyType::PrepareElementsForPrinting)
          ;
          
          class_< ContinuumExplicitSolverStrategyType, bases< ExplicitSolverStrategyType >,  boost::noncopyable>
          (
          "ContinuumExplicitSolverStrategy", init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, DEMIntegrationScheme::Pointer, SpatialSearch::Pointer>())
                  .def("PrepareContactElementsForPrinting", &ContinuumExplicitSolverStrategyType::PrepareContactElementsForPrinting)                  
          ;
         
          class_< IterativeSolverStrategy, bases< ExplicitSolverStrategyType >,  boost::noncopyable>
          (
          "IterativeSolverStrategy", init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, DEMIntegrationScheme::Pointer, SpatialSearch::Pointer>())
          
          ;
         
          class_< VerletVelocitySolverStrategy<ExplicitSolverStrategyType>,  boost::noncopyable>
          (
          "VerletVelocitySolverStrategy", init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, DEMIntegrationScheme::Pointer, SpatialSearch::Pointer>())
                                         
          ;

          class_< VerletVelocitySolverStrategy<ContinuumExplicitSolverStrategyType>,  boost::noncopyable>
          (
          "ContinuumVerletVelocitySolverStrategy", init< ExplicitSolverSettings&, double, double, double, int, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, DEMIntegrationScheme::Pointer, SpatialSearch::Pointer>())

          ;


        }

    }  // namespace Python.

} // Namespace Kratos
