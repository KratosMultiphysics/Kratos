//
// Author:
// Guillermo Casas gcasas@cimne.upc.edu
//
#if !defined(KRATOS_ADAMS_BASHFORTH_STRATEGY)
#define KRATOS_ADAMS_BASHFORTH_STRATEGY

// /* External includes */

// System includes

// Project includes
#include "utilities/timer.h"
#include "custom_elements/Particle_Contact_Element.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <time.h>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

//#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "../DEM_application/custom_strategies/strategies/explicit_solver_strategy.h"

#ifdef USING_CGAL
#include <CGAL/spatial_sort.h>
#endif

/* Timer defines */
#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

namespace Kratos {

    class AdamsBashforthStrategy: public ExplicitSolverStrategy {
    public:

        typedef ExplicitSolverStrategy BaseType;
        typedef typename BaseType::NodesArrayType NodesArrayType;
        typedef typename BaseType::ElementsArrayType ElementsArrayType;
        typedef typename BaseType::ElementsIterator ElementsIterator;
        typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
        typedef WeakPointerVector<Element> ParticleWeakVectorType;
        typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
        typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;

        using BaseType::mpInlet_model_part;
        using BaseType::mpCluster_model_part;
        using BaseType::mpContact_model_part;
        using BaseType::GetModelPart;
        using BaseType::GetFemModelPart;
        using BaseType::mNumberOfThreads;
        using BaseType::mListOfSphericParticles;
        using BaseType::mListOfGhostSphericParticles;
        using BaseType::SearchNeighbours;
        using BaseType::SetSearchRadiiOnAllParticles;
        /// Pointer definition of AdamsBashforthStrategy
        KRATOS_CLASS_POINTER_DEFINITION(AdamsBashforthStrategy);

        AdamsBashforthStrategy(){}

        AdamsBashforthStrategy(ExplicitSolverSettings& settings,
                const double max_delta_time,
                const int n_step_search,
                const double safety_factor,
                const int delta_option,
                typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                typename DEM_FEM_Search::Pointer p_dem_fem_search,
                typename DEMIntegrationScheme::Pointer pScheme,
                typename SpatialSearch::Pointer pSpSearch,
                const bool do_search_balls = true):
                ExplicitSolverStrategy(settings,
                                       max_delta_time,
                                       n_step_search,
                                       safety_factor,
                                       delta_option,
                                       p_creator_destructor,
                                       p_dem_fem_search,
                                       pScheme,
                                       pSpSearch,
                                       do_search_balls)
        {
            mFirstStep = true;
            ExplicitSolverStrategy::GetParticleCreatorDestructor() = p_creator_destructor;
        }
        /// Destructor.
        virtual ~AdamsBashforthStrategy() {
            Timer::SetOuputFile("TimesPartialRelease");
            Timer::PrintTimingInformation();
        }

        double Solve();

    protected:

        void ReconstructForces(ModelPart& r_model_part);

    private:
        bool mFirstStep;
    }; // Class AdamsBashforthStrategy
} // namespace Kratos.

#endif // KRATOS_ADAMS_BASHFORTH_STRATEGY  defined
