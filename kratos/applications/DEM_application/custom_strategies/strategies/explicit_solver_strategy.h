//
// Authors:
// Miguel Angel Celigueta maceli@cimne.upc.edu
// Miquel Santasusana msantasusana@cimne.upc.edu
//


#if !defined(KRATOS_EXPLICIT_SOLVER_STRATEGY)
#define KRATOS_EXPLICIT_SOLVER_STRATEGY


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

#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/dem_integration_scheme.h"
#include "custom_utilities/create_and_destroy.h"
#include "custom_utilities/dem_fem_utilities.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/inlet.h"

#include "custom_elements/cluster3D.h"

////Cfeng
#include "custom_utilities/dem_fem_search.h"
#include "custom_utilities/discrete_particle_configure.h"
#include "custom_utilities/rigid_face_geometrical_object_configure.h"

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

    class ExplicitSolverSettings {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverSettings);

        ExplicitSolverSettings() {
        }

        ~ExplicitSolverSettings() {
        }
        ModelPart::Pointer r_model_part;
        ModelPart::Pointer contact_model_part;
        ModelPart::Pointer fem_model_part;
        ModelPart::Pointer cluster_model_part;
        ModelPart::Pointer inlet_model_part;
    };
    
    class KRATOS_API(DEM_APPLICATION) ExplicitSolverStrategy {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;
        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ElementsArrayType::iterator ElementsIterator;
        typedef ModelPart::ConditionsContainerType ConditionsArrayType;
        typedef ModelPart::NodesContainerType::ContainerType NodesContainerType;
        typedef ModelPart::ElementsContainerType::ContainerType ElementsContainerType;
        typedef ModelPart::ConditionsContainerType::ContainerType ConditionsContainerType;
        typedef SpatialSearch::ResultElementsContainerType ResultElementsContainerType;
        typedef SpatialSearch::VectorResultElementsContainerType VectorResultElementsContainerType;
        typedef SpatialSearch::RadiusArrayType RadiusArrayType;
        typedef SpatialSearch::DistanceType DistanceType;
        typedef SpatialSearch::VectorDistanceType VectorDistanceType;
        typedef SpatialSearch::ResultConditionsContainerType ResultConditionsContainerType;
        typedef SpatialSearch::VectorResultConditionsContainerType VectorResultConditionsContainerType;
        typedef PointerVectorSet<Properties, IndexedObject> PropertiesContainerType;
        typedef typename PropertiesContainerType::iterator PropertiesIterator;
        typedef DiscreteParticleConfigure<3> ElementConfigureType;
        typedef RigidFaceGeometricalObjectConfigure<3> RigidFaceGeometricalConfigureType;
        typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3ul> > > ComponentOf3ComponentsVariableType;

        /// Pointer definition of ExplicitSolverStrategy
        KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverStrategy);

        ExplicitSolverStrategy() {
        }

        ExplicitSolverStrategy(ExplicitSolverSettings& settings,
                const double max_delta_time,
                const int n_step_search,
                const double safety_factor,
                const int delta_option,
                typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                typename DEM_FEM_Search::Pointer p_dem_fem_search,
                typename DEMIntegrationScheme::Pointer pScheme,
                typename SpatialSearch::Pointer pSpSearch,
                const bool do_search_balls = true)
        /*:
        BaseType(*(settings.r_model_part), true)*/ {
            mDeltaOption = delta_option;
            mpParticleCreatorDestructor = p_creator_destructor;
            mpDemFemSearch = p_dem_fem_search;
            mpScheme = pScheme;
            mpSpSearch = pSpSearch;
            mDoSearchNeighbourElements = do_search_balls;
            p_creator_destructor->SetDoSearchNeighbourElements(mDoSearchNeighbourElements);
            mMaxTimeStep = max_delta_time;
            mNStepSearch = n_step_search;
            mSafetyFactor = safety_factor;
            mpDem_model_part = &(*(settings.r_model_part));
            
            if (mpDem_model_part == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "Undefined settings.r_model_part in ExplicitSolverStrategy constructor", "")

                mpContact_model_part = &(*(settings.contact_model_part));
            if (mpContact_model_part == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "Undefined settings.contact_model_part in ExplicitSolverStrategy constructor", "")

                mpFem_model_part = &(*(settings.fem_model_part));
            if (mpFem_model_part == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "Undefined settings.fem_model_part in ExplicitSolverStrategy constructor", "")

                mpCluster_model_part = &(*(settings.cluster_model_part));
            if (mpCluster_model_part == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "Undefined settings.cluster_model_part in ExplicitSolverStrategy constructor", "")

                mpInlet_model_part = &(*(settings.inlet_model_part));
            if (mpInlet_model_part == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "Undefined settings.inlet_model_part in ExplicitSolverStrategy constructor", "")

            }

        /// Destructor.
        virtual ~ExplicitSolverStrategy() {
            Timer::SetOuputFile("TimesPartialRelease");
            Timer::PrintTimingInformation();
        }

        struct LessX {
            bool operator()(const SphericParticle* p, const SphericParticle* q) const {return p->GetGeometry()[0].Coordinates()[0] < q->GetGeometry()[0].Coordinates()[0];}
        };
        struct LessY {
            bool operator()(const SphericParticle* p, const SphericParticle* q) const {return p->GetGeometry()[0].Coordinates()[1] < q->GetGeometry()[0].Coordinates()[1];}
        };
        struct LessZ {
            bool operator()(const SphericParticle* p, const SphericParticle* q) const {return p->GetGeometry()[0].Coordinates()[2] < q->GetGeometry()[0].Coordinates()[2];}
        };
        struct SpatialSortingTraits {
            typedef SphericParticle* Point_2;
            typedef LessX Less_x_2;
            typedef LessY Less_y_2;
            typedef LessZ Less_z_2;
            Less_x_2 less_x_2_object() const {return Less_x_2();}
            Less_y_2 less_y_2_object() const {return Less_y_2();}
            Less_z_2 less_z_2_object() const { return Less_z_2();}
        };

#ifdef USING_CGAL
        void ReorderParticles() {
            SpatialSortingTraits sst;
            CGAL::spatial_sort(mListOfSphericParticles.begin(), mListOfSphericParticles.end(), sst);
        }
#endif
        template <class T>
        void RebuildListOfSphericParticles(ElementsArrayType& pElements, std::vector<T*>& rCustomListOfParticles){
            KRATOS_TRY
            rCustomListOfParticles.resize(pElements.size());

            #pragma omp parallel for
            for (int k = 0; k < (int)pElements.size(); k++){
              typename ElementsArrayType::iterator particle_pointer_it = pElements.ptr_begin() + k;
              T* spheric_particle = dynamic_cast<T*>(&(*particle_pointer_it));
              rCustomListOfParticles[k] = spheric_particle;
            }
            return;
            KRATOS_CATCH("")
        }
        void RebuildListOfDiscontinuumSphericParticles() {
            RebuildListOfSphericParticles<SphericParticle>(GetModelPart().GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        }
       
        void RebuildPropertiesProxyPointers(std::vector<SphericParticle*>& rCustomListOfSphericParticles);
        void SendProcessInfoToClustersModelPart();
        void UpdateMaxIdOfCreatorDestructor();
        void RepairPointersToNormalProperties(std::vector<SphericParticle*>& rCustomListOfSphericParticles);
        virtual void Initialize();
        virtual void CalculateMaxTimeStep();
        virtual void InitializeClusters();
        virtual void GetClustersForce();
        virtual double Solve();
        void SearchDEMOperations(ModelPart& r_model_part, bool has_mpi = true);
        void SearchFEMOperations(ModelPart& r_model_part, bool has_mpi = true) ;
        virtual void ForceOperations(ModelPart& r_model_part);
        void InitialTimeStepCalculation(); //TODO: remove this one
        void GetForce();
        void FastGetForce();
        virtual void PerformTimeIntegrationOfMotion(int StepFlag = 0);
        void InitializeSolutionStep();
        virtual void BoundingBoxUtility(bool is_time_to_mark_and_remove = true);
        virtual void FinalizeSolutionStep();
        void InitializeElements();
        void InitializeDEMElements();
        void InitializeFEMElements();
        virtual void CalculateConditionsRHSAndAdd();
        void ClearFEMForces();
        void CalculateNodalPressuresAndStressesOnWalls();        
        void SetFlagAndVariableToNodes(const Kratos::Flags& r_flag_name, ComponentOf3ComponentsVariableType& r_variable_to_set, const double value, NodesArrayType& r_nodes_array);
        void SetVariableToNodes(ComponentOf3ComponentsVariableType& r_variable_to_set, const double value, NodesArrayType& r_nodes_array);
        void ResetPrescribedMotionFlagsRespectingImposedDofs();
        void ApplyPrescribedBoundaryConditions();
        void ApplyInitialConditions();
        void SetSearchRadiiOnAllParticles(ModelPart& r_model_part, const double added_search_distance = 0.0, const double amplification = 1.0);
        void SetNormalRadiiOnAllParticles(ModelPart& r_model_part);
        void SetSearchRadiiWithFemOnAllParticles(ModelPart& r_model_part, const double added_search_distance = 0.0, const double amplification = 1.0);
        virtual void SearchNeighbours();
        virtual void ComputeNewNeighboursHistoricalData();
        virtual void ComputeNewRigidFaceNeighboursHistoricalData();
        virtual void SearchRigidFaceNeighbours();
        void DoubleHierarchyMethod();
        /* This should work only with one iteration, but it with mpi does not */
        void CalculateInitialMaxIndentations();
        void PrepareContactModelPart(ModelPart& r_model_part, ModelPart& mcontacts_model_part);
        void PrepareElementsForPrinting();
        void SynchronizeHistoricalVariables(ModelPart& r_model_part);
        void SynchronizeRHS(ModelPart& r_model_part);
        void CleanEnergies();
        void GlobalDamping();

        
        ModelPart& GetModelPart() { return (*mpDem_model_part);}
        ModelPart& GetFemModelPart() { return (*mpFem_model_part);}
        ModelPart& GetContactModelPart() { return (*mpContact_model_part);}
        ModelPart& GetClusterModelPart() { return (*mpCluster_model_part);}
        ModelPart& GetInletModelPart() { return (*mpInlet_model_part);}
        
        VectorResultElementsContainerType& GetResults() { return (mResults);}
        VectorDistanceType& GetResultsDistances() { return (mResultsDistances);}
        RadiusArrayType& GetArrayOfAmplifiedRadii() { return (mArrayOfAmplifiedRadii);}
        int& GetNStepSearch() { return (mNStepSearch);}
        int& GetSearchControl() { return mSearchControl;}
        int& GetNumberOfThreads() { return (mNumberOfThreads);}
        double& GetMaxTimeStep() { return (mMaxTimeStep);}
        double& GetSafetyFactor() { return (mSafetyFactor);}
        int& GetDeltaOption() { return (mDeltaOption);}
        vector<unsigned int>& GetElementPartition() { return (mElementPartition);}
        typename ParticleCreatorDestructor::Pointer& GetParticleCreatorDestructor() { return (mpParticleCreatorDestructor);}
        typename DEMIntegrationScheme::Pointer& GetScheme() { return (mpScheme);}
        typename SpatialSearch::Pointer& GetSpSearch() { return (mpSpSearch);}
        VectorResultConditionsContainerType& GetRigidFaceResults() { return (mRigidFaceResults);}
        VectorDistanceType& GetRigidFaceResultsDistances() { return (mRigidFaceResultsDistances);}
        vector<unsigned int>& GetConditionPartition() { return (mConditionPartition);}
        typename DEM_FEM_Search::Pointer& GetDemFemSearch() { return (mpDemFemSearch);}
        std::vector<PropertiesProxy> mFastProperties;
        virtual ElementsArrayType& GetElements(ModelPart& r_model_part) { return r_model_part.GetCommunicator().LocalMesh().Elements();}

    protected:

        VectorResultElementsContainerType mResults;
        VectorDistanceType mResultsDistances;
        RadiusArrayType mArrayOfAmplifiedRadii;
        int mNStepSearch;
        int mSearchControl;
        int mNumberOfThreads;
        double mMaxTimeStep;
        double mSafetyFactor;
        int mDeltaOption;
        vector<unsigned int> mElementPartition;
        typename ParticleCreatorDestructor::Pointer mpParticleCreatorDestructor;
        typename DEM_FEM_Search::Pointer mpDemFemSearch;
        typename DEMIntegrationScheme::Pointer mpScheme;
        typename SpatialSearch::Pointer mpSpSearch;
        bool mDoSearchNeighbourElements;
        VectorResultConditionsContainerType mRigidFaceResults;
        VectorDistanceType mRigidFaceResultsDistances;
        vector<unsigned int> mConditionPartition;
        ModelPart *mpFem_model_part;
        ModelPart *mpDem_model_part;
        ModelPart *mpInlet_model_part;
        ModelPart *mpContact_model_part;
        ModelPart *mpCluster_model_part;
        std::vector<SphericParticle*> mListOfSphericParticles;
        std::vector<SphericParticle*> mListOfGhostSphericParticles;

    }; // Class ExplicitSolverStrategy

} // namespace Kratos.

#endif // KRATOS_EXPLICIT_SOLVER_STRATEGY  defined
