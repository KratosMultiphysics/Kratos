//
// Authors:
// Miguel Angel Celigueta maceli@cimne.upc.edu
// Miquel Santasusana msantasusana@cimne.upc.edu
//


#if !defined(KRATOS_EXPLICIT_SOLVER_STRATEGY)
#define KRATOS_EXPLICIT_SOLVER_STRATEGY

#ifdef   SINGLE
#define  REAL float
#else    // not SINGLE
#define  REAL double
#endif   // not SINGLE

#ifndef TRILIBRARY
#define TRILIBRARY
#endif

#include "triangle.h"

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

#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/dem_integration_scheme.h"
#include "custom_utilities/create_and_destroy.h"
#include "custom_utilities/dem_fem_utilities.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/inlet.h"

#include "custom_elements/cluster3D.h"
#include "custom_elements/rigid_body_element.h"
////Cfeng
#include "custom_utilities/dem_fem_search.h"
#include "custom_utilities/discrete_particle_configure.h"
#include "custom_utilities/node_configure.h"
#include "custom_utilities/rigid_face_geometrical_object_configure.h"

#ifdef USING_CGAL
#include <CGAL/spatial_sort.h>
#endif

namespace Kratos {

    extern "C" {
      void triangulate(char*, struct triangulateio*, struct triangulateio*, struct triangulateio*);
      void trifree(void*);
    }

    class ExplicitSolverSettings {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverSettings);

        ExplicitSolverSettings() {
        }

        ~ExplicitSolverSettings() {
        }
        ModelPart* r_model_part;
        ModelPart* contact_model_part;
        ModelPart* fem_model_part;
        ModelPart* cluster_model_part;
        ModelPart* inlet_model_part;
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
        typedef PropertiesContainerType::iterator PropertiesIterator;
        typedef DiscreteParticleConfigure<3> ElementConfigureType;
        typedef NodeConfigure<3> NodeConfigureType;
        typedef RigidFaceGeometricalObjectConfigure<3> RigidFaceGeometricalConfigureType;
        typedef Variable<double> ComponentOf3ComponentsVariableType;

        /// Pointer definition of ExplicitSolverStrategy
        KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverStrategy);

        ExplicitSolverStrategy() {
        }

        ExplicitSolverStrategy(ExplicitSolverSettings& settings,
                                const double max_delta_time,
                                const int n_step_search,
                                const double safety_factor,
                                const int delta_option,
                                ParticleCreatorDestructor::Pointer p_creator_destructor,
                                DEM_FEM_Search::Pointer p_dem_fem_search,
                                SpatialSearch::Pointer pSpSearch,
                                Parameters strategy_parameters) {

            mParameters = strategy_parameters;
            mDeltaOption = delta_option;
            mpParticleCreatorDestructor = p_creator_destructor;
            mpDemFemSearch = p_dem_fem_search;
            mpSpSearch = pSpSearch;

            //Also checks old flag name for backward compatibility issues.
            if(mParameters["do_search_dem_neighbours"].GetBool()) {
                mDoSearchNeighbourElements = true;
            } else mDoSearchNeighbourElements = false;
            p_creator_destructor->SetDoSearchNeighbourElements(mDoSearchNeighbourElements);

            if(mParameters["do_search_fem_neighbours"].GetBool()) mDoSearchNeighbourFEMElements = true;
            else mDoSearchNeighbourFEMElements = false;

            mMaxTimeStep = max_delta_time;
            mNStepSearch = n_step_search;
            mSafetyFactor = safety_factor;

            mpDem_model_part = &(*(settings.r_model_part));
            KRATOS_ERROR_IF(mpDem_model_part == NULL) << "Undefined settings.r_model_part in ExplicitSolverStrategy constructor" << std::endl;

            mpContact_model_part = &(*(settings.contact_model_part));
            KRATOS_ERROR_IF(mpContact_model_part == NULL) << "Undefined settings.contact_model_part in ExplicitSolverStrategy constructor" << std::endl;

            mpFem_model_part = &(*(settings.fem_model_part));
            KRATOS_ERROR_IF(mpFem_model_part == NULL) << "Undefined settings.fem_model_part in ExplicitSolverStrategy constructor" << std::endl;

            mpCluster_model_part = &(*(settings.cluster_model_part));
            KRATOS_ERROR_IF(mpCluster_model_part == NULL) << "Undefined settings.cluster_model_part in ExplicitSolverStrategy constructor" << std::endl;

            mpInlet_model_part = &(*(settings.inlet_model_part));
            KRATOS_ERROR_IF(mpInlet_model_part == NULL) << "Undefined settings.inlet_model_part in ExplicitSolverStrategy constructor" << std::endl;

            if(mParameters["RemoveBallsInitiallyTouchingWalls"].GetBool()) mRemoveBallsInitiallyTouchingWallsOption = true;
            else mRemoveBallsInitiallyTouchingWallsOption = false;

        }

        /// Destructor.
        virtual ~ExplicitSolverStrategy() {
            //Timer::SetOuputFile("TimesPartialRelease");
            //Timer::PrintTimingInformation();
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
                ElementsArrayType::iterator particle_pointer_it = pElements.ptr_begin() + k;
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
        virtual void AttachSpheresToStickyWalls();
        virtual void DisplayThreadInfo();
        double CalculateMaxInletTimeStep();
        virtual void InitializeClusters();
        virtual void GetClustersForce();
        virtual void GetRigidBodyElementsForce();
        virtual double SolveSolutionStep();
        void SearchDEMOperations(ModelPart& r_model_part, bool has_mpi = true);
        void SearchFEMOperations(ModelPart& r_model_part, bool has_mpi = true) ;
        virtual void ForceOperations(ModelPart& r_model_part);
        void GetForce();
        void FastGetForce();
        virtual void PerformTimeIntegrationOfMotion(int StepFlag = 0);
        void InitializeSolutionStep();
        virtual void BoundingBoxUtility(bool is_time_to_mark_and_remove = true);
        virtual void FinalizeSolutionStep();
        void InitializeElements();
        void InitializeDEMElements();
        void InitializeFEMElements();
        //void InitializeRigidBodyElements();
        void InitializeFEMWallsAsRigidBodyElements(ModelPart::SubModelPartsContainerType::iterator& sub_model_part);
        void MarkToDeleteAllSpheresInitiallyIndentedWithFEM(ModelPart& rSpheresModelPart);
        void ComputeNodalArea();
        void ComputeNormalPressureVectorField();
        virtual void CalculateConditionsRHSAndAdd();
        void ClearFEMForces();
        void CalculateNodalPressuresAndStressesOnWalls();
        void SetFlagAndVariableToNodes(const Kratos::Flags& r_flag_name, ComponentOf3ComponentsVariableType& r_variable_to_set, const double value, NodesArrayType& r_nodes_array);
        void SetVariableToNodes(ComponentOf3ComponentsVariableType& r_variable_to_set, const double value, NodesArrayType& r_nodes_array);
        void ResetPrescribedMotionFlagsRespectingImposedDofs();
        void ApplyPrescribedBoundaryConditions();
        void ApplyInitialConditions();
        virtual void SetSearchRadiiOnAllParticles(ModelPart& r_model_part, const double added_search_distance = 0.0, const double amplification = 1.0);
        void SetNormalRadiiOnAllParticlesBeforeInitilization(ModelPart& r_model_part);
        void SetNormalRadiiOnAllParticles(ModelPart& r_model_part);
        virtual void SetSearchRadiiWithFemOnAllParticles(ModelPart& r_model_part, const double added_search_distance = 0.0, const double amplification = 1.0);
        virtual void SearchNeighbours();
        virtual void ComputeNewNeighboursHistoricalData();
        virtual void CreateContactElements();
        void InitializeContactElements();
        // void ContactInitializeSolutionStep();
        void PrepareContactElementsForPrinting();
        virtual void ComputeNewRigidFaceNeighboursHistoricalData();
        virtual void SearchRigidFaceNeighbours();
        void CheckHierarchyWithCurrentNeighbours();
        /* This should work only with one iteration, but it with mpi does not */
        void CalculateInitialMaxIndentations(const ProcessInfo& r_process_info);
        void PrepareContactModelPart(ModelPart& r_model_part, ModelPart& mcontacts_model_part);
        void PrepareElementsForPrinting();
        void SynchronizeHistoricalVariables(ModelPart& r_model_part);
        void SynchronizeRHS(ModelPart& r_model_part);
        void Check_MPI(bool& has_mpi);
        virtual double ComputeCoordinationNumber(double& standard_dev);

        ModelPart& GetModelPart() { return (*mpDem_model_part);}
        ModelPart& GetFemModelPart() { return (*mpFem_model_part);}
        ModelPart& GetContactModelPart() { return (*mpContact_model_part);}
        ModelPart& GetClusterModelPart() { return (*mpCluster_model_part);}
        ModelPart& GetInletModelPart() { return (*mpInlet_model_part);}
        ModelPart& GetRigidBodyModelPart() { return (*mpRigidBody_model_part);}

        VectorResultElementsContainerType& GetResults() { return (mResults);}
        VectorDistanceType& GetResultsDistances() { return (mResultsDistances);}
        RadiusArrayType& GetArrayOfAmplifiedRadii() { return (mArrayOfAmplifiedRadii);}
        int& GetNStepSearch() { return (mNStepSearch);}
        int& GetSearchControl() { return mSearchControl;}
        int& GetNumberOfThreads() { return (mNumberOfThreads);}
        double& GetMaxTimeStep() { return (mMaxTimeStep);}
        double& GetSafetyFactor() { return (mSafetyFactor);}
        int& GetDeltaOption() { return (mDeltaOption);}
        ParticleCreatorDestructor::Pointer& GetParticleCreatorDestructor() { return (mpParticleCreatorDestructor);}
        SpatialSearch::Pointer& GetSpSearch() { return (mpSpSearch);}
        VectorResultConditionsContainerType& GetRigidFaceResults() { return (mRigidFaceResults);}
        VectorDistanceType& GetRigidFaceResultsDistances() { return (mRigidFaceResultsDistances);}
        DEM_FEM_Search::Pointer& GetDemFemSearch() { return (mpDemFemSearch);}
        virtual ElementsArrayType& GetElements(ModelPart& r_model_part) { return r_model_part.GetCommunicator().LocalMesh().Elements();}
        virtual ElementsArrayType& GetAllElements(ModelPart& r_model_part) {
            return r_model_part.Elements();
        }

        //==========================================================================================================================================
        // HIERARCHICAL MULTISCALE RVE - START
        //==========================================================================================================================================

        // Properties
        bool   mRVE_FlatWalls;          // Flag for flat rigid walls
        bool   mRVE_Solve;              // Flag for evaluating RVE in current step
        bool   mRVE_Compress;           // Flag for compressing RVE
        bool   mRVE_Equilibrium;        // Flag for static equilibrium of particles
        int    mRVE_Dimension;          // RVE dimension: 2D or 3D
        int    mRVE_EqSteps;            // Number of RVE solution steps in equilibrium
        int    mRVE_NumParticles;       // Total number of particles inside RVE (does not consider wall particles)
        int    mRVE_NumParticlesInner;  // Total number of inner particles (not in contact with walls)
        int    mRVE_NumParticlesWalls;  // Total number of wall particles
        int    mRVE_NumContacts;        // Total number of contacts in RVE (all contacts)
        int    mRVE_NumContactsInner;   // Total number of inner contacts in RVE (considers only contacts involving inner particles)
        double mRVE_MeanRadius;         // Mean radius of all particles
        double mRVE_AvgCoordNum;        // Average coordination number per particle
        double mRVE_AvgCoordNumInner;   // Average coordination number of inner particles (not in contact with walls)
        double mRVE_VolSolid;           // Volume of solid (particles) in RVE discounting overlaps
        double mRVE_VolSolidInner;      // Volume of solid (particles) in convex hull discounting overlaps
        double mRVE_VolTotal;           // RVE total volume (volume inside flat walls)
        double mRVE_VolInner;           // RVE inner volume (considering only inner particles)
        double mRVE_Porosity;           // RVE porosity considering full RVE volume
        double mRVE_PorosityInner;      // RVE porosity considering convex hull volume
        double mRVE_VoidRatio;          // RVE void ratio considering full RVE volume
        double mRVE_VoidRatioInner;     // RVE void ratio considering convex hull volume
        double mRVE_Anisotropy;         // Fabric anisotropy (all particles)
        double mRVE_AnisotropyInner;    // Fabric anisotropy (inner particles)
        double mRVE_EffectStress;       // Mean effective stress (all particles)
        double mRVE_EffectStressInner;  // Mean effective stress (inner particles)
        double mRVE_DevStress;          // Deviatoric stress (all particles)
        double mRVE_DevStressInner;     // Deviatoric stress (inner particles)
        double mRVE_WallForces;         // Total force applied by walls (normal only)
        double mRVE_WallStress;         // Total stress applied by walls (normal only)
        double mRVE_StdDevRoseXYAll;    // Std dev of rose values in XY plane for all particles (divided by number of particles)
        double mRVE_StdDevRoseAzAll;    // Std dev of rose values in Az plane for all particles (divided by number of particles)
        double mRVE_StdDevRoseXYInn;    // Std dev of rose values in XY plane for inner particles (divided by number of particles)
        double mRVE_StdDevRoseAzInn;    // Std dev of rose values in Az plane for inner particles (divided by number of particles)

        std::vector<double> mRVE_CornerCoordsX; // X coordinates of 4 RVE corners (low-left, low-right, up-right, up-left)
        std::vector<double> mRVE_CornerCoordsY; // X coordinates of 4 RVE corners (low-left, low-right, up-right, up-left)

        std::vector<DEMWall*> mRVE_WallXMin; // Vector of RVE flat walls in negative X direction
        std::vector<DEMWall*> mRVE_WallXMax; // Vector of RVE flat walls in positive X direction
        std::vector<DEMWall*> mRVE_WallYMin; // Vector of RVE flat walls in negative Y direction
        std::vector<DEMWall*> mRVE_WallYMax; // Vector of RVE flat walls in positive Y direction
        std::vector<DEMWall*> mRVE_WallZMin; // Vector of RVE flat walls in negative Z direction
        std::vector<DEMWall*> mRVE_WallZMax; // Vector of RVE flat walls in positive Z direction

        std::vector<SphericParticle*> mRVE_InnerVolParticles; // Vector of particles composing the inner volume

        std::vector<SphericParticle*> mRVE_WallParticleXMin; // Vector of RVE particle walls in negative X direction
        std::vector<SphericParticle*> mRVE_WallParticleXMax; // Vector of RVE particle walls in positive X direction
        std::vector<SphericParticle*> mRVE_WallParticleYMin; // Vector of RVE particle walls in negative Y direction
        std::vector<SphericParticle*> mRVE_WallParticleYMax; // Vector of RVE particle walls in positive Y direction
        std::vector<SphericParticle*> mRVE_WallParticleZMin; // Vector of RVE particle walls in negative Z direction
        std::vector<SphericParticle*> mRVE_WallParticleZMax; // Vector of RVE particle walls in positive Z direction

        std::vector<double> mRVE_ForceChain; // Vector of force chains coordinates: [x1,y1,z1,x2,y2,z2,F, x1,y1,z1,x2,y2,z2,F, x1,y1,z1,x2,y2,z2,F, ...]
        Matrix mRVE_RoseDiagram;             // Rose diagram of contacts: Row 1 = angle ranges in plane XY; Row 2 = azimute ranges wrt to plane XY;
        Matrix mRVE_RoseDiagramInner;        // Rose diagram of contacts for inner particles (not in contact with walls)
        Matrix mRVE_FabricTensor;            // Fabric tensor (all particles)
        Matrix mRVE_FabricTensorInner;       // Fabric tensor (inner particles)
        Matrix mRVE_CauchyTensor;            // Cauchy stress tensor (all particles)
        Matrix mRVE_CauchyTensorInner;       // Cauchy stress tensor (inner particles)
        Matrix mRVE_TangentTensor;           // Tangent operator tensor (all particles)
        Matrix mRVE_TangentTensorInner;      // Tangent operator tensor (inner particles)
        Matrix mRVE_ConductivityTensor;      // Effective thermal conductivity tensor (all particles)
        Matrix mRVE_ConductivityTensorInner; // Effective thermal conductivity tensor (inner particles)

        std::ofstream mRVE_FileCoordinatesHistory;
        std::ofstream mRVE_FileCoordinatesLast;
        std::ofstream mRVE_FilePorosity;
        std::ofstream mRVE_FileContactNumber;
        std::ofstream mRVE_FileCoordNumber;
        std::ofstream mRVE_FileInnerVolumeParticles;
        std::ofstream mRVE_FileForceChain;
        std::ofstream mRVE_FileElasticContactForces;
        std::ofstream mRVE_FileElasticContactForcesWalls;
        std::ofstream mRVE_FileRoseDiagram;
        std::ofstream mRVE_FileRoseDiagramInner;
        std::ofstream mRVE_FileRoseDiagramUniformity;
        std::ofstream mRVE_FileAnisotropy;
        std::ofstream mRVE_FileFabricTensor;
        std::ofstream mRVE_FileFabricTensorInner;
        std::ofstream mRVE_FileStress;
        std::ofstream mRVE_FileCauchyTensor;
        std::ofstream mRVE_FileCauchyTensorInner;
        std::ofstream mRVE_FileTangentTensor;
        std::ofstream mRVE_FileTangentTensorInner;
        std::ofstream mRVE_FileConductivityTensor;
        std::ofstream mRVE_FileConductivityTensorInner;
        std::ofstream mRVE_FileFKS;

        // Methods
        void RVEInitialize             (void);
        void RVEInitializeSolutionStep (void);
        void RVEExecuteParticlePre     (SphericParticle* p_particle);
        void RVEExecuteParticlePos     (SphericParticle* p_particle);
        void RVEFinalizeSolutionStep   (void);
        void Finalize                  (void);
        void RVEFinalize               (void);

        void RVEAssembleWallVectors             (void);
        void RVEAssembleWallVectors2D_Flat      (void);
        void RVEAssembleWallVectors3D_Flat      (void);
        void RVEAssembleWallVectors2D_Particles (void);
        void RVEAssembleWallVectors3D_Particles (void);

        void   RVEComputeCorners        (void);
        double RVEComputeTotalSurface   (void);
        double RVEComputeTotalVolume    (void);
        double RVEComputeInnerVolume    (void);
        double RVEComputeParticleVolume (SphericParticle* p_particle);
        void   RVEComputePorosity       (void);
        void   RVEComputeInnerPorosity  (void);
        void   RVEHomogenization        (void);
        void   RVEComputeRoseUniformity (void);
        void   RVEStopCompression       (void);

        void RVEReadOldForces       (void);
        void RVEWriteFiles          (void);
        void RVEWriteCoords         (void);
        void RVEWriteForceParticles (void);
        void RVEWriteForceContacts  (void);
        void RVEOpenFiles           (void);
        void RVECloseFiles          (void);

        void ClearTriangle (struct triangulateio& rTr);
        void FreeTriangle  (struct triangulateio& rTr);

        //==========================================================================================================================================
        // HIERARCHICAL MULTISCALE RVE - FINISH
        //==========================================================================================================================================

    protected:

        Parameters mParameters;
        bool mRemoveBallsInitiallyTouchingWallsOption;
        VectorResultElementsContainerType mResults;
        VectorDistanceType mResultsDistances;
        RadiusArrayType mArrayOfAmplifiedRadii;
        int mNStepSearch;
        int mSearchControl;
        int mNumberOfThreads;
        double mMaxTimeStep;
        double mSafetyFactor;
        int mDeltaOption;
        ParticleCreatorDestructor::Pointer mpParticleCreatorDestructor;
        DEM_FEM_Search::Pointer mpDemFemSearch;
        SpatialSearch::Pointer mpSpSearch;
        bool mDoSearchNeighbourElements;
        bool mDoSearchNeighbourFEMElements;
        VectorResultConditionsContainerType mRigidFaceResults;
        VectorDistanceType mRigidFaceResultsDistances;
        ModelPart *mpFem_model_part;
        ModelPart *mpDem_model_part;
        ModelPart *mpInlet_model_part;
        ModelPart *mpContact_model_part;
        ModelPart *mpCluster_model_part;
        ModelPart *mpRigidBody_model_part;
        std::vector<SphericParticle*> mListOfSphericParticles;
        std::vector<SphericParticle*> mListOfGhostSphericParticles;

    }; // Class ExplicitSolverStrategy

} // namespace Kratos.

#endif // KRATOS_EXPLICIT_SOLVER_STRATEGY  defined
