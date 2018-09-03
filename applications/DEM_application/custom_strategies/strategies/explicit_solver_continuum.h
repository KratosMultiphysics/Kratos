//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_CONTINUUM_EXPLICIT_SOLVER_STRATEGY)
#define  KRATOS_CONTINUUM_EXPLICIT_SOLVER_STRATEGY
#include "custom_strategies/strategies/explicit_solver_strategy.h"
#include "custom_elements/spheric_continuum_particle.h"
#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) ContinuumExplicitSolverStrategy : public ExplicitSolverStrategy {
    public:

        typedef ExplicitSolverStrategy BaseType;
        typedef BaseType::NodesArrayType NodesArrayType;
        typedef BaseType::ElementsArrayType ElementsArrayType;
        typedef BaseType::ElementsIterator ElementsIterator;
        typedef BaseType::ConditionsArrayType ConditionsArrayType;
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

        /// Pointer definition of ExplicitSolverStrategy
        KRATOS_CLASS_POINTER_DEFINITION(ContinuumExplicitSolverStrategy);

        /// Default constructor.

        ContinuumExplicitSolverStrategy() {
        }

        ContinuumExplicitSolverStrategy(
                ExplicitSolverSettings& settings,
                const double max_delta_time,
                const int n_step_search,
                const double safety_factor,
                const int delta_option,
                ParticleCreatorDestructor::Pointer p_creator_destructor,
                DEM_FEM_Search::Pointer p_dem_fem_search,
                SpatialSearch::Pointer pSpSearch)
        : ExplicitSolverStrategy(settings, max_delta_time, n_step_search, safety_factor, delta_option, p_creator_destructor, p_dem_fem_search, pSpSearch) {
            BaseType::GetParticleCreatorDestructor() = p_creator_destructor;
        }

        /// Destructor.

        virtual ~ContinuumExplicitSolverStrategy() {
            //Timer::SetOuputFile("TimesPartialRelease");
            //Timer::PrintTimingInformation();
        }

        virtual void Initialize() override;
        virtual double Solve() override;
        void SearchFEMOperations(ModelPart& r_model_part, bool has_mpi);
        void SearchDEMOperations(ModelPart& r_model_part, bool has_mpi);
        void ComputeNewNeighboursHistoricalData() override;
        void ComputeNewRigidFaceNeighboursHistoricalData() override;
        void CreateContactElements();
        void InitializeContactElements();
        void ContactInitializeSolutionStep();
        void PrepareContactElementsForPrinting();
        void SetCoordinationNumber(ModelPart& r_model_part);
        double ComputeCoordinationNumber(double& standard_dev);
        void BoundingBoxUtility(bool is_time_to_mark_and_remove = true) override;
        void Check_MPI(bool& has_mpi);
        virtual void CalculateMaxSearchDistance();
        virtual void MeshRepairOperations();
        virtual void DestroyMarkedParticlesRebuildLists();
        void CalculateMeanContactArea();
        void SetInitialDemContacts();
        void SetInitialFemContacts();
        void FinalizeSolutionStep() override;
        void FinalizeSolutionStepFEM();
        void MarkNewSkinParticles();

        virtual void Add_As_Own(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element) {
            KRATOS_TRY
            mcontacts_model_part.Elements().push_back(p_contact_element);
            KRATOS_CATCH("")
        }
        //En aquest cas m'afegeixo jo al local i a la interface local corresponent amb la particio del vei ghost

        virtual void Add_As_Local(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element) {
        }
        //I aqui m'afegeixio yo com a ghost de la particio del vei local

        virtual void Add_As_Ghost(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element) {
        }

        virtual void Sort_Contact_Modelpart(ModelPart& mcontacts_model_part) {
        }

        virtual void Reassign_Ids(ModelPart& mcontacts_model_part) {
        }

        virtual ElementsArrayType& GetAllElements(ModelPart& r_model_part) {
            return r_model_part.Elements();
        }

        virtual ElementsArrayType& GetElements(ModelPart& r_model_part) override {
            return r_model_part.GetCommunicator().LocalMesh().Elements();
        }

    protected:

        bool mcontinuum_simulating_option;
        int mFixSwitch;
        //bool   mDempackOption;
        std::vector<SphericContinuumParticle*> mListOfSphericContinuumParticles;
        std::vector<SphericContinuumParticle*> mListOfGhostSphericContinuumParticles;

    }; // Class ContinuumExplicitSolverStrategy

} // namespace Kratos

#endif // KRATOS_FILENAME_H_INCLUDED  defined
