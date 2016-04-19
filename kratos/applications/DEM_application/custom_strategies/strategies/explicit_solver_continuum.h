//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_CONTINUUM_EXPLICIT_SOLVER_STRATEGY)
#define  KRATOS_CONTINUUM_EXPLICIT_SOLVER_STRATEGY

#include "custom_strategies/strategies/explicit_solver_strategy.h"
#include "DEM_definitions.h"

#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

namespace Kratos
{
    template<
    class TSparseSpace,
    class TDenseSpace,
    class TLinearSolver>
    class ContinuumExplicitSolverStrategy: public ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
    {
        public:
        ///@name Type Definitions
        ///@{

        typedef ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>   BaseType;

        typedef typename BaseType::NodesArrayType                             NodesArrayType;
        typedef typename BaseType::ElementsArrayType                          ElementsArrayType;
        typedef typename BaseType::ElementsIterator                           ElementsIterator;
        typedef typename BaseType::ConditionsArrayType                        ConditionsArrayType;

        /*  Revisar charlie */

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
        using BaseType::mNeighbourCounter;

        /// Pointer definition of ExplicitSolverStrategy

        KRATOS_CLASS_POINTER_DEFINITION(ContinuumExplicitSolverStrategy);

        /// Default constructor.
        ContinuumExplicitSolverStrategy(){}

        ContinuumExplicitSolverStrategy(
                               ExplicitSolverSettings& settings,
                               const double max_delta_time,
                               const int n_step_search,
                               const double safety_factor,
                               const int delta_option,
                               typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                               typename DEM_FEM_Search::Pointer p_dem_fem_search,
                               typename DEMIntegrationScheme::Pointer pScheme,
                               typename SpatialSearch::Pointer pSpSearch)
        :ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(settings, max_delta_time, n_step_search, safety_factor, delta_option, p_creator_destructor, p_dem_fem_search, pScheme, pSpSearch)
        {                    
            BaseType::GetParticleCreatorDestructor()   = p_creator_destructor;                            
        }

        /// Destructor.
        virtual ~ContinuumExplicitSolverStrategy()
        {
           Timer::SetOuputFile("TimesPartialRelease");
           Timer::PrintTimingInformation();
        }      
        
        virtual void Initialize()
        {
            KRATOS_TRY

            std::cout << "------------------CONTINUUM SOLVER STRATEGY---------------------" << "\n" <<std::endl;

            ModelPart& r_model_part     = GetModelPart();
            ModelPart& fem_model_part   = GetFemModelPart();
            ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

            // Omp initializations
            mNumberOfThreads = OpenMPUtils::GetNumThreads();

            std::cout << "          **************************************************" << std::endl;
            std::cout << "            Parallelism Info:  MPI number of nodes: " << r_model_part.GetCommunicator().TotalProcesses() <<std::endl;
            if( r_model_part.GetCommunicator().TotalProcesses() > 1 )
                std::cout << "            Parallelism Info:  MPI node Id: " << r_model_part.GetCommunicator().MyPID() <<std::endl;
            std::cout << "            Parallelism Info:  OMP number of processors: " << mNumberOfThreads <<std::endl;
            std::cout << "          **************************************************" << std::endl << std::endl;

            this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles);
            this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
            this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
            this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

            r_process_info[SEARCH_CONTROL_VECTOR].resize(mNumberOfThreads);
            for (int i=0; i<mNumberOfThreads; i++) r_process_info[SEARCH_CONTROL_VECTOR][i] = 0;
            BaseType::GetNeighbourCounter().resize(mNumberOfThreads);

            CreatePropertiesProxies(BaseType::mFastProperties, r_model_part, *mpInlet_model_part, *mpCluster_model_part);

            BaseType::RepairPointersToNormalProperties(mListOfSphericParticles); 
            BaseType::RepairPointersToNormalProperties(mListOfGhostSphericParticles); 

            BaseType::RebuildPropertiesProxyPointers(mListOfSphericParticles);
            BaseType::RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);        

            BaseType::GetSearchControl()         = r_process_info[SEARCH_CONTROL];

            BaseType::InitializeDEMElements();
            BaseType::InitializeFEMElements();                
            BaseType::InitializeSolutionStep();
            BaseType::ApplyInitialConditions();

            this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles);
            this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
            this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
            this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

            // 0. Set search radius.
            unsigned int number_of_elements = r_model_part.GetCommunicator().LocalMesh().Elements().size();
            if(BaseType::GetResults().size() != number_of_elements) BaseType::GetResults().resize(number_of_elements);        
            BaseType::GetResultsDistances().resize(number_of_elements);                  
            if(fem_model_part.Nodes().size()>0) {
                BaseType::GetRigidFaceResults().resize(number_of_elements);
                BaseType::GetRigidFaceResultsDistances().resize(number_of_elements);
            }                       

            // 3. Search Neighbors with tolerance (after first repartition process)
            SetSearchRadiiOnAllParticles(r_model_part, r_process_info[SEARCH_TOLERANCE], 0.0);
            SearchNeighbours();

            if (BaseType::GetDeltaOption() == 2) {
                SetCoordinationNumber(r_model_part);          
            }

            this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles);
            this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
            this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
            this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

            bool has_mpi = false;
            Check_MPI(has_mpi);

            if (has_mpi) {
                BaseType::RepairPointersToNormalProperties(mListOfSphericParticles);
                BaseType::RepairPointersToNormalProperties(mListOfGhostSphericParticles);
            }

            BaseType::RebuildPropertiesProxyPointers(mListOfSphericParticles);
            BaseType::RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

            if (has_mpi) {
                //RebuildListsOfPointersOfEachParticle(); //Serialized pointers are lost, so we rebuild them using Id's        
            }

            // 4. Set Initial Contacts        
            if (r_process_info[CASE_OPTION] != 0) {            
                SetInitialDemContacts();
                CalculateMeanContactArea();
            }   
            
            ComputeNewNeighboursHistoricalData();                    

            if (fem_model_part.Nodes().size()>0) {  
                SetSearchRadiiOnAllParticles(r_model_part, 0.0, 0.0); //Strict radius, not amplified (0.0 added distance, 1.0 factor of amplification)
                BaseType::SearchRigidFaceNeighbours(r_process_info[LOCAL_RESOLUTION_METHOD]);
                SetInitialFemContacts();
                BaseType::ComputeNewRigidFaceNeighboursHistoricalData();        
            }

            if (r_process_info[CONTACT_MESH_OPTION] == 1) {   
                CreateContactElements();
                InitializeContactElements();
            }

            // 5. Finalize Solution Step        
            //BaseType::FinalizeSolutionStep();

            //KRATOS_WATCH(r_model_part.GetNodalSolutionStepVariablesList())        
            KRATOS_CATCH("")        
        }// Initialize()

        virtual double Solve()
        {

            KRATOS_TRY

            ModelPart& r_model_part            = GetModelPart();


            bool has_mpi = false;
            VariablesList r_modelpart_nodal_variables_list = r_model_part.GetNodalSolutionStepVariablesList();
            if(r_modelpart_nodal_variables_list.Has(PARTITION_INDEX) )  has_mpi = true;

            BaseType::InitializeSolutionStep();
            SearchOperations(r_model_part, has_mpi);
            BaseType::ForceOperations(r_model_part);
            BaseType::PerformTimeIntegrationOfMotion();           
            FinalizeSolutionStep();

            KRATOS_CATCH("") 

            return 0.0;

        }//Solve()        
            

        void SearchOperations(ModelPart& r_model_part, bool has_mpi)
        {
            ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

            if(r_process_info[SEARCH_CONTROL]==0) {  
                for (int i = 0; i < mNumberOfThreads; i++) {
                    if(r_process_info[SEARCH_CONTROL_VECTOR][i]==1) {                
                        r_process_info[SEARCH_CONTROL]=1;                          
                        std::cout << "From now on, the search is activated because some failure occurred " <<std::endl;   
                        break;                
                    }
                }             
            }

            const int time_step = r_process_info[TIME_STEPS];
            const double time = r_process_info[TIME];
            const bool is_time_to_search_neighbours = (time_step + 1) % BaseType::GetNStepSearch() == 0 && (time_step > 0); //Neighboring search. Every N times.

            if (r_process_info[SEARCH_CONTROL] > 0) {

                if (is_time_to_search_neighbours) {

                    if (r_process_info[BOUNDING_BOX_OPTION] && time >= r_process_info[BOUNDING_BOX_START_TIME] && time <= r_process_info[BOUNDING_BOX_STOP_TIME]) {

                        BoundingBoxUtility();                        
                    }
                    else {
                        ParticleCreatorDestructor::Pointer& p_creator_destructor=BaseType::GetParticleCreatorDestructor();
                        p_creator_destructor->DestroyParticles(r_model_part);
                        p_creator_destructor->DestroyContactElements(*mpContact_model_part);
                    }

                    this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles); //These lists are necessary for the loop in SearchNeighbours
                    this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles); 

                    
                    SetSearchRadiiOnAllParticles(r_model_part, r_process_info[SEARCH_TOLERANCE], r_process_info[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION]);
                    SearchNeighbours(); //the amplification factor has been modified after the first search.

                    this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles); //These lists are necessary because the elements in this partition might have changed.
                    this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles); 
                    this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
                    this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

                    if(has_mpi){
                      BaseType::RepairPointersToNormalProperties(mListOfSphericParticles);
                      BaseType::RepairPointersToNormalProperties(mListOfGhostSphericParticles);
                      //RebuildListsOfPointersOfEachParticle(); //Serialized member variables which are pointers are lost, so we rebuild them using Id's
                    }

                    BaseType::RebuildPropertiesProxyPointers(mListOfSphericParticles);
                    BaseType::RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

                    ComputeNewNeighboursHistoricalData();

                    BaseType::SetSearchRadiiOnAllParticles(r_model_part, 0.0, 0.0); //Strict radius, not amplified (0.0 added distance, 1.0 factor of amplification)
                    BaseType::SearchRigidFaceNeighbours(r_process_info[LOCAL_RESOLUTION_METHOD]);
                    BaseType::ComputeNewRigidFaceNeighboursHistoricalData();
                    r_process_info[SEARCH_CONTROL] = 2;
                    //if (r_process_info[BOUNDING_BOX_OPTION] == 1 && has_mpi) {  //This block rebuilds all the bonds between continuum particles
                        if (r_process_info[CONTACT_MESH_OPTION] == 1) {                            
                            CreateContactElements();
                            InitializeContactElements();
                        }
                    //}
                }
                else{
                    r_process_info[SEARCH_CONTROL] = 1;
                }
            }
            //Synch this var.
            r_model_part.GetCommunicator().MaxAll(r_process_info[SEARCH_CONTROL]);
        }
      
    void ComputeNewNeighboursHistoricalData() {
        KRATOS_TRY  

        #pragma omp parallel
        {
            std::vector<int>                  TempNeighboursIds; //We are passing all these temporal vectors as arguments because creating them inside the function is slower (memory allocation and deallocation)
            std::vector<array_1d<double, 3> > TempNeighbourElasticContactForces;
            std::vector<array_1d<double, 3> > TempNeighbourTotalContactForces;
            std::vector<SphericParticle*>     TempNeighbourElements;
            
            const int number_of_particles = (int)mListOfSphericContinuumParticles.size();
         
            #pragma omp for
            for(int i=0; i<number_of_particles; i++){
                mListOfSphericContinuumParticles[i]->ReorderAndRecoverInitialPositionsAndFilter(TempNeighbourElements);
                mListOfSphericContinuumParticles[i]->ComputeNewNeighboursHistoricalData(TempNeighboursIds,
                                                                                        TempNeighbourElasticContactForces,
                                                                                        TempNeighbourTotalContactForces);
            }
        }

        KRATOS_CATCH("")
    }
    
    void ComputeNewRigidFaceNeighboursHistoricalData()
    {
        KRATOS_TRY
        const int number_of_particles = (int)mListOfSphericParticles.size();

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++){
            mListOfSphericContinuumParticles[i]->ReorderFEMneighbours();
            mListOfSphericParticles[i]->ComputeNewRigidFaceNeighboursHistoricalData();
        }

        KRATOS_CATCH("")
    }    

    void CreateContactElements() 
    {                
        KRATOS_TRY        
                
        std::string ElementName;
        ElementName = std::string("ParticleContactElement");
        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);
        
        //Here we are going to create contact elements when we are on a target particle and we see a neighbor whose id is higher than ours.
        //We create also a pointer from the node to the element, after creating it.
        //When our particle has a higher ID than the neighbor we also create a pointer to the (previously) created contact element.
        //We proceed in this way because we want to have the pointers to contact elements in a list in the same order than the initial elements order.
                             
        const int number_of_particles = (int)mListOfSphericContinuumParticles.size();        
        int used_bonds_counter = 0; 
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                unsigned int continuous_initial_neighbors_size = mListOfSphericContinuumParticles[i]->mContinuumInitialNeighborsSize;
                mListOfSphericContinuumParticles[i]->mBondElements.resize(continuous_initial_neighbors_size);
                for (unsigned int j=0; j<mListOfSphericContinuumParticles[i]->mBondElements.size(); j++) {
                    mListOfSphericContinuumParticles[i]->mBondElements[j] = NULL;
                }            
            }                
            
            int private_counter = 0;            
            Element::Pointer p_new_contact_element;
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                bool add_new_bond = true;
                std::vector<SphericParticle*>& neighbour_elements = mListOfSphericContinuumParticles[i]->mNeighbourElements;
                unsigned int continuous_initial_neighbors_size = mListOfSphericContinuumParticles[i]->mContinuumInitialNeighborsSize;

                for (unsigned int j = 0; j < continuous_initial_neighbors_size; j++) {
                    SphericContinuumParticle* neighbour_element = dynamic_cast<SphericContinuumParticle*>(neighbour_elements[j]);
                    if (neighbour_element == NULL) continue; //The initial neighbor was deleted at some point in time!!
                    if (mListOfSphericContinuumParticles[i]->Id() > neighbour_element->Id()) continue;                                               
                    
                    #pragma omp critical
                    {
                        if (used_bonds_counter < (int)(*mpContact_model_part).Elements().size()) {
                            add_new_bond = false;
                            private_counter = used_bonds_counter;
                            used_bonds_counter++;
                        }      
                    }
                    if(!add_new_bond) {
                        Element::Pointer& p_old_contact_element = (*mpContact_model_part).Elements().GetContainer()[private_counter];
                        p_old_contact_element->GetGeometry()(0) = mListOfSphericContinuumParticles[i]->GetGeometry()(0);
                        p_old_contact_element->GetGeometry()(1) = neighbour_element->GetGeometry()(0); 
                        p_old_contact_element->SetId(used_bonds_counter);
                        p_old_contact_element->SetProperties(mListOfSphericContinuumParticles[i]->pGetProperties());
                        ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*>( p_old_contact_element.get() );     
                        mListOfSphericContinuumParticles[i]->mBondElements[j] = p_bond;  
                    }
                    else {
                        Geometry<Node<3> >::PointsArrayType  NodeArray(2);
                        NodeArray.GetContainer()[0] = mListOfSphericContinuumParticles[i]->GetGeometry()(0);
                        NodeArray.GetContainer()[1] = neighbour_element->GetGeometry()(0);                                                
                        const Properties::Pointer& properties = mListOfSphericContinuumParticles[i]->pGetProperties();
                        p_new_contact_element = rReferenceElement.Create(used_bonds_counter+1, NodeArray, properties);
                        
                        #pragma omp critical
                        {                            
                            (*mpContact_model_part).Elements().push_back(p_new_contact_element);
                            used_bonds_counter++;
                        } 
                        ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*>( p_new_contact_element.get() );     
                        mListOfSphericContinuumParticles[i]->mBondElements[j] = p_bond;  
                    }
                                                                                    
                }            
            }

            #pragma omp single
            {
                if((int)(*mpContact_model_part).Elements().size() > used_bonds_counter){
                    (*mpContact_model_part).Elements().erase((*mpContact_model_part).Elements().ptr_begin() + used_bonds_counter, (*mpContact_model_part).Elements().ptr_end());  
                }
            } 

            #pragma omp for 
            for (int i = 0; i < number_of_particles; i++) {
                std::vector<SphericParticle*>& neighbour_elements = mListOfSphericContinuumParticles[i]->mNeighbourElements;
                unsigned int continuous_initial_neighbors_size = mListOfSphericContinuumParticles[i]->mContinuumInitialNeighborsSize;

                for (unsigned int j = 0; j < continuous_initial_neighbors_size; j++) {
                    SphericContinuumParticle* neighbour_element = dynamic_cast<SphericContinuumParticle*>(neighbour_elements[j]);
                    if (neighbour_element == NULL) continue; //The initial neighbor was deleted at some point in time!!
                        //ATTENTION: Ghost nodes do not have mContinuumIniNeighbourElements in general, so this bond will remain as NULL!! 
                    if (mListOfSphericContinuumParticles[i]->Id() < neighbour_element->Id()) continue;
                        //In all functions using mBondElements we must check that this bond is not used.

                    for (unsigned int k = 0; k < neighbour_element->mContinuumInitialNeighborsSize; k++ ) {
                        //ATTENTION: Ghost nodes do not have mContinuumIniNeighbourElements in general, so this bond will remain as NULL!! 
                        //In all functions using mBondElements we must check that this bond is not used.
                        if (neighbour_element->mNeighbourElements[k] == NULL) continue; //The initial neighbor was deleted at some point in time!!
                        if (neighbour_element->mNeighbourElements[k]->Id() == mListOfSphericContinuumParticles[i]->Id()) {
                            ParticleContactElement* bond = neighbour_element->mBondElements[k];
                            mListOfSphericContinuumParticles[i]->mBondElements[j] = bond; 
                            break;
                        }
                    }
                }
            }           
        } //#pragma omp parallel        
        KRATOS_CATCH("")                     
    } //CreateContactElements    
    
    void InitializeContactElements()
    {
        KRATOS_TRY
        //CONTACT MODEL PART      
        ElementsArrayType& pContactElements = GetAllElements(*mpContact_model_part);      
        vector<unsigned int> contact_element_partition;
        OpenMPUtils::CreatePartition(mNumberOfThreads, pContactElements.size(), contact_element_partition);

        #pragma omp parallel for        
        for (int k = 0; k < mNumberOfThreads; k++) {
            typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
            typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];
            
            for (typename ElementsArrayType::iterator it_contact= it_contact_begin; it_contact!=it_contact_end; ++it_contact) {            
                (it_contact)->Initialize();                 
            } //loop over CONTACT ELEMENTS
        }// loop threads OpenMP

      KRATOS_CATCH("")       
    }
    
    void ContactInitializeSolutionStep()
    {
        ElementsArrayType& pContactElements = GetAllElements(*mpContact_model_part);      
        ProcessInfo& r_process_info  = (*mpContact_model_part).GetProcessInfo();

        vector<unsigned int> contact_element_partition;

        OpenMPUtils::CreatePartition(mNumberOfThreads, pContactElements.size(), contact_element_partition);
        #pragma omp parallel for

        for (int k = 0; k < mNumberOfThreads; k++)
        {
            typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
            typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];

            for (typename ElementsArrayType::iterator it_contact= it_contact_begin; it_contact!=it_contact_end; ++it_contact)               
            {                             
                (it_contact)->InitializeSolutionStep(r_process_info); 
            } //loop over CONTACT ELEMENTS

        }// loop threads OpenMP
          
    } //Contact_InitializeSolutionStep
        
    void PrepareContactElementsForPrinting() {
       
        ElementsArrayType& pContactElements = GetAllElements(*mpContact_model_part);
               
        vector<unsigned int> contact_element_partition;

        OpenMPUtils::CreatePartition(mNumberOfThreads, pContactElements.size(), contact_element_partition);

        #pragma omp parallel for

        for (int k = 0; k < mNumberOfThreads; k++)
        { 
            typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
            typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];

            for (typename ElementsArrayType::iterator it_contact= it_contact_begin; it_contact!=it_contact_end; ++it_contact) {
                Element* raw_p_contact_element = &(*it_contact);
                ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*>( raw_p_contact_element );    
                p_bond->PrepareForPrinting();
            } //loop over CONTACT ELEMENTS
        }// loop threads OpenMP        
        //Important TODO: renumber all id's to avoid repetition across partitions
    } //PrepareContactElementsForPrinting      
      
    void SetCoordinationNumber(ModelPart& r_model_part)
    {
        ProcessInfo& r_process_info      = r_model_part.GetProcessInfo();

        const double in_coordination_number = r_process_info[COORDINATION_NUMBER];
        double out_coordination_number = ComputeCoordinationNumber();
        int iteration = 0;
        int maxiteration = 100;
        double& added_search_distance = r_process_info[SEARCH_TOLERANCE];

        std::cout<<"Setting up Coordination Number by increasing or decreasing the search radius... "<<std::endl;

        if(in_coordination_number <= 0.0) {
            KRATOS_THROW_ERROR(std::runtime_error, "The specified Coordination Number is less or equal to zero, N.C. = ",in_coordination_number)
        }

        while (fabs(out_coordination_number/in_coordination_number - 1.0) > 1e-3) {
            if (iteration>=maxiteration) break;
            iteration++;
            if(out_coordination_number == 0.0) {
                std::cout << "Coordination Number method not supported in this case" << "\n" <<std::endl;
                KRATOS_THROW_ERROR(std::runtime_error, "The specified tangency method is not supported for this problem, please use absolute value instead"," ")
                break;
            }
            added_search_distance *= in_coordination_number/out_coordination_number;
            SetSearchRadiiOnAllParticles(r_model_part, added_search_distance, 0.0);
            SearchNeighbours(); //r_process_info[SEARCH_TOLERANCE] will be used inside this function, and it's the variable we are updating in this while
            out_coordination_number = ComputeCoordinationNumber();
        }//while

        if (iteration<maxiteration) std::cout<< "Coordination Number iteration converged after "<<iteration<< " iterations, to value " <<out_coordination_number<< " using an extension of " << added_search_distance <<". "<<"\n"<<std::endl;
        else {
            std::cout << "Coordination Number iteration did NOT converge after "<<iteration<<" iterations. Coordination number reached is "<<out_coordination_number<<". "<<"\n"<<std::endl;
            KRATOS_THROW_ERROR(std::runtime_error,"Please use a Absolute tolerance instead "," ")
            //NOTE: if it doesn't converge, problems occur with contact mesh and rigid face contact.
        }

    } //SetCoordinationNumber

    double ComputeCoordinationNumber()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        unsigned int total_contacts = 0;
        const int number_of_particles = (int)mListOfSphericParticles.size();
        
        mNumberOfThreads = OpenMPUtils::GetNumThreads();
        mNeighbourCounter.resize(mNumberOfThreads);

        #pragma omp parallel
        {
            mNeighbourCounter[OpenMPUtils::ThisThread()] = 0;
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                  mNeighbourCounter[OpenMPUtils::ThisThread()] += mListOfSphericParticles[i]->mNeighbourElements.size();
            }
        }
        for (int i = 0; i < mNumberOfThreads; i++)
        {
          total_contacts += mNeighbourCounter[i];
        }

        int global_total_contacts = total_contacts;
        r_model_part.GetCommunicator().SumAll(global_total_contacts);
        int global_number_of_elements = (int)pElements.size();
        r_model_part.GetCommunicator().SumAll(global_number_of_elements);

        double coord_number = double(global_total_contacts)/double(global_number_of_elements);

        return coord_number;

        KRATOS_CATCH("")
    }
    
    void BoundingBoxUtility(bool is_time_to_mark_and_remove = true) {
        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info    = r_model_part.GetProcessInfo();
        ParticleCreatorDestructor::Pointer& p_creator_destructor=BaseType::GetParticleCreatorDestructor();

        p_creator_destructor->MarkDistantParticlesForErasing(r_model_part);

        if(r_process_info[CONTACT_MESH_OPTION] == 1) {
            p_creator_destructor->MarkContactElementsForErasing(r_model_part, *mpContact_model_part);
            p_creator_destructor->DestroyContactElements(*mpContact_model_part);
        }
        p_creator_destructor->DestroyParticles(r_model_part);        

        KRATOS_CATCH("")
    }
    
    void Check_MPI(bool& has_mpi){
        VariablesList r_modelpart_nodal_variables_list = GetModelPart().GetNodalSolutionStepVariablesList();
        if(r_modelpart_nodal_variables_list.Has(PARTITION_INDEX) )  has_mpi = true;
    }
       

    void CalculateMeanContactArea() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        bool has_mpi = false;
        Check_MPI(has_mpi);

        const int number_of_particles = (int) mListOfSphericContinuumParticles.size();

        #pragma omp parallel for 
        for (int i = 0; i < number_of_particles; i++) { //Do not do this for the ghost particles!
            mListOfSphericContinuumParticles[i]->CalculateMeanContactArea(has_mpi, r_process_info);
        }

        KRATOS_CATCH("")
    } 
 
    void SetInitialDemContacts()
    {           
        KRATOS_TRY

        ProcessInfo& r_process_info  = GetModelPart().GetProcessInfo();
        const int number_of_particles = (int)mListOfSphericContinuumParticles.size();

        #pragma omp parallel 
        {
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++){
                mListOfSphericContinuumParticles[i]->SetInitialSphereContacts(r_process_info);
                mListOfSphericContinuumParticles[i]->CreateContinuumConstitutiveLaws();
            }
            
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++){
                mListOfSphericContinuumParticles[i]->SearchSkinParticles(r_process_info); 
            }
            
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++){
                mListOfSphericContinuumParticles[i]->ContactAreaWeighting();
            }
        }

        KRATOS_CATCH("")

    } 

    
    void SetInitialFemContacts() {           
        KRATOS_TRY

        const int number_of_particles = (int)mListOfSphericContinuumParticles.size();

        #pragma omp parallel for 
        for ( int i = 0; i<number_of_particles; i++){
            mListOfSphericContinuumParticles[i]->SetInitialFemContacts();
        }          

        KRATOS_CATCH("")
    }    
    
    void FinalizeSolutionStep()
    {
        BaseType::FinalizeSolutionStep();
        FinalizeSolutionStepFEM();
        
    }    
         
    void FinalizeSolutionStepFEM() {      
        KRATOS_TRY

        ConditionsArrayType& pConditions      = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();     
        ProcessInfo& r_process_info  = GetFemModelPart().GetProcessInfo();
        Vector rhs_cond;
        vector<unsigned int> condition_partition;
        OpenMPUtils::CreatePartition(mNumberOfThreads, pConditions.size(), condition_partition);

        #pragma omp parallel for private (rhs_cond)            
        for(int k=0; k<mNumberOfThreads; k++) {
            typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
            typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];

            for (typename ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it){
                it->FinalizeSolutionStep(r_process_info);              
            }
        }

        KRATOS_CATCH("")
    }
     
        
    virtual void Add_As_Own(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element)
    {
        KRATOS_TRY
        
        mcontacts_model_part.Elements().push_back(p_contact_element);
        
        KRATOS_CATCH("")
    }
    
    //En aquest cas m'afegeixo jo al local i a la interface local corresponent amb la particio del vei ghost
    virtual void Add_As_Local(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element)
    {
        /* */
    }
    
    //I aqui m'afegeixio yo com a ghost de la particio del vei local
    virtual void Add_As_Ghost(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element)
    {
        /* */
    }
    
    virtual void Sort_Contact_Modelpart(ModelPart& mcontacts_model_part)
    {
        /* */
    }
    
    virtual void Reassign_Ids(ModelPart& mcontacts_model_part)
    {
        /* */
    }
    
    virtual ElementsArrayType& GetAllElements(ModelPart& r_model_part)
    {
        return r_model_part.Elements();
    }
    
    virtual ElementsArrayType& GetElements(ModelPart& r_model_part)
    {
        return r_model_part.GetCommunicator().LocalMesh().Elements();
    }

    protected:
    
    bool   mcontinuum_simulating_option;
    int    mFixSwitch;
    //bool   mDempackOption;
    std::vector<SphericContinuumParticle*>  mListOfSphericContinuumParticles;
    std::vector<SphericContinuumParticle*>  mListOfGhostSphericContinuumParticles;

    }; // Class ContinuumExplicitSolverStrategy

}  // namespace Kratos

#endif // KRATOS_FILENAME_H_INCLUDED  defined

