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
      
      /*
      virtual void RebuildListsOfPointersOfEachParticle(){
          
          ModelPart& r_model_part             = GetModelPart();
          ElementsArrayType& r_local_elems = r_model_part.GetCommunicator().LocalMesh().Elements();
          r_local_elems.Sort(); //if Sort() is done before parallel region, local sorts won't do anything. If this Sort() is not done, the code crashes in the next parallel region!!)
          const int number_of_particles = (int)mListOfSphericContinuumParticles.size();
          
          #pragma omp parallel for
          for ( int i=0; i<number_of_particles; i++) {   
              mListOfSphericContinuumParticles[i]->mContinuumIniNeighbourElements.resize(mListOfSphericContinuumParticles[i]->mIniNeighbourIds.size());
              for (int j=0; j<(int)mListOfSphericContinuumParticles[i]->mIniNeighbourIds.size(); j++) {
                  typename ElementsArrayType::const_iterator elem_iterator = r_local_elems.find(mListOfSphericContinuumParticles[i]->mIniNeighbourIds[j]);  //using const_iterator saves doing another Sort() inside find()
                  if (elem_iterator == r_local_elems.end()) {
                      mListOfSphericContinuumParticles[i]->mContinuumIniNeighbourElements[j] = NULL; //points to NULL if it does not find the particle with this Id
                      continue;
                  }
                  mListOfSphericContinuumParticles[i]->mContinuumIniNeighbourElements[j] = dynamic_cast<SphericContinuumParticle*>( &(*elem_iterator) );
              }  
          }                               
      }
      */
      
      void SearchNeighboursInContinuum(const bool has_mpi) {
          KRATOS_TRY
                  
          ModelPart& r_model_part           = GetModelPart();
          ProcessInfo& rcurrent_process_info  = r_model_part.GetProcessInfo();
          
          this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles); //These lists are necessary for the loop in SearchNeighbours
          this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles); 
          
          BaseType::SearchNeighbours(rcurrent_process_info[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION]); //the amplification factor has been modified after the first search.
          
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

          
          
          KRATOS_CATCH("")
      }
      
      virtual void Initialize()
      {
        
        KRATOS_TRY
         
        std::cout << "------------------CONTINUUM SOLVER STRATEGY---------------------" << "\n" <<std::endl;

        ModelPart& r_model_part             = GetModelPart();
        ModelPart& fem_model_part           = GetFemModelPart();
        ProcessInfo& rcurrent_process_info    = r_model_part.GetProcessInfo();
        
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
        
        rcurrent_process_info[SEARCH_CONTROL_VECTOR].resize(mNumberOfThreads);
        for (int i=0; i<mNumberOfThreads; i++) rcurrent_process_info[SEARCH_CONTROL_VECTOR][i] = 0;
        BaseType::GetNeighbourCounter().resize(mNumberOfThreads);

        CreatePropertiesProxies(BaseType::mFastProperties, r_model_part, *mpInlet_model_part, *mpCluster_model_part);
        
        BaseType::RepairPointersToNormalProperties(mListOfSphericParticles); 
        BaseType::RepairPointersToNormalProperties(mListOfGhostSphericParticles); 
        
        BaseType::RebuildPropertiesProxyPointers(mListOfSphericParticles);
        BaseType::RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);        
        
        //mDempackOption    = bool(rcurrent_process_info[DEMPACK_OPTION]);
        
        BaseType::GetBoundingBoxOption()     = rcurrent_process_info[BOUNDING_BOX_OPTION];
        BaseType::GetSearchControl()         = rcurrent_process_info[SEARCH_CONTROL];
              
        //BaseType::InitializeSolutionStep();  //SLS      
        BaseType::InitializeDEMElements();
        BaseType::InitializeFEMElements();                
        BaseType::InitializeSolutionStep(); //SLS
        //BaseType::ApplyPrescribedBoundaryConditions();

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
        BaseType::SearchNeighbours();
        
        if(BaseType::GetDeltaOption() == 2) {
            BaseType::SetCoordinationNumber(r_model_part);          
        }
        
        this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles);
        this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
        this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

        bool has_mpi = false;
        Check_MPI(has_mpi);
        
        if(has_mpi){
                BaseType::RepairPointersToNormalProperties(mListOfSphericParticles);
                BaseType::RepairPointersToNormalProperties(mListOfGhostSphericParticles);
        }

        BaseType::RebuildPropertiesProxyPointers(mListOfSphericParticles);
        BaseType::RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

        if(has_mpi){
            //RebuildListsOfPointersOfEachParticle(); //Serialized pointers are lost, so we rebuild them using Id's        
        }
        
        // 4. Set Initial Contacts        
        if( rcurrent_process_info[CASE_OPTION] !=0 ) {            
          SetInitialDemContacts();
        }   
        
        ComputeNewNeighboursHistoricalData();                    
        
        if(fem_model_part.Nodes().size()>0) {        
          BaseType::SearchRigidFaceNeighbours(rcurrent_process_info[LOCAL_RESOLUTION_METHOD]);
          SetInitialFemContacts();
          BaseType::ComputeNewRigidFaceNeighboursHistoricalData();        
        }
        
        if(rcurrent_process_info[CONTACT_MESH_OPTION] == 1) {   
            CreateContactElements();
            InitializeContactElements();
            ParticleAreaCalculate(true); //first time;
            ContactCalculateArea();
            ParticleAreaCalculate(false); //2nd time
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
                   
          InitializeSolutionStep();
          SearchOperations(r_model_part, has_mpi);
          BaseType::ForceOperations(r_model_part);
          BaseType::PerformTimeIntegrationOfMotion();           
          FinalizeSolutionStep();

          KRATOS_CATCH("") 
          
          return 0.0;
          
      }//Solve()
      
      void InitializeSolutionStep()
      {
          BaseType::InitializeSolutionStep();
          
      }
          //ContactInitializeSolutionStep();
      
      void SearchOperations(ModelPart& r_model_part, bool has_mpi)
      {
          ProcessInfo& rcurrent_process_info = r_model_part.GetProcessInfo();
                                  
          if(rcurrent_process_info[SEARCH_CONTROL]==0) {  
            for (int i = 0; i < mNumberOfThreads; i++) {
              if(rcurrent_process_info[SEARCH_CONTROL_VECTOR][i]==1) {                
                rcurrent_process_info[SEARCH_CONTROL]=1;                          
                std::cout << "From now on, the search is activated because some failure occurred " <<std::endl;   
                break;                
               }
             }             
          }

          const int time_step = rcurrent_process_info[TIME_STEPS];
          const double time = r_process_info[TIME];
          
            if (rcurrent_process_info[SEARCH_CONTROL] > 0) {

                if ((time_step + 1) % BaseType::GetNStepSearch() == 0 && time_step > 0) {

                    if (BaseType::GetBoundingBoxOption() == 1 && ((time >= this->GetBoundingBoxStartTime()) && (time <= this->GetBoundingBoxStopTime()))) {
                        BoundingBoxUtility();                        
                    }
                    else {
                        ParticleCreatorDestructor::Pointer& p_creator_destructor=BaseType::GetParticleCreatorDestructor();
                        p_creator_destructor->DestroyParticles(r_model_part);
                        p_creator_destructor->DestroyContactElements(*mpContact_model_part);
                    }

                    SearchNeighboursInContinuum(has_mpi); 

                    ComputeNewNeighboursHistoricalData();

                    BaseType::SearchRigidFaceNeighbours(rcurrent_process_info[LOCAL_RESOLUTION_METHOD]);
                    BaseType::ComputeNewRigidFaceNeighboursHistoricalData();
                    rcurrent_process_info[SEARCH_CONTROL] = 2;
                    if (BaseType::GetBoundingBoxOption() == 1 && has_mpi) {  //This block rebuilds all the bonds between continuum particles
                        if (rcurrent_process_info[CONTACT_MESH_OPTION] == 1) {                            
                            CreateContactElements();
                            InitializeContactElements();
                            ParticleAreaCalculate(true); //first time;
                            ContactCalculateArea();
                            ParticleAreaCalculate(false); //2nd time
                        }
                    }
                }
                else{
                    rcurrent_process_info[SEARCH_CONTROL] = 1;
                }
          }
          //Synch this var.
          r_model_part.GetCommunicator().MaxAll(rcurrent_process_info[SEARCH_CONTROL]);
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
    
    void CreateContactElements() // TODO: Re-use existing bonds, because allocating and de-allocating all of them takes a lot of time and it is not parallel!!
    {                
        KRATOS_TRY        
                
        ElementsArrayType rElements;
        ((*mpContact_model_part).Elements()).swap(rElements);
        
        int index_new_ids = 1;                    
        std::string ElementName;
        ElementName = std::string("ParticleContactElement");
        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);
        
          //Here we are going to create contact elements when we are on a target particle and we see a neighbor whose id is higher than ours.
          //We create also a pointer from the node to the element, after creating it.
          //When our particle has a higher ID than the neighbor we also create a pointer to the (previously) created contact element.
          //We proceed in this way because we want to have the pointers to contact elements in a list in the same order than the initial elements order.
                     
        Properties::Pointer properties = (*mpContact_model_part).pGetProperties(0); //Needed for the creation. It is arbitrary since there are non meaningful properties in this application.
        
        const int number_of_particles = (int)mListOfSphericContinuumParticles.size();
        
        #pragma omp parallel for 
        for (int i = 0; i < number_of_particles; i++) {
            unsigned int continuous_initial_neighbors_size = mListOfSphericContinuumParticles[i]->mContinuumInitialNeighborsSize;
            mListOfSphericContinuumParticles[i]->mBondElements.resize(continuous_initial_neighbors_size);
            for (unsigned int j=0; j<mListOfSphericContinuumParticles[i]->mBondElements.size(); j++) {
                
                mListOfSphericContinuumParticles[i]->mBondElements[j] = NULL;
            }            
        }                
        
        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {
             
            std::vector<SphericParticle*>& neighbour_elements = mListOfSphericContinuumParticles[i]->mNeighbourElements;
            unsigned int continuous_initial_neighbors_size = mListOfSphericContinuumParticles[i]->mContinuumInitialNeighborsSize;

            for (unsigned int j = 0; j < continuous_initial_neighbors_size; j++) {
                
                SphericContinuumParticle* neighbour_element = dynamic_cast<SphericContinuumParticle*>(neighbour_elements[j]);
                
                if (neighbour_element == NULL) continue; //The initial neighbor was deleted at some point in time!!

                if (mListOfSphericContinuumParticles[i]->Id() > neighbour_element->Id()) continue;                                
                
                Geometry<Node<3> >::PointsArrayType  NodeArray(2);
                
                NodeArray.GetContainer()[0] = mListOfSphericContinuumParticles[i]->GetGeometry()(0);
                NodeArray.GetContainer()[1] = neighbour_element->GetGeometry()(0);
                
                Element::Pointer p_contact_element;
                                
                #pragma omp critical
                {
                    p_contact_element = rReferenceElement.Create(index_new_ids, NodeArray, properties);
                    (*mpContact_model_part).Elements().push_back(p_contact_element);
                    index_new_ids++;
                }
                
                Element* raw_p_contact_element = p_contact_element.get();
                ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*>( raw_p_contact_element );     
                mListOfSphericContinuumParticles[i]->mBondElements[j] = p_bond;                                                                  
            }            
      
        } //for all Spheric Continuum Particles
        
        #pragma omp parallel for
        for (int i = 0; i<number_of_particles; i++) {

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
                    //if (r_continuum_ini_neighbours[j]->mContinuumIniNeighbourElements[k] == NULL) continue; //The initial neighbor was deleted at some point in time!!
                    if (neighbour_element->mNeighbourElements[k] == NULL) continue; //The initial neighbor was deleted at some point in time!!
                    
                    //if (r_continuum_ini_neighbours[j]->mContinuumIniNeighbourElements[k]->Id() == mListOfSphericContinuumParticles[i]->Id()) {
                    if (neighbour_element->mNeighbourElements[k]->Id() == mListOfSphericContinuumParticles[i]->Id()) {
                        
                        ParticleContactElement* bond = neighbour_element->mBondElements[k];
                        mListOfSphericContinuumParticles[i]->mBondElements[j] = bond; 
                        break;
                    }
                }
            }
        }                                
        
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
          ProcessInfo& rcurrent_process_info  = (*mpContact_model_part).GetProcessInfo();
         
          vector<unsigned int> contact_element_partition;
          
          OpenMPUtils::CreatePartition(mNumberOfThreads, pContactElements.size(), contact_element_partition);
          #pragma omp parallel for
          
          for (int k = 0; k < mNumberOfThreads; k++)
          {
              typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
              typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];
              
              for (typename ElementsArrayType::iterator it_contact= it_contact_begin; it_contact!=it_contact_end; ++it_contact)               
              {                             
                (it_contact)->InitializeSolutionStep(rcurrent_process_info); 
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
      
    
    void BoundingBoxUtility() {
          KRATOS_TRY
                            
          ModelPart& r_model_part = GetModelPart();
          ProcessInfo& rcurrent_process_info    = r_model_part.GetProcessInfo();
          ParticleCreatorDestructor::Pointer& p_creator_destructor=BaseType::GetParticleCreatorDestructor();
                    
          p_creator_destructor->MarkDistantParticlesForErasing(r_model_part);
          //p_creator_destructor->MarkInitialNeighboursThatAreBeingRemoved(r_model_part);
          if(rcurrent_process_info[CONTACT_MESH_OPTION] == 1)
              p_creator_destructor->MarkContactElementsForErasing(r_model_part, *mpContact_model_part);
          p_creator_destructor->DestroyParticles(r_model_part);
          p_creator_destructor->DestroyContactElements(*mpContact_model_part);

          KRATOS_CATCH("")
    }
    
    void Check_MPI(bool& has_mpi){
        VariablesList r_modelpart_nodal_variables_list = GetModelPart().GetNodalSolutionStepVariablesList();
        if(r_modelpart_nodal_variables_list.Has(PARTITION_INDEX) )  has_mpi = true;
    }
    
    void ContactCalculateArea()
    {        
        bool has_mpi = false;
        Check_MPI(has_mpi);                
        
        ElementsArrayType& pContactElements = GetAllElements(*mpContact_model_part);
              
        vector<unsigned int> contact_element_partition;
          
        OpenMPUtils::CreatePartition(mNumberOfThreads, pContactElements.size(), contact_element_partition);

        #pragma omp parallel for                   
        for(int k=0; k<mNumberOfThreads; k++) {
            typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
            typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];
              
            for (typename ElementsArrayType::iterator it= it_contact_begin; it!=it_contact_end; ++it) {
                Element* raw_p_contact_element = &(*it);
                ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*>( raw_p_contact_element );    
                p_bond->CalculateMeanContactArea(has_mpi);                
             } //loop over CONTACT ELEMENTS
          }// loop threads OpenMP
    } //ContactCalculateArea

    void ParticleAreaCalculate(const bool first) {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& rcurrent_process_info = r_model_part.GetProcessInfo();

        bool has_mpi = false;
        Check_MPI(has_mpi);

        const int number_of_particles = (int) mListOfSphericContinuumParticles.size();

        #pragma omp parallel for 
        for (int i = 0; i < number_of_particles; i++) { //Do not do this for the ghost particles!
            mListOfSphericContinuumParticles[i]->CalculateMeanContactArea(has_mpi, rcurrent_process_info, first);
        }

        KRATOS_CATCH("")
    } //ParticleAreaCalculate
 
    void SetInitialDemContacts()
    {           
      KRATOS_TRY
      
      ProcessInfo& rcurrent_process_info  = GetModelPart().GetProcessInfo();
      const int number_of_particles = (int)mListOfSphericContinuumParticles.size();

      #pragma omp parallel for
      for ( int i = 0; i<number_of_particles; i++){
          mListOfSphericContinuumParticles[i]->SetInitialSphereContacts(rcurrent_process_info);
          mListOfSphericContinuumParticles[i]->CreateContinuumConstitutiveLaws(rcurrent_process_info);
          mListOfSphericContinuumParticles[i]->ContactAreaWeighting();
      }            
           
      KRATOS_CATCH("")

    } //SetInitialDemContacts
    
    void SetInitialFemContacts() {           
      KRATOS_TRY
      
      const int number_of_particles = (int)mListOfSphericContinuumParticles.size();
      
      #pragma omp parallel for 
      for ( int i = 0; i<number_of_particles; i++){
          mListOfSphericContinuumParticles[i]->SetInitialFemContacts();
      }          
      
      KRATOS_CATCH("")
    } //SetInitialDemContacts     
    
     void FinalizeSolutionStep()
      {
        
        BaseType::FinalizeSolutionStep();
        FinalizeSolutionStepFEM();
        
    } //SetInitialDemContacts            
         
     void FinalizeSolutionStepFEM() {      
      KRATOS_TRY
      
      ConditionsArrayType& pConditions      = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();     
      ProcessInfo& rcurrent_process_info  = GetFemModelPart().GetProcessInfo();
      Vector rhs_cond;
      vector<unsigned int> condition_partition;
      OpenMPUtils::CreatePartition(mNumberOfThreads, pConditions.size(), condition_partition);
            
      #pragma omp parallel for private (rhs_cond)            
      for(int k=0; k<mNumberOfThreads; k++) {
          typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
          typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];

          for (typename ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it){
              it->FinalizeSolutionStep(rcurrent_process_info);              
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


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined




