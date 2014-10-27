//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
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

      
      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(ContinuumExplicitSolverStrategy);

      /// Default constructor.
      ContinuumExplicitSolverStrategy(){}

      ContinuumExplicitSolverStrategy(
                             ModelPart& model_part,
                             ModelPart& fem_model_part,
                             ModelPart& cluster_model_part,
                             ModelPart& contacts_model_part,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const bool MoveMeshFlag, //TODO: is this variable used??
                             const int    delta_option,
                             const double search_tolerance,
                             const double coordination_number,
                             //const bool delta_option,
                             //const bool continuum_simulating_option,
                             typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                             typename IntegrationScheme::Pointer pScheme,
                             typename SpatialSearch::Pointer pSpSearch
      ): ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, fem_model_part, cluster_model_part, max_delta_time, n_step_search, safety_factor, MoveMeshFlag, delta_option, search_tolerance, coordination_number, p_creator_destructor, pScheme, pSpSearch), mcontacts_model_part(contacts_model_part)
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
         
        std::cout << "------------------CONTINUUM EXPLICIT SOLVER STRATEGY---------------------" << "\n" <<std::endl;

        ModelPart& r_model_part             = BaseType::GetModelPart();
        ModelPart& fem_model_part           = BaseType::GetFemModelPart();
        ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
        
        // Omp initializations
        this->GetNumberOfThreads() = OpenMPUtils::GetNumThreads();
        
        std::cout << "          **************************************************" << std::endl;
        std::cout << "            Parallelism Info:  MPI number of nodes: " << r_model_part.GetCommunicator().TotalProcesses() <<std::endl;
        if( r_model_part.GetCommunicator().TotalProcesses() > 1 )
            std::cout << "            Parallelism Info:  MPI node Id: " << r_model_part.GetCommunicator().MyPID() <<std::endl;
        std::cout << "            Parallelism Info:  OMP number of processors: " << this->GetNumberOfThreads() <<std::endl;
        std::cout << "          **************************************************" << std::endl << std::endl;
        
        this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles);
        this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
        this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), BaseType::mListOfSphericParticles);
        this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().GhostMesh().Elements(), BaseType::mListOfGhostSphericParticles);
        
        rCurrentProcessInfo[ACTIVATE_SEARCH_VECTOR].resize(this->GetNumberOfThreads());
        this->GetNeighbourCounter().resize(this->GetNumberOfThreads());

        BaseType::CreatePropertiesProxies();
        
        mDempackOption    = bool(rCurrentProcessInfo[DEMPACK_OPTION]); 
        
        unsigned int number_of_elements = r_model_part.GetCommunicator().LocalMesh().Elements().size();
                  
        if(this->GetResults().size() != number_of_elements)
          this->GetResults().resize(number_of_elements);
        
        this->GetResultsDistances().resize(number_of_elements);
                  
        if(fem_model_part.Nodes().size()>0)
        {
          this->GetRigidFaceResults().resize(number_of_elements);
          this->GetRigidFaceResultsDistances().resize(number_of_elements);
        }       
        
        this->GetBoundingBoxOption()     = rCurrentProcessInfo[BOUNDING_BOX_OPTION];

        BaseType::InitializeSolutionStep();
        
        BaseType::InitializeElements();
        BaseType::InitializeFEMElements();
        
        this->GetInitializeWasPerformed() = true;
        
        BaseType::ApplyPrescribedBoundaryConditions();
        
        // 0. Set search radius.
        
        BaseType::SetOriginalRadius(r_model_part);
        BaseType::SetSearchRadius(r_model_part, 1.0);

        // 3. Search Neighbours with tolerance (after first repartition process)
        BaseType::SearchNeighbours();
        
          if(this->GetDeltaOption() == 2)
        {
        
          BaseType::SetCoordinationNumber(r_model_part);
          
        }
        // 4. Set Initial Contacts
        
        if( rCurrentProcessInfo[CASE_OPTION] !=0 ) 
        {            
          this->SetInitialDemContacts();
        }   
        
        BaseType::ComputeNewNeighboursHistoricalData();                    
        
        if(fem_model_part.Nodes().size()>0)
        {
        
          BaseType::SearchRigidFaceNeighbours();
          this->SetInitialFemContacts();
          BaseType::ComputeNewRigidFaceNeighboursHistoricalData();
        
        }
        
        BaseType::RebuildPropertiesProxyPointers(BaseType::mListOfSphericParticles);
        BaseType::RebuildPropertiesProxyPointers(BaseType::mListOfGhostSphericParticles);
        
        //the search radius is modified for the next steps.
        BaseType::SetSearchRadius(r_model_part, rCurrentProcessInfo[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION]);

        if(rCurrentProcessInfo[CONTACT_MESH_OPTION] == 1)
          
        {   
            this->CreateContactElements();
            this->InitializeContactElements();
            rCurrentProcessInfo[AREA_CALCULATED_FLAG] = false;
            this->Particle_Area_Calculate(); //first time;
            rCurrentProcessInfo[AREA_CALCULATED_FLAG] = true;
            this->Contact_Calculate_Area();
            this->Particle_Area_Calculate(); //2nd time
        }
    
        // 5. Finalize Solution Step
        
        BaseType::FinalizeSolutionStep();
        //KRATOS_WATCH(r_model_part.GetNodalSolutionStepVariablesList())
        
        KRATOS_CATCH("")
        
        //KRATOS_TIMER_STOP("INITIALIZE")

      }// Initialize()

      virtual double Solve()
      {
            
          KRATOS_TRY
   
          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
          int time_step                     = rCurrentProcessInfo[TIME_STEPS];
          mFixSwitch                        = rCurrentProcessInfo[FIX_VELOCITIES_FLAG];
          
          this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles);
          this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
          this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), BaseType::mListOfSphericParticles);
          this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().GhostMesh().Elements(), BaseType::mListOfGhostSphericParticles);
          
          BaseType::RebuildPropertiesProxyPointers(BaseType::mListOfSphericParticles);
          BaseType::RebuildPropertiesProxyPointers(BaseType::mListOfGhostSphericParticles);
          
          // 1. Initialize step   /////////////////////////////////            
          BaseType::InitializeSolutionStep();
          
          this->ContactInitializeSolutionStep();
          
                                  
          if(rCurrentProcessInfo[ACTIVATE_SEARCH]==0)
          {            
            for (int i = 0; i < this->GetNumberOfThreads(); i++)
            {
              if(rCurrentProcessInfo[ACTIVATE_SEARCH_VECTOR][i]==1)
              {
                rCurrentProcessInfo[ACTIVATE_SEARCH]=1;
                std::cout << "From now on, the search is activated because some failure occurred " <<std::endl;   
                break;
                
               }
             }             
          }

          // 2. Calculate forces   /////////////////////////////////            
          BaseType::GetForce();
          //DEM_FEM..... "should be gathered into one single RHS for both particle and FEM nodes          
          BaseType::Clear_forces_FEM();
          BaseType::Calculate_Conditions_RHS_and_Add();
          if(mDempackOption) { this->GlobalDamping(); }       
          
           
          // 3. Move particles   /////////////////////////////////  
          BaseType::SynchronizeSolidMesh(r_model_part); // Should be... just TOTAL_FORCES (and ELASTIC_FORCES) and PARTICLE_MOMENT
          
          this->PerformTimeIntegrationOfMotion(rCurrentProcessInfo); 
                                      
          //Synch this var.
          r_model_part.GetCommunicator().MaxAll(rCurrentProcessInfo[ACTIVATE_SEARCH]);

          
          // 4. Search Neighbours /////////////////////////////////     
          if(rCurrentProcessInfo[ACTIVATE_SEARCH]==1)
          {
              if ((time_step + 1)%this->GetNStepSearch() == 0 && time_step > 0){

                  if (this->GetBoundingBoxOption() == 1)
                  {
                      this->BoundingBoxUtility();  
                      this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles);
                      this->template RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
                      this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), BaseType::mListOfSphericParticles);
                      this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().GhostMesh().Elements(), BaseType::mListOfGhostSphericParticles);
                  }                  
                   BaseType::SetSearchRadius(r_model_part,rCurrentProcessInfo[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION]);
                   BaseType::SearchNeighbours(); //the amplification factor has been modified after the first search.
                   BaseType::ComputeNewNeighboursHistoricalData();
    
                   BaseType::SetOriginalRadius(r_model_part);  
                   BaseType::SearchRigidFaceNeighbours();
                   BaseType::ComputeNewRigidFaceNeighboursHistoricalData();                   
              }

          }                          

          // 5. Finalize step   /////////////////////////////////            
          BaseType::FinalizeSolutionStep();
          FinalizeSolutionStepFEM();
          
          //this->DebugOperations(r_model_part);
		
          KRATOS_CATCH("") 
          
          return 0.0;
          
      }//Solve()

      
    void DebugOperations(ModelPart& r_model_part)
    
    {
          
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();
          ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
         
          int time_step = rCurrentProcessInfo[TIME_STEPS];
         
          if(time_step == 2)
          {
            
            std::vector<double> total_area_vector (this->GetNumberOfThreads());
            
            double total_area = 0.0;
            
            this->AreaDebugging(total_area_vector);

            for (int i=0 ; i<this->GetNumberOfThreads(); i++)
            {
              total_area += total_area_vector[i];
              
            //std::cout<<"The total measured area is: "<<total_area_vector[OpenMPUtils::ThisThread()]<<std::endl;
            }
            
            std::cout<<"The total measured area is: "<<total_area<<std::endl;
            std::cout<<" "<<std::endl;
          }          

    } 
    
    
    void CreateContactElements() //better not to apply OMP paralelization since it is creation of spheres
    {                
        KRATOS_TRY

        int index_new_ids = 1;
                    
        std::string ElementName;
        ElementName = std::string("ParticleContactElement");
        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);
        
          //Here we are going to create contact elements when we are on a target particle and we see a neighbour whose id is higher than ours.
          //We create also a pointer from the node to the element, after creating it.
          //When our particle has a higher ID than the neighbour we also create a pointer to the (previously) created contact element.
          //We proceed in this way because we want to have the pointers to contact elements in a list in the same order than the initial elements order.
                     
        Properties::Pointer properties =  mcontacts_model_part.pGetProperties(0); //Needed for the creation. It is arbitrary since there are non meaningful properties in this application.
        
        #pragma omp parallel for 
        for ( int i = 0; i<(int)mListOfSphericContinuumParticles.size(); i++){
                       
            std::vector<SphericContinuumParticle*> & r_continuum_ini_neighbours = mListOfSphericContinuumParticles[i]->mContinuumIniNeighbourElements;
            mListOfSphericContinuumParticles[i]->mBondElements.resize(r_continuum_ini_neighbours.size());
                       
            for ( unsigned int j=0; j<r_continuum_ini_neighbours.size(); j++ ) {
                
                if ( mListOfSphericContinuumParticles[i]->GetGeometry()(0)->Id() > r_continuum_ini_neighbours[j]->GetGeometry()(0)->Id() ) continue;

                if ( mListOfSphericContinuumParticles[i]->mContinuumGroup != r_continuum_ini_neighbours[j]->mContinuumGroup ) continue;                                                 
                
                Geometry<Node<3> >::PointsArrayType  NodeArray(2);
                
                NodeArray.GetContainer()[0] = mListOfSphericContinuumParticles[i]->GetGeometry()(0);
                NodeArray.GetContainer()[1] = r_continuum_ini_neighbours[j]->GetGeometry()(0);

                Element::Pointer p_contact_element = rReferenceElement.Create(index_new_ids, NodeArray, properties);
                
                #pragma omp critical
                {
                    mcontacts_model_part.Elements().push_back(p_contact_element);
                    index_new_ids++;
                }
                
                Element* raw_p_contact_element = p_contact_element.get();
                Particle_Contact_Element* p_bond = dynamic_cast<Particle_Contact_Element*>( raw_p_contact_element );     
                mListOfSphericContinuumParticles[i]->mBondElements[j] = p_bond;                                                                  
            }            
      
        } //for all Spheric Continuum Particles
        
        #pragma omp parallel for
        for ( int i = 0; i<(int)mListOfSphericContinuumParticles.size(); i++){
            
            std::vector<SphericContinuumParticle*> & r_continuum_ini_neighbours = mListOfSphericContinuumParticles[i]->mContinuumIniNeighbourElements;
            //mListOfSphericContinuumParticles[i]->mBondElements.clear();                       
            for ( unsigned int j = 0; j<r_continuum_ini_neighbours.size(); j++ ) {
                
                if ( mListOfSphericContinuumParticles[i]->GetGeometry()(0)->Id() < r_continuum_ini_neighbours[j]->GetGeometry()(0)->Id() ) continue;
                
                if ( mListOfSphericContinuumParticles[i]->mContinuumGroup != r_continuum_ini_neighbours[j]->mContinuumGroup ) continue;                
                
                for ( unsigned int k=0; k<r_continuum_ini_neighbours[j]->mContinuumIniNeighbourElements.size(); k++ )
                {                        
                    if( r_continuum_ini_neighbours[j]->mContinuumIniNeighbourElements[k]->GetGeometry()(0)->Id() == mListOfSphericContinuumParticles[i]->GetGeometry()(0)->Id() ) {
                        Particle_Contact_Element* bond = r_continuum_ini_neighbours[j]->mBondElements[k];
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
      
        ElementsArrayType& pContactElements = GetAllElements(mcontacts_model_part);
      
        vector<unsigned int> contact_element_partition;

        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pContactElements.size(), contact_element_partition);

        #pragma omp parallel for
        
        for (int k = 0; k < this->GetNumberOfThreads(); k++)
        
        {

            typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
            typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];
            
            for (typename ElementsArrayType::iterator it_contact= it_contact_begin; it_contact!=it_contact_end; ++it_contact)
            
            {
            
                (it_contact)->Initialize(); 
                
            } //loop over CONTACT ELEMENTS

        }// loop threads OpenMP

      KRATOS_CATCH("")
        
    }
   void ContactInitializeSolutionStep()
    {
       
          ElementsArrayType& pContactElements = GetAllElements(mcontacts_model_part);      
          ProcessInfo& rCurrentProcessInfo  = mcontacts_model_part.GetProcessInfo();
         
          vector<unsigned int> contact_element_partition;
          
          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pContactElements.size(), contact_element_partition);
          #pragma omp parallel for
          
          for(int k=0; k<this->GetNumberOfThreads(); k++)
          {
              typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
              typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];
              
              for (typename ElementsArrayType::iterator it_contact= it_contact_begin; it_contact!=it_contact_end; ++it_contact)               
              {                             
                (it_contact)->InitializeSolutionStep(rCurrentProcessInfo); 
              } //loop over CONTACT ELEMENTS
              
          }// loop threads OpenMP
          
   } //Contact_InitializeSolutionStep
        
   void PrepareContactElementsForPrinting()
    {
       
        ElementsArrayType& pContactElements = GetAllElements(mcontacts_model_part);
               
        vector<unsigned int> contact_element_partition;

        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pContactElements.size(), contact_element_partition);

        #pragma omp parallel for

        for(int k=0; k<this->GetNumberOfThreads(); k++)

        { 
            typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
            typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];

            for (typename ElementsArrayType::iterator it_contact= it_contact_begin; it_contact!=it_contact_end; ++it_contact)               
            {

                Element* raw_p_contact_element = &(*it_contact);
                Particle_Contact_Element* p_bond = dynamic_cast<Particle_Contact_Element*>( raw_p_contact_element );    
                p_bond->PrepareForPrinting();

            } //loop over CONTACT ELEMENTS

        }// loop threads OpenMP

    } //PrepareContactElementsForPrinting
      
    
    void BoundingBoxUtility()
      {
          KRATOS_TRY
                            
          ModelPart& r_model_part = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
          ParticleCreatorDestructor::Pointer& p_creator_destructor=BaseType::GetParticleCreatorDestructor();
                    
          p_creator_destructor->MarkDistantParticlesForErasing(r_model_part);
          if(rCurrentProcessInfo[CONTACT_MESH_OPTION] == 1)
              p_creator_destructor->MarkContactElementsForErasing(r_model_part, mcontacts_model_part);
          p_creator_destructor->DestroyParticles(r_model_part);
          p_creator_destructor->DestroyContactElements(mcontacts_model_part);

          KRATOS_CATCH("")
      }
    
    void Contact_Calculate_Area()
    {
        ModelPart& r_model_part           = BaseType::GetModelPart();
        bool has_mpi = false;
        VariablesList r_modelpart_nodal_variables_list = r_model_part.GetNodalSolutionStepVariablesList();
        if(r_modelpart_nodal_variables_list.Has(PARTITION_INDEX) )  has_mpi = true;
        
        ElementsArrayType& pContactElements = GetAllElements(mcontacts_model_part);
              
        vector<unsigned int> contact_element_partition;
          
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pContactElements.size(), contact_element_partition);

        #pragma omp parallel for         
          
        for(int k=0; k<this->GetNumberOfThreads(); k++)
        {
            typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
            typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];
              
            for (typename ElementsArrayType::iterator it= it_contact_begin; it!=it_contact_end; ++it)               
            {
                Element* raw_p_contact_element = &(*it);
                Particle_Contact_Element* p_bond = dynamic_cast<Particle_Contact_Element*>( raw_p_contact_element );    
                p_bond->CalculateMeanContactArea(has_mpi);
                
             } //loop over CONTACT ELEMENTS
          }// loop threads OpenMP
    } //Contact_calculate_Area
    
      
    void Contact_Debug()
    {
       
      ElementsArrayType& pContactElements = GetAllElements(mcontacts_model_part);
      
      ProcessInfo& rCurrentProcessInfo  = mcontacts_model_part.GetProcessInfo();
      
      double Output = 0.0;
      
          vector<unsigned int> contact_element_partition;
          
          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pContactElements.size(), contact_element_partition);

          #pragma omp parallel for          
          
          for(int k=0; k<this->GetNumberOfThreads(); k++)

          {
              typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
              typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];
              
              for (typename ElementsArrayType::iterator it_contact= it_contact_begin; it_contact!=it_contact_end; ++it_contact)               
              {

                (it_contact)->Calculate(DUMMY_DEBUG_DOUBLE,Output,rCurrentProcessInfo); 

              } //loop over CONTACT ELEMENTS

          }// loop threads OpenMP

      } //Contact_InitializeSolutionStep
         
   
     void virtual PerformTimeIntegrationOfMotion(ProcessInfo& rCurrentProcessInfo)
     {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();               
        
        BaseType::GetScheme()->Calculate(r_model_part);
        
        KRATOS_CATCH("")
      }
   
    void GlobalDamping()
        {
              
          KRATOS_TRY

          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements      = GetElements(r_model_part);
          
          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
          
          double Output = 0.0;
          
          #pragma omp parallel for          
          for (int k = 0; k < this->GetNumberOfThreads(); k++){
              typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (typename ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
              {
                
                  (it)->Calculate(DEMPACK_DAMPING,Output,rCurrentProcessInfo); 
                
              } //loop over particles

          }// loop threads OpenMP
          
          KRATOS_CATCH("")

        } //GlobalDamping

    void Particle_Area_Calculate()
    {
           
      KRATOS_TRY

      ModelPart& r_model_part           = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
      
      bool has_mpi = false;
      VariablesList r_modelpart_nodal_variables_list = r_model_part.GetNodalSolutionStepVariablesList();
      if(r_modelpart_nodal_variables_list.Has(PARTITION_INDEX) )  has_mpi = true;                        
      
      #pragma omp parallel for 
      for ( int i = 0; i<(int)mListOfSphericContinuumParticles.size(); i++){
         mListOfSphericContinuumParticles[i]->CalculateMeanContactArea(has_mpi,rCurrentProcessInfo); 
      }
      
      KRATOS_CATCH("")

    } //Particle_Area_Calculate
 
    void SetInitialDemContacts()
    {           
      KRATOS_TRY
      
      ProcessInfo& rCurrentProcessInfo  = BaseType::GetModelPart().GetProcessInfo();

      #pragma omp parallel for 
      for ( int i = 0; i<(int)mListOfSphericContinuumParticles.size(); i++){
          mListOfSphericContinuumParticles[i]->SetInitialSphereContacts(rCurrentProcessInfo);
          mListOfSphericContinuumParticles[i]->CreateContinuumConstitutiveLaws();
      }
           
      KRATOS_CATCH("")

    } //SetInitialDemContacts
    
    void SetInitialFemContacts()
    {
           
      KRATOS_TRY
      
      #pragma omp parallel for 
      for ( int i = 0; i<(int)mListOfSphericContinuumParticles.size(); i++){
          mListOfSphericContinuumParticles[i]->SetInitialFemContacts();
      }
      
      KRATOS_CATCH("")

    } //SetInitialDemContacts            
         
     void FinalizeSolutionStepFEM()
    {
      
      KRATOS_TRY
      
      ConditionsArrayType& pConditions      = BaseType::GetFemModelPart().GetCommunicator().LocalMesh().Conditions();     

      ProcessInfo& rCurrentProcessInfo  = BaseType::GetFemModelPart().GetProcessInfo();

      Vector rhs_cond;

      vector<unsigned int> condition_partition;
      OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pConditions.size(), condition_partition);
            
      #pragma omp parallel for private (rhs_cond)            
      for(int k=0; k<this->GetNumberOfThreads(); k++)
      {
          typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
          typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];

          for (typename ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it){
              it->FinalizeSolutionStep(rCurrentProcessInfo);              
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
    
    ModelPart& mcontacts_model_part;
    //bool   mdelta_option;
    bool   mcontinuum_simulating_option;
    int    mFixSwitch;
    bool   mDempackOption;
    std::vector<SphericContinuumParticle*>  mListOfSphericContinuumParticles;
    std::vector<SphericContinuumParticle*>  mListOfGhostSphericContinuumParticles;

  }; // Class ContinuumExplicitSolverStrategy


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined




