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
                             ModelPart& contacts_model_part,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const bool MoveMeshFlag,
                             const int    delta_option,
                             const double search_tolerance,
                             const double coordination_number,
                             //const bool delta_option,
                             //const bool continuum_simulating_option,
                             typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                             typename IntegrationScheme::Pointer pScheme,
                             typename SpatialSearch::Pointer pSpSearch
      ): ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, fem_model_part, max_delta_time, n_step_search, safety_factor, MoveMeshFlag, delta_option, search_tolerance, coordination_number, p_creator_destructor, pScheme, pSpSearch), mcontacts_model_part(contacts_model_part)
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
         
        std::cout << "---------------------CONTINUUM EXPLICIT SOLVER STRATEGY-------------------------------" << "\n" <<std::endl;

        ModelPart& r_model_part             = BaseType::GetModelPart();
        ModelPart& fem_model_part           = BaseType::GetFemModelPart();
        ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
                
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

        // Omp initializations
        this->GetNumberOfThreads() = OpenMPUtils::GetNumThreads();
    
        rCurrentProcessInfo[ACTIVATE_SEARCH_VECTOR].resize(this->GetNumberOfThreads());
        this->GetNeighbourCounter().resize(this->GetNumberOfThreads());
        
        this->GetBoundingBoxOption()     = rCurrentProcessInfo[BOUNDING_BOX_OPTION];

        BaseType::InitializeSolutionStep();
        BaseType::InitializeElements();
        this->GetInitializeWasPerformed() = true;
        
        this->ApplyPrescribedBoundaryConditions();
        
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
        KRATOS_WATCH(r_model_part.GetNodalSolutionStepVariablesList())
        
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
          
          // 1. Initialize step   /////////////////////////////////            
          BaseType::InitializeSolutionStep();
          this->ContactInitializeSolutionStep();
          

          // 2. Calculate forces   /////////////////////////////////            
          BaseType::GetForce();
          //DEM_FEM..... "should be gathered into one single RHS for both particle and FEM nodes          
          Clear_forces_FEM();
          Calculate_Conditions_RHS_and_Add();
          if(mDempackOption) { this->GlobalDamping(); }       
          
           
          // 3. Move particles   /////////////////////////////////  
          BaseType::SynchronizeSolidMesh(r_model_part); // Should be... just TOTAL_FORCES (and ELASTIC_FORCES) and PARTICLE_MOMENT
          
          this->PerformTimeIntegrationOfMotion(rCurrentProcessInfo); 
           
          
          // 4. Search Neighbours /////////////////////////////////                   
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
         
          //Synch this var.
          r_model_part.GetCommunicator().MaxAll(rCurrentProcessInfo[ACTIVATE_SEARCH]);

          if(rCurrentProcessInfo[ACTIVATE_SEARCH]==1)
          {
              if ((time_step + 1)%this->GetNStepSearch() == 0 && time_step > 0){

                  if (this->GetBoundingBoxOption() == 1)
                  {
                      this->BoundingBoxUtility();                                            
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
          //DEBUG OPERATIONS
          
          //this->CheckPairWiseBreaking();
          
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
          
          /*      
           if(time_step == 1000)
           {
             this->Contact_Debug();
           }
          */
          
          //KRATOS_TIMER_STOP("SOLVE")

    } 
    
      
    void CreateContactElements() //better not to apply OMP paralelization since it is creation of spheres
    {                

        KRATOS_TRY

        typedef Node < 3 > NodeType;
        typedef Geometry<NodeType> GeometryType;

        ModelPart& r_sphere_model_part          = BaseType::GetModelPart();
        ElementsArrayType& pSphereElements      = GetAllElements(r_sphere_model_part);

        int index_new_ids = 1;
                    
        std::string ElementName;
        ElementName = std::string("ParticleContactElement");
        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);
        
        PrepareContactModelPart(r_sphere_model_part,mcontacts_model_part);

          //Here we are going to create contact elements when we are on a target particle and we see a neighbour whose id is higher than ours.
          //We create also a pointer from the node to the element, after creating it.
          //When our particle has a higher ID than the neighbour we also create a pointer to the (previously) created contact element.
          //We proceed in this way because we want to have the pointers to contact elements in a list in the same order than the initial elements order.
                 
        
        for (typename ElementsArrayType::ptr_iterator it= pSphereElements.ptr_begin(); it!=pSphereElements.ptr_end(); ++it)
        {
            ParticleWeakVectorType& r_continuum_ini_neighbours = (*it)->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);   

            size_t neighbour_index = 0;

            for(ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator = r_continuum_ini_neighbours.ptr_begin();
                continuum_ini_neighbour_iterator != r_continuum_ini_neighbours.ptr_end(); continuum_ini_neighbour_iterator++)
            {
                
                if ( (*it)->Id() < (*continuum_ini_neighbour_iterator).lock()->Id() || ContactElementsParallelCondition(it,continuum_ini_neighbour_iterator))                           //to avoid repetition
                {
 
                    Properties::Pointer properties =  mcontacts_model_part.pGetProperties(0);                   //Needed for the creation. It is arbitrary since there are non meaningful properties in this application.
                    Geometry<Node<3> >::PointsArrayType  NodeArray(2);
                    

                    if( (*it)->Id() < (*continuum_ini_neighbour_iterator).lock()->Id() )
                    {
                        NodeArray.GetContainer()[0] = (*it)->GetGeometry()(0);
                        NodeArray.GetContainer()[1] = (*continuum_ini_neighbour_iterator).lock()->GetGeometry()(0);
                    }    
                    else
                    {
                        NodeArray.GetContainer()[1] = (*it)->GetGeometry()(0);
                        NodeArray.GetContainer()[0] = (*continuum_ini_neighbour_iterator).lock()->GetGeometry()(0);          
                    }

                    Element::Pointer p_contact_element = rReferenceElement.Create(index_new_ids, NodeArray, properties);
                    Element::WeakPointer p_weak = Element::WeakPointer(p_contact_element);
                    
                     //generating the elements
                    
                    if(ContactElementsParallelCondition(it,continuum_ini_neighbour_iterator))  //false, only for MPI
                    {
                        ((*continuum_ini_neighbour_iterator).lock())->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER).push_back(p_weak);
                        (*it)->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(neighbour_index) = p_weak;
                        
                        //If ghost element is in a different partition and out local element has lower id add it as local, otherwise as ghost.
                        if( (*it)->Id() < (*continuum_ini_neighbour_iterator).lock()->Id() )
                        {
                            
                            Add_As_Local(r_sphere_model_part,mcontacts_model_part,continuum_ini_neighbour_iterator,p_contact_element);
                                                  
                        }
                        
                        else
                        {   
        
                            Add_As_Ghost(r_sphere_model_part,mcontacts_model_part,continuum_ini_neighbour_iterator,p_contact_element);
                        
                          
                        }
                    }
                    
                    else 
                    {
                      
                      Add_As_Own(r_sphere_model_part,mcontacts_model_part,continuum_ini_neighbour_iterator,p_contact_element);
                        
                      (*it)->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(neighbour_index) = p_weak;
                     
                    }
                    
                    //copiar el weak a la variable nodal punters a barres
                    
                    index_new_ids++;    
 
                } 

                neighbour_index++;             
            }
      
        } //for (ElementsArrayType::ptr_iterator it= pSphereElements.ptr_begin(); it!=pSphereElements.ptr_end(); ++it)
    
        for (typename ElementsArrayType::ptr_iterator it= pSphereElements.ptr_begin(); it!=pSphereElements.ptr_end(); ++it)
        {
            ParticleWeakVectorType& r_continuum_ini_neighbours = (*it)->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
            
            size_t neighbour_index = 0;
                 
            for(ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator = r_continuum_ini_neighbours.ptr_begin();
                   continuum_ini_neighbour_iterator != r_continuum_ini_neighbours.ptr_end(); continuum_ini_neighbour_iterator++)
            {
                int neigh_size_ini_cont_neigh = (*continuum_ini_neighbour_iterator).lock()->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS).size(); //this is the size of the initial continuum neighbours of the neighbour of the particle where we are focused on.

                if (!((*it)->Id() < (*continuum_ini_neighbour_iterator).lock()->Id() ) || ContactElementsParallelCondition(it,continuum_ini_neighbour_iterator))                   //to avoid repetition
                {   
                    int index = -1;
                     
                    for (int iii=0; iii< neigh_size_ini_cont_neigh; iii++)
                    {
                        int neigh_neigh_ID = (*continuum_ini_neighbour_iterator).lock()->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[iii];
                                               
                        if( neigh_neigh_ID == int((*it)->Id()))
                        { 
                            index = iii; //we keep the last iii of the iteration and this is the one to do pushback         
                            break; 
                        }
                    } // for each ini continuum neighbour's ini continuum neighbour.

                    if (index == -1)
                    {
                        std::cout << "Wrong index!!!!" << std::endl;
                    }
                    else
                    {
                        if(index >= int(((*continuum_ini_neighbour_iterator).lock())->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER).size()))
                        {
                            std::cout << "ERROR: " << (*it)->Id() << " " << ((*continuum_ini_neighbour_iterator).lock())->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER).size() << " " << index << std::endl;
                        }

                        (*it)->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(neighbour_index) = ((*continuum_ini_neighbour_iterator).lock())->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(index);
                    }
                } //if target id > neigh id
                
                neighbour_index++;
                
            } // for every ini continuum neighbour   
            
        } //loop over particles
        
        Sort_Contact_Modelpart(mcontacts_model_part);

        Reassign_Ids(mcontacts_model_part);

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

                (it_contact)->Calculate(MEAN_CONTACT_AREA,Output,rCurrentProcessInfo); 

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
        
        /*        
        if (mFixSwitch)
        {
            
            FixHorizontalVelocities();
            rCurrentProcessInfo[FIX_VELOCITIES_FLAG] = 0;
          
        } // mFixSwitch
        */
        
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
      ElementsArrayType& pElements      = GetElements(r_model_part);
      
      OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
      
      double Output = 0.0;
      
      #pragma omp parallel for
      
      for (int k = 0; k < this->GetNumberOfThreads(); k++){
          typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
          typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

          for (typename ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
          {
            
              (it)->Calculate(MEAN_CONTACT_AREA,Output,rCurrentProcessInfo); 
            
          } //loop over particles

      }// loop threads OpenMP
      
      KRATOS_CATCH("")

    } //Particle_Area_Calculate
 
    void SetInitialDemContacts()
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
            
              (it)->Calculate(CALCULATE_SET_INITIAL_DEM_CONTACTS,Output,rCurrentProcessInfo); 
            
          } //loop over particles

      }// loop threads OpenMP
      
      KRATOS_CATCH("")

    } //SetInitialDemContacts
    
    void SetInitialFemContacts()
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
            
              (it)->Calculate(CALCULATE_SET_INITIAL_FEM_CONTACTS,Output,rCurrentProcessInfo); 
            
          } //loop over particles

      }// loop threads OpenMP
      
      KRATOS_CATCH("")

    } //SetInitialDemContacts
 

     void ApplyPrescribedBoundaryConditions()
    {
      
      KRATOS_TRY

      ModelPart& r_model_part           = BaseType::GetModelPart();

      for (ModelPart::MeshesContainerType::iterator mesh_it = r_model_part.GetMeshes().begin(); mesh_it != r_model_part.GetMeshes().end(); ++mesh_it)
      {

          bool fix_x = bool((*mesh_it)[IMPOSED_VELOCITY_X]);
          bool fix_y = bool((*mesh_it)[IMPOSED_VELOCITY_Y]);
          bool fix_z = bool((*mesh_it)[IMPOSED_VELOCITY_Z]);
          
         if( fix_x || fix_y || fix_z )
         {
         
          double vel_x = (*mesh_it)[IMPOSED_VELOCITY_X_VALUE];
          double vel_y = (*mesh_it)[IMPOSED_VELOCITY_Y_VALUE];  
          double vel_z = (*mesh_it)[IMPOSED_VELOCITY_Z_VALUE];  

          NodesArrayType& pNodes = mesh_it->Nodes();
        
          vector<unsigned int> node_partition;
          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pNodes.size(), node_partition);

          #pragma omp parallel for
          
          for(int k=0; k<this->GetNumberOfThreads(); k++)
          {
              typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
              typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

              for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
              {

                array_1d<double, 3>& velocity = i->FastGetSolutionStepValue(VELOCITY);
                
                if(fix_x)
                {
                  velocity[0] = vel_x;
                  //unsigned int pos_x = i->FastGetSolutionStepValue(VELOCITY_X_DOF_POS);  
                  //i->GetDof(VELOCITY_X, pos_x).FixDof(); 
                  //i->GetDof(VELOCITY_X).FixDof(); 
                  i->Set(DEMFlags::FIXED_VEL_X,true);
                }
                
                if(fix_y)
                {
                  velocity[1] = vel_y;
                  //unsigned int pos_y = i->FastGetSolutionStepValue(VELOCITY_Y_DOF_POS);  
                  //i->GetDof(VELOCITY_Y, pos_y).FixDof(); 
                  //i->GetDof(VELOCITY_Y).FixDof(); 
                  i->Set(DEMFlags::FIXED_VEL_Y,true);

                }
                
                if(fix_z)
                {
                  velocity[2] = vel_z;
                  //unsigned int pos_z = i->FastGetSolutionStepValue(VELOCITY_Z_DOF_POS);  
                  //i->GetDof(VELOCITY_Z, pos_z).FixDof(); 
                  //i->GetDof(VELOCITY_Z).FixDof();
                  i->Set(DEMFlags::FIXED_VEL_Z,true);
                
                }
     
              } //loop over particles

            }// loop threads OpenMP
            
          } //if(fix_x || fix_y || fix_z)
        
      } //for each mesh
      
      KRATOS_CATCH("")

    }
    

   void CheckPairWiseBreaking()
     {
          KRATOS_TRY

          // SPHERE MODEL PART

          ModelPart& r_model_part             = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();

          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

          double dummy = 0.0;
          
          #pragma omp parallel for
          
          for (int k = 0; k < this->GetNumberOfThreads(); k++){
              typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (typename ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                    
                     it->Calculate(DUMMY_DEBUG_DOUBLE, dummy , rCurrentProcessInfo);
                    
              } // loop over particles

          } // loop threads OpenMP

        KRATOS_CATCH("")
      }
      
      void AreaDebugging(std::vector<double>& total_area_vector)
        {
            KRATOS_TRY

            // SPHERE MODEL PART

            ModelPart& r_model_part             = BaseType::GetModelPart();
            ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
            
            rCurrentProcessInfo[AREA_VERTICAL_CENTRE] = 0.0;
            
            ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();

            OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

            #pragma omp parallel for

            for (int k = 0; k < this->GetNumberOfThreads(); k++){
                typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
                typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

                for (typename ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                      
                        double partial_area = 0.0;

                        it->Calculate(LOCAL_CONTACT_AREA_HIGH, partial_area , rCurrentProcessInfo);
                        
                        total_area_vector[OpenMPUtils::ThisThread()] += partial_area;
                        
                      
                } // loop over particles

            } // loop threads OpenMP

          KRATOS_CATCH("")
        }
   
    
    //DEMFFEM
    
    void Calculate_Conditions_RHS_and_Add()
    {
      
      KRATOS_TRY
      
      ConditionsArrayType& pConditions      = BaseType::GetFemModelPart().GetCommunicator().LocalMesh().Conditions();     

      ProcessInfo& CurrentProcessInfo  = BaseType::GetFemModelPart().GetProcessInfo();

      Vector rhs_cond;

      vector<unsigned int> condition_partition;
      OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pConditions.size(), condition_partition);
      unsigned int index;
      
      #pragma omp parallel for private (index, rhs_cond)     
      
      for(int k=0; k<this->GetNumberOfThreads(); k++)
      {
          typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
          typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];

          for (typename ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
          {
              Condition::GeometryType& geom = it->GetGeometry();
              
              it->CalculateRightHandSide(rhs_cond,CurrentProcessInfo);
            
              const unsigned int& dim = geom.WorkingSpaceDimension();
              for (unsigned int i = 0; i <geom.size(); i++)
              {
                  index = i*dim;//*2;
                  array_1d<double,3>& node_rhs = geom(i)->FastGetSolutionStepValue(ELASTIC_FORCES);//TOTAL_FORCES
                 // array_1d<double,3>& node_elastic_rhs = geom(i)->FastGetSolutionStepValue(ELASTIC_FORCES);
                  
                  for(unsigned int kk=0; kk<dim; kk++)
                  {
                      geom(i)->SetLock();

                      node_rhs[kk] = node_rhs[kk] + rhs_cond[index+kk];
                      //node_elastic_rhs[kk] = node_elastic_rhs[kk] + rhs_cond[index+kk+3];
                      geom(i)->UnSetLock();
                  }
                  
              }                   
              
          }
          
      }

      KRATOS_CATCH("")
    }
    
    
    void Clear_forces_FEM()

    {
        KRATOS_TRY

        ModelPart& fem_model_part  = BaseType::GetFemModelPart();
        NodesArrayType& pNodes   = fem_model_part.Nodes();

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pNodes.size(), node_partition);

        #pragma omp parallel for
        
        for(int k=0; k<this->GetNumberOfThreads(); k++)
        {
            typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
            typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

            for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {
                array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(ELASTIC_FORCES));
                //array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(TOTAL_FORCES));
                noalias(node_rhs)             = ZeroVector(3);
            }
        }

        KRATOS_CATCH("")
    }

    
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

          for (typename ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
          {
              
            it->FinalizeSolutionStep(rCurrentProcessInfo);
              
          }
      }

      KRATOS_CATCH("")
    }

    virtual void PrepareContactModelPart(ModelPart& r_model_part, ModelPart& mcontacts_model_part)
    {  
        /* */
    }
   
    
    virtual bool ContactElementsParallelCondition(typename ElementsArrayType::ptr_iterator it, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator)
    {
        return false;
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
    

  }; // Class ContinuumExplicitSolverStrategy


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined




