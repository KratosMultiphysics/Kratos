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
                             ModelPart& contacts_model_part,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const bool MoveMeshFlag,
                             const bool delta_option,
                             const bool continuum_simulating_option,
                             typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                             typename IntegrationScheme::Pointer pScheme,
                             typename SpatialSearch::Pointer pSpSearch
      ): ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, max_delta_time, n_step_search, safety_factor, MoveMeshFlag, p_creator_destructor, pScheme, pSpSearch), mcontacts_model_part(contacts_model_part)
      {
          mdelta_option                 = delta_option;
          mcontinuum_simulating_option  = continuum_simulating_option;
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
        
         KRATOS_TIMER_START("INITIALIZE")
          KRATOS_WATCH("---------------------CONTINUUM EXPLICIT SOLVER STRATEGY-------------------------------")
               
          KRATOS_TRY

          ModelPart& rModelPart            = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
                
          int NumberOfElements = rModelPart.GetCommunicator().LocalMesh().ElementsArray().end() - rModelPart.GetCommunicator().LocalMesh().ElementsArray().begin();
          
          this->GetResults().resize(NumberOfElements);
          this->GetResultsDistances().resize(NumberOfElements);
          this->GetRadius().resize(NumberOfElements);
          
          // Omp initializations
          this->GetNumberOfThreads() = OpenMPUtils::GetNumThreads();
     
          // 0. Set search radius
          BaseType::SetSearchRadius(rModelPart,rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION]);
          
          // 1. Search Neighbours with tolerance (Not in mpi.)
          this->GetBoundingBoxOption()     = rCurrentProcessInfo[BOUNDING_BOX_OPTION];

          // 2. Initializing elements and perform the repartition
          if (this->GetElementsAreInitialized() == false){
              BaseType::InitializeElements();
          }
  
          this->GetInitializeWasPerformed() = true;

          // 3. Search Neighbours with tolerance (after first repartition process)
          this->SearchInitialNeighbours(rModelPart);
          
          //the search radius is modified for the next steps.
        
          double amplification = rCurrentProcessInfo[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION];
          
          ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();
          for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
          {
           
              this->GetRadius()[particle_pointer_it - pElements.begin()] *= (amplification); 
             
          }
         
          // 4. Set Initial Contacts
           BaseType::InitializeSolutionStep();

            if(rCurrentProcessInfo[CONTACT_MESH_OPTION] == 1)
          {   
              this->CreateContactElements();
              this->InitializeContactElements();
          }    

          KRATOS_CATCH("")
          
           KRATOS_TIMER_STOP("INITIALIZE")

      }// Initialize()

      virtual double Solve()
      {
            
          KRATOS_TRY
          
          KRATOS_TIMER_START("SOLVE")

          KRATOS_TIMER_START("BEGIN")
          ModelPart& rModelPart            = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

          int time_step = rCurrentProcessInfo[TIME_STEPS];
          KRATOS_TIMER_STOP("BEGIN")

          //STRATEGY:
          // 1. Here we initialize member variables that depend on the rCurrentProcessInfo
          KRATOS_TIMER_START("InitializeSolutionStep")
          BaseType::InitializeSolutionStep();
          KRATOS_TIMER_STOP("InitializeSolutionStep")
          
          KRATOS_TIMER_START("Contact_InitializeSolutionStep")
          this->Contact_InitializeSolutionStep();
          KRATOS_TIMER_STOP("Contact_InitializeSolutionStep")
          
          // 2. Get and Calculate the forces
          KRATOS_TIMER_START("GetForce")
          BaseType::GetForce();
          KRATOS_TIMER_STOP("GetForce")
           
          // 3. Motion Integration
          KRATOS_TIMER_START("PerformTimeIntegrationOfMotion")
          BaseType::PerformTimeIntegrationOfMotion(); //llama al scheme, i aquesta ja fa el calcul dels despaçaments i tot
          KRATOS_TIMER_STOP("PerformTimeIntegrationOfMotion")
          
          // 4. Synchronize  
          KRATOS_TIMER_START("SynchronizeSolidMesh")
          BaseType::SynchronizeSolidMesh(rModelPart);
          KRATOS_TIMER_STOP("SynchronizeSolidMesh")
          

          // 5. Neighbouring search. Every N times. + destruction of particles outside the bounding box
          KRATOS_TIMER_START("SearchNeighbours")
          if (rCurrentProcessInfo[ACTIVATE_SEARCH] == 1){

              if ((time_step + 1)%this->GetNStepSearch() == 0 && time_step > 0){

                  if (this->GetBoundingBoxOption() == 1){
                      BaseType::BoundingBoxUtility();
                  }

                   this->SearchNeighbours(rModelPart); //the amplification factor has benn modified after the first search.
              }

          }
          KRATOS_TIMER_STOP("SearchNeighbours")
          
          KRATOS_TIMER_START("FinalizeSolutionStep")
          BaseType::FinalizeSolutionStep();
          KRATOS_TIMER_STOP("FinalizeSolutionStep")
          
          
      KRATOS_TIMER_STOP("SOLVE")
          
          
          return 0.00;
                
          
          KRATOS_CATCH("") 
          
          
      }//Solve()

    void SearchInitialNeighbours(ModelPart& rModelPart)
    {
        KRATOS_TRY

        ModelPart& rModelPart               = BaseType::GetModelPart();
        ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();
        ProcessInfo& rCurrentProcessInfo    = rModelPart.GetProcessInfo();

        double extension = rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];
        
        for (SpatialSearch::ElementsContainerType::iterator i = pElements.begin(); i != pElements.end(); i++){
            i->GetValue(NEIGHBOUR_ELEMENTS).clear();
        }

        this->GetSpSearch()->SearchElementsInRadiusExclusive(rModelPart,this->GetRadius(),this->GetResults(),this->GetResultsDistances());

        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
        
        #pragma omp parallel for
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            size_t ResultCounter = this->GetElementPartition()[k];

            for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it,++ResultCounter){

                for (SpatialSearch::ResultElementsContainerType::iterator neighbour_it = this->GetResults()[ResultCounter].begin(); neighbour_it != this->GetResults()[ResultCounter].end(); ++neighbour_it)
                {
                  
                   double particle_radius  = (particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                   double neigh_radius     = (*neighbour_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                        
                   array_1d<double,3> other_to_me_vect      = (particle_pointer_it)->GetGeometry()(0)->Coordinates() - (*neighbour_it)->GetGeometry()(0)->Coordinates();
                   double distance                          = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                              other_to_me_vect[1] * other_to_me_vect[1] +
                                                              other_to_me_vect[2] * other_to_me_vect[2]);

                   double neighbour_search_radius           = (1.0 + extension) * neigh_radius;
                                                                       
                    if( (distance - particle_radius) < neighbour_search_radius )
                    {
                        
                      particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS).push_back(*neighbour_it);
                        
                    }
                  
                }

                this->GetResults()[ResultCounter].clear();
                this->GetResultsDistances()[ResultCounter].clear();
            }

        }

        KRATOS_CATCH("")
    } //search initial neighbours
      
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

          //Here we are going to create contact elements when we are on a target particle and we see a neighbour which id is higher than us.
          //We create also a pointer from the node to the element, after creating it.
          //When our particle has a higher ID than the neighbour we also create a pointer to the (previously) created contact element.
          //We proced in this way becouse we want to have the pointers to contact elements in a list in the same order than the initial elements order.
                 
        
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
                    
                    if(ContactElementsParallelCondition(it,continuum_ini_neighbour_iterator))
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
                    } // for each ini continuum neighbour's ini continuum neigbour.

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
        
    void Contact_InitializeSolutionStep()
    {
       
      ElementsArrayType& pContactElements = GetAllElements(mcontacts_model_part);
      
      ProcessInfo& rCurrentProcessInfo  = mcontacts_model_part.GetProcessInfo();
         
          vector<unsigned int> contact_element_partition;
          
          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pContactElements.size(), contact_element_partition);

          #pragma omp parallel for //private(index)
          
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
        return r_model_part.Elements();
    }

    protected:
    
    ModelPart& mcontacts_model_part;
    bool   mdelta_option;
    bool   mcontinuum_simulating_option;
    

  }; // Class ContinuumExplicitSolverStrategy


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined




