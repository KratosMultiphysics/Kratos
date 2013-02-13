//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//

#if !defined(KRATOS_EXPLICIT_SOLVER_STRATEGY)
#define  KRATOS_EXPLICIT_SOLVER_STRATEGY

// /* External includes */
// #include "boost/smart_ptr.hpp"

// System includes

// Project includes
#include "utilities/timer.h"
#include "custom_utilities/neighbours_calculator.h"
#include "custom_utilities/create_and_destroy.h"

#include "custom_elements/spheric_particle.h"
#include "custom_elements/Particle_Contact_Element.h"

#include "includes/variables.h"
#include "DEM_application.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

//#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "geometries/line_3d_2.h"
#include "custom_utilities/neighbours_calculator.h"
#include "custom_strategies/schemes/integration_scheme.h"

/* REMOVE THIS */
//#include "custom_utilities/mpi_discrete_particle_configure.h" //CARLOS REVISA-HO

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

  ///@name Kratos Globals
  ///@{ 
  
  ///@} 
  ///@name Type Definitions
  ///@{ 
  
  ///@} 
  ///@name  Enum's
  ///@{
      
  ///@}
  ///@name  Functions 
  ///@{
      
  ///@}
  ///@name Kratos Classes
  ///@{
  
  /// Short class definition.
  /** Detail class definition.
  */
  template<
  class TSparseSpace,
  class TDenseSpace, 
  class TLinearSolver> 
  class ExplicitSolverStrategy : public  SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
     {
      public:
      ///@name Type Definitions
      ///@{
 
      typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
      typedef typename BaseType::TDataType TDataType;
      typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
      typedef typename BaseType::TSchemeType TSchemeType;
      typedef typename BaseType::DofsArrayType DofsArrayType;
      typedef typename Element::DofsVectorType DofsVectorType;

      typedef ModelPart::NodesContainerType NodesArrayType;
      typedef ModelPart::ElementsContainerType ElementsArrayType;
      typedef ModelPart::ConditionsContainerType ConditionsArrayType;     

      typedef ModelPart::NodesContainerType::ContainerType      NodesContainerType;
      typedef ModelPart::ElementsContainerType::ContainerType   ElementsContainerType;
      typedef ModelPart::ConditionsContainerType::ContainerType ConditionsContainerType;
      
      typedef DiscreteElement                     ParticleType;
      
      typedef WeakPointerVector<Element> ParticleWeakVectorType; 
      typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;

      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverStrategy);
 
     
      
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ExplicitSolverStrategy(){}
      
      ExplicitSolverStrategy(ModelPart& model_part,
                             ModelPart& contacts_model_part, 
                             const int dimension,
                             const double enlargement_factor,
                             const double damping_ratio,
                             const double fraction_delta_time,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const bool MoveMeshFlag,
                             const bool delta_option,
                             const bool continuum_simulating_option,
                             typename IntegrationScheme::Pointer pScheme
      ) : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag),mcontacts_model_part(contacts_model_part),mdimension(dimension) //inicialitzacio de variables const. no poden inicialitzarse a l'esquerra d'un igual. //les referencies tambe aqui
      {
          mdelta_option                = delta_option;
          mcontinuum_simulating_option = continuum_simulating_option;
          mvirtual_mass                = false;  //M: it has to be implemented.
          mElementsAreInitialized      = false;
          mConditionsAreInitialized    = false;
          mCalculateOldTime            = false;
          mSolutionStepIsInitialized   = false;
          mInitializeWasPerformed      = false;
          mComputeTime                 = false;
          mInitialConditions           = false;
          mEnlargementFactor           = enlargement_factor;
          mdamping_ratio               = damping_ratio; //not in use
          mfraction_delta_time         = fraction_delta_time;
          mmax_delta_time              = max_delta_time;
          molddelta_time               = 0.00;
          mtimestep                    = 0.00;
          mpScheme                     = pScheme;
          
          mtimestep                    = max_delta_time;
          mnstepsearch                 = n_step_search;
          msafety_factor               = safety_factor;
      }

      /// Destructor.
      virtual ~ExplicitSolverStrategy(){ 
 	Timer my_timer;
 	my_timer.PrintTimingInformation();
      }
      
      void Initialized()
      {
          KRATOS_TRY
          
#ifdef CUSTOMTIMER
 	  Timer::SetOuputFile("customTimer.time");
	  
 	  Timer::Start("INITIALIZED");
#endif
          //M: faig una primera búsqueda abans de inicialitzar elements pk allí guardaré veins inicials i altres coses.
          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();

          //1. Search Neighbours with tolerance (Not in mpi.)
          bool extension_option = false;

          if(rCurrentProcessInfo[CASE_OPTION] != 0)
          {
              extension_option = true;
          }

          //2. Initializing elements and perform the repartition
          if(mElementsAreInitialized == false)
          {
              Repart(r_model_part,1);

              InitializeElements();
          }
          
          //3. Search Neighbours with tolerance (afther first repartition process)
          SearchIniNeighbours(r_model_part,extension_option);  
          mInitializeWasPerformed   = true;

          // 4. Set Initial Contacts
          if(mdelta_option || mcontinuum_simulating_option)
          {
              Set_Initial_Contacts(mdelta_option, mcontinuum_simulating_option);  //delta option no fa falta i fer el continuu
          }

          // 5. Create the contact elements.
          //KRATOS_WATCH(rCurrentProcessInfo[CONTACT_MESH_OPTION])
          if(rCurrentProcessInfo[CONTACT_MESH_OPTION] == 1)
          {   
              CreateContactElements();
          }
          
          if(rCurrentProcessInfo[CONCRETE_TEST_OPTION] == 1)
          {
              rCurrentProcessInfo[ACTIVATE_SEARCH] = 0;
          }
         
          if( rCurrentProcessInfo[INT_DUMMY_6] )
		  {

		  FreeVelocities();
		  
		  }		    
          
                    

          //6.Final operations
          FinalizeSolutionStep();
#ifdef CUSTOMTIMER
 	  Timer::Stop("INITIALIZED");
#endif
          KRATOS_CATCH("")
      }
       
      double Solve()
      {
          KRATOS_TRY

          std::cout<<std::fixed<<std::setw(15)<<std::scientific<<std::setprecision(5);
          
          ModelPart& r_model_part          = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();
          
          bool extension_option = false;

          if(rCurrentProcessInfo[CASE_OPTION] )
          {
              extension_option = true;
          }

          int time_step = rCurrentProcessInfo[TIME_STEPS];
         
          //STRATEGY:
          
          //1.0
          InitializeSolutionStep();

          //1.1
          Synchronize(r_model_part,mcontacts_model_part);
#ifdef CUSTOMTIMER
 	  Timer::Start("SOLVEFORCE");
#endif
	  //1. Get and Calculate the forces
          GetForce();
                //C.1
                //TransferDataContactElements();
  
          //1.1. Calculate Local Dampings
          int rota_damp_id            = rCurrentProcessInfo[ROTA_DAMP_TYPE];
          int rotation_OPTION         = rCurrentProcessInfo[ROTATION_OPTION];
      
          if ( (rotation_OPTION == 1) && (rota_damp_id == 1) ) 
          {
              ApplyRotationalDampings();
          }
#ifdef CUSTOMTIMER

 	  Timer::Stop("SOLVEFORCE");

 	  Timer::Start("SOLVEMOTION");
#endif
	  //2. Motion Integration
          ComputeIntermedialVelocityAndNewDisplacement(); //llama al scheme, i aquesta ja fa el calcul dels despaçaments i tot
          
          if( time_step == 1)
          {
              mParticle_Creator_Destructor.CalculateSurroundingBoundingBox(r_model_part, mEnlargementFactor);
          }
#ifdef CUSTOMTIMER
 	  Timer::Stop("SOLVEMOTION");
        
          //***************************///
          
 	  Timer::Start("SOLVENEIGHBORS");
#endif
          //3. Neighbouring search. Every N times. +bounding box destruction
          
          if(rCurrentProcessInfo[ACTIVATE_SEARCH]==1)
          {

              if ( (time_step + 1)%mnstepsearch == 0 && time_step >0 )
              {

                  if(rCurrentProcessInfo[BOUNDING_BOX_OPTION]==1)
                  {
                      BoundingBoxUtility(mEnlargementFactor);
                  }
                  
                  SearchNeighbours(r_model_part,extension_option); //extension option false;
              }
          }
#ifdef CUSTOMTIMER
 	  Timer::Stop("SOLVENEIGHBORS");
          
          //ONLY A DEBUG TES.T

	  Timer::Start("COMPRESIVECHECK");
#endif
          if ( ( 4<3 ) && ( (time_step + 1)%150 == 0 && time_step >0 ))
          {
          
			CompressiveCheck();

			double mean = rCurrentProcessInfo[DOUBLE_DUMMY_1]/rCurrentProcessInfo[INT_DUMMY_5];
		
			rCurrentProcessInfo.SetValue(DOUBLE_DUMMY_2, mean) ;
			rCurrentProcessInfo.SetValue(DOUBLE_DUMMY_1, 0.0) ;
			rCurrentProcessInfo.SetValue(INT_DUMMY_5, 0) ;
			
		  }
#ifdef CUSTOMTIMER
	Timer::Stop("COMPRESIVECHECK");

	
	//4.Final operations
	Timer::Start("SOLVEFINAL");
#endif
        FinalizeSolutionStep();
#ifdef CUSTOMTIMER
	Timer::Stop("SOLVEFINAL");
#endif
          return 0.00;
          KRATOS_CATCH("")
      }
      
    
        
      void InitialCriticalTime()
      { 
          KRATOS_TRY

          //COMPUTE CRITICAL DELTA TIME

          if(mComputeTime==false)
          {
              ComputeCriticalTime();
              mComputeTime = true;
          }

          KRATOS_CATCH("")
      }
      
//M: A IMPLEMENTAR PER PROBLEMES STATICS
       /*
         void CalculateVirtualMass()
        {
            KRATOS_TRY

            if(mvirtual_mass == true)
            {
              ModelPart& r_model_part          = BaseType::GetModelPart();
              ElementsArrayType& pElements     = r_model_part.Elements();

              ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();

              //NodesArrayType &pNode            = r_model_part.Nodes();
              typename NodesArrayType::iterator inode;
              for(inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++)
              {
                  inode->FastGetSolutionStepValue(NODAL_MASS) = 0.0;
              }

                typename ElementsArrayType::iterator it_begin=pElements.ptr_begin();
                typename ElementsArrayType::iterator it_end=pElements.ptr_end();
                for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
                {

                   double Young  = it->GetProperties()[YOUNG_MODULUS];
                   double Length = it->GetGeometry().Length();
                   double Volume = 0.0;
                   double VirtualMass = 0.0;

                   Element::GeometryType& geom = it->GetGeometry();

                   if (geom.size() == 1)
                   {
                      VirtualMass = Young * M_PI * it->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                      if(rCurrentProcessInfo[PARTICLE_IF_CAL_ROTATE] == 1)
                      {
                          VirtualMass = VirtualMass * 2.5;
                      }
                   }
                   else if (it->GetGeometry().Dimension() == 2 && geom.size() > 2)
                   {
                       Volume = it->GetGeometry().Area();

                       VirtualMass = Young / (Length * Length) * Volume;
                   }
                   else if (it->GetGeometry().Dimension() == 3 && geom.size() > 3 )
                   {
                       Volume = it->GetGeometry().Volume();

                       VirtualMass = Young / (Length * Length) * Volume;
                   }


                    for (unsigned int i = 0; i <geom.size(); i++)
                     {
                        double& mass = geom(i)->FastGetSolutionStepValue(NODAL_MASS);
                        mass  = mass + VirtualMass / (double)geom.size();
                     }
                }
            }

            mIfHaveCalVirtualMass = true;

            KRATOS_CATCH("")
        }
      */
    
      void GetForce()
      {
          KRATOS_TRY

          //dummy variable
          //const Variable<double>& rDUMMY_FORCES = DUMMY_FORCES;
          Vector rhs_cond;
          //M: aixo es una xapuza

          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();  //M: ho necesitu aki per algoo?? per treure la tolerancia porser    
          ElementsArrayType& pElements      = GetElements(r_model_part);
          
          #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
          #endif

          vector<unsigned int> element_partition;
          OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

  //         unsigned int index = 0;

          #pragma omp parallel for //private(index)
          for(int k=0; k<number_of_threads; k++)
          {
              typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
              typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
              
              for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
              {
                  (it)->CalculateRightHandSide(rhs_cond, rCurrentProcessInfo);
                  //we use this function to call the calculate forces in general funct.
              } //loop over particles
          }// loop threads OpenMP

          KRATOS_CATCH("")
      }

      void ComputeIntermedialVelocityAndNewDisplacement()
      {
		 ModelPart& r_model_part = BaseType::GetModelPart();  
		 ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
		
		 //KRATOS_WATCH(rCurrentProcessInfo[INT_DUMMY_6])
		//KRATOS_WATCH(rCurrentProcessInfo[TIME_STEPS])
		//KRATOS_WATCH(int(rCurrentProcessInfo[DOUBLE_DUMMY_3]*rCurrentProcessInfo[FINAL_SIMULATION_TIME]/rCurrentProcessInfo[DELTA_TIME] ))
		
          if( rCurrentProcessInfo[INT_DUMMY_6] && ( rCurrentProcessInfo[TIME_STEPS] == int(0.01*rCurrentProcessInfo[DOUBLE_DUMMY_3]*rCurrentProcessInfo[FINAL_SIMULATION_TIME]/rCurrentProcessInfo[DELTA_TIME] ) ) )
		  
		  {
	
		  FixVelocities();
		  }		  
		  mpScheme->Calculate(r_model_part);
      }
      
      void ComputeCriticalTime()
      {
          KRATOS_TRY

          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements      = GetElements(r_model_part);
          
          typename ElementsArrayType::iterator it_begin = pElements.ptr_begin();
          typename ElementsArrayType::iterator it_end   = pElements.ptr_end();

          double TimeStepTemp = 0.0;
          double mfactor = 1;
              
          for(ElementsArrayType::iterator it = it_begin; it!= it_end; it++)
          {
              it->Calculate(DELTA_TIME, TimeStepTemp, rCurrentProcessInfo);

              if(mtimestep > TimeStepTemp)
              {
                  mtimestep = TimeStepTemp;
                  mfactor = msafety_factor;
              }
          }
              
          mtimestep = mtimestep/mfactor;
          
          if (mtimestep < rCurrentProcessInfo[DELTA_TIME])
          {
              rCurrentProcessInfo[DELTA_TIME] = mtimestep;
          }
      
          std::cout<<"******************Calculating TimeStep Is "<<mtimestep<<  "******************" <<std::endl;
                                
          KRATOS_CATCH("")
      }

      void InitializeSolutionStep()
      {
          KRATOS_TRY

          //SPHERE MODEL PART
                  
          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements      = GetElements(r_model_part);

          #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
          #endif

          vector<unsigned int> element_partition;
          OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

          #pragma omp parallel for //private(index)
          for(int k=0; k<number_of_threads; k++)

          {

            typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
              {

                (it)->InitializeSolutionStep(rCurrentProcessInfo); //we use this function to call the set initial contacts and the add continuum contacts.

              } //loop over particles

          }// loop threads OpenMP

       
    if ( rCurrentProcessInfo[CONTACT_MESH_OPTION] == 1 )
    {
        //CONTACT MODEL PART
                
        //ProcessInfo& rCurrentProcessInfo  = contacts_model_part.GetProcessInfo();
//      ElementsArrayType& pContactElements = GetElements(mcontacts_model_part);
        ElementsArrayType& pContactElements = GetAllElements(mcontacts_model_part);

    
        vector<unsigned int> contact_element_partition;
        //OpenMPUtils::CreatePartition(number_of_threads, pContactElements.size(), contact_element_partition);
        OpenMPUtils::CreatePartition(number_of_threads, pContactElements.size(), contact_element_partition);

        #pragma omp parallel for //private(index)
        for(int k=0; k<number_of_threads; k++)

        {

          typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
          typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];
          
          for (ElementsArrayType::iterator it_contact= it_contact_begin; it_contact!=it_contact_end; ++it_contact)
          
        {
        
          (it_contact)->InitializeSolutionStep(rCurrentProcessInfo); 

        } //loop over CONTACT ELEMENTS

        }// loop threads OpenMP

    }
    

        KRATOS_CATCH("")
      }

      void BoundingBoxUtility(double enlargement_factor)
      {
          KRATOS_TRY

          ModelPart& r_model_part = BaseType::GetModelPart();
          mParticle_Creator_Destructor.DestroyDistantParticles( r_model_part );

          KRATOS_CATCH("")
      }
      
    void CreateContactElements() //better not to apply OMP paralelization since it is creation of spheres
    {                
        KRATOS_TRY
        
        typedef Node < 3 > NodeType;
        typedef Geometry<NodeType> GeometryType;
        
        
        ModelPart& r_sphere_model_part          = BaseType::GetModelPart();
        //ProcessInfo& rCurrentProcessInfo        = r_sphere_model_part.GetProcessInfo();
        ElementsArrayType& pSphereElements      = GetAllElements(r_sphere_model_part);
        
        //ModelPart& r_contacts_model_part        = BaseType::GetModelPart(); //NOOOOOOOOOOOOOOOO
        //ElementsArrayType& pContactElements     = GetElements(r_contacts_model_part);

        //int index_new_ids = 1; //J.Cotela says it starts on 1. Is it 0?
        int index_new_ids = 1; //Using mpi we can't simply set this to 1
                    
        std::string ElementName;
        ElementName = std::string("ParticleContactElement");
        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);
        
        PrepareContactModelPart(r_sphere_model_part,mcontacts_model_part);
        
         /*
         * 
         * Here we are going to create contact elements when we are on a target particle and we see a neighbour which id is higher than us.
         * We create also a pointer from the node to the element, after creating it.
         * When our particle has a higher ID than the neighbour we also create a pointer to the (previously) created contact element.
         * We proced in this way becouse we want to have the pointers to contact elements in a list in the same order than the initial elements order.
         *
        */
        
        for (ElementsArrayType::ptr_iterator it= pSphereElements.ptr_begin(); it!=pSphereElements.ptr_end(); ++it)
        {
            ParticleWeakVectorType& r_continuum_ini_neighbours = (*it)->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);   
            size_t neighbour_index = 0;
            
            for(ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator = r_continuum_ini_neighbours.ptr_begin();
                   continuum_ini_neighbour_iterator != r_continuum_ini_neighbours.ptr_end(); continuum_ini_neighbour_iterator++)
            {
                if ( (*it)->Id() < (*continuum_ini_neighbour_iterator).lock()->Id() || ContactElementsParallelCondition(it,continuum_ini_neighbour_iterator))                           //to avoid repetition
                {
                    Properties::Pointer properties =  mcontacts_model_part.pGetProperties(0);                   // It is arbitrary since there are non meaningful properties in this application.
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
                    
  /*                  if( (*continuum_ini_neighbour_iterator).lock()->Id() == 3569 || (*continuum_ini_neighbour_iterator).lock()->Id() == 3404 ) 
                        std::cout << "Tracking neighbour " << (*it)->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX) << " -- " << (*continuum_ini_neighbour_iterator).lock()->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX) << " Not yet added with id " << (*it)->Id() << " - "  << (*continuum_ini_neighbour_iterator).lock()->Id() << std::endl;
                    
                    if( (*it)->Id() == 3569 || (*it)->Id() == 3404) 
                        std::cout << "Tracking parent " << (*it)->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX) << " -- " << (*continuum_ini_neighbour_iterator).lock()->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX)  << " Not yet added with id " << (*it)->Id() << " - "  << (*continuum_ini_neighbour_iterator).lock()->Id() << " Mask: " << (*it)->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH) << std::endl;
*/
                    //Notice that "(*it)->GetValue(PARTITION_INDEX) != currentPartition" it's an impossible situation
                    if(ContactElementsParallelCondition(it,continuum_ini_neighbour_iterator))
                    {
                        ((*continuum_ini_neighbour_iterator).lock())->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER).push_back(p_weak);
                        (*it)->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(neighbour_index) = p_weak;
                        
                        //If ghost element is in a different partition and out local element has lower id add it as local, otherwise as ghost.
                        if( (*it)->Id() < (*continuum_ini_neighbour_iterator).lock()->Id() )
                        {
       /*                     if( ((*it)->Id() == 3404 && (*continuum_ini_neighbour_iterator).lock()->Id() == 3569) || 
                                ((*it)->Id() == 3569 && (*continuum_ini_neighbour_iterator).lock()->Id() == 3404)    )
                                std::cout << "Partition " << (*it)->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX) << " Added as local with id " << (*it)->Id() << std::endl;
             */
             Add_As_Local(r_sphere_model_part,mcontacts_model_part,continuum_ini_neighbour_iterator,p_contact_element);
                        }
                        else
                        {   
           /*                 if( ((*it)->Id() == 3404 && (*continuum_ini_neighbour_iterator).lock()->Id() == 3569) || 
                                ((*it)->Id() == 3569 && (*continuum_ini_neighbour_iterator).lock()->Id() == 3404)    )
                                std::cout << "Partition " << (*it)->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX) << " Added as ghost with id " << (*it)->Id() << std::endl;
               */
               Add_As_Ghost(r_sphere_model_part,mcontacts_model_part,continuum_ini_neighbour_iterator,p_contact_element);
                        }
                    }
                    else 
                    {
          /*    
						if( ((*it)->Id() == 3404 && (*continuum_ini_neighbour_iterator).lock()->Id() == 3569) || 
                            ((*it)->Id() == 3569 && (*continuum_ini_neighbour_iterator).lock()->Id() == 3404)    )
                            std::cout << "Partition " << (*it)->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX) << " Added as own with id " << (*it)->Id() << std::endl;
              */
              Add_As_Own(r_sphere_model_part,mcontacts_model_part,continuum_ini_neighbour_iterator,p_contact_element);
                        (*it)->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(neighbour_index) = p_weak;
                    }
                    //copiar el weak a la variable nodal punters a barres
                    
                    index_new_ids++;        
                } 
                neighbour_index++;
            }
        }
            
        for (ElementsArrayType::ptr_iterator it= pSphereElements.ptr_begin(); it!=pSphereElements.ptr_end(); ++it)
        {
            ParticleWeakVectorType& r_continuum_ini_neighbours = (*it)->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
            size_t neighbour_index = 0;
                 
            for(ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator = r_continuum_ini_neighbours.ptr_begin();
                   continuum_ini_neighbour_iterator != r_continuum_ini_neighbours.ptr_end(); continuum_ini_neighbour_iterator++)
            {
                int neigh_size_ini_cont_neigh = (*continuum_ini_neighbour_iterator).lock()->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS).size(); //this is the size of the initial continuum neighbours of the neighbour of the particle where we are focused on.

                if (!((*it)->Id() < (*continuum_ini_neighbour_iterator).lock()->Id() || ContactElementsParallelCondition(it,continuum_ini_neighbour_iterator)))                   //to avoid repetition
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

      void ApplyRotationalDampings()
      {
          KRATOS_TRY

          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements      = GetElements(r_model_part);

          #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
          #endif

          vector<unsigned int> element_partition;
          OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);
          
          #pragma omp parallel for //private(index)
          for(int k=0; k<number_of_threads; k++)
          {
              typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
              typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

              double dummy = 0.0;

              for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
              {
                  
                  it->Calculate(PARTICLE_ROTATION_DAMP_RATIO, dummy, rCurrentProcessInfo);
              } //loop over particles
          }// loop threads OpenMP

          KRATOS_CATCH("")
      }//Apply local damps

      void MoveMesh()
      {
      }

      void FinalizeSolutionStep()
      {
          KRATOS_TRY

          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements      = GetElements(r_model_part);
          
          int trihedron_OPTION = rCurrentProcessInfo[TRIHEDRON_OPTION];

          #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
          #endif

          vector<unsigned int> element_partition;
          OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

          #pragma omp parallel for //private(index)

          for(int k=0; k<number_of_threads; k++)
          {
              typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
              typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

              for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
              {
                  (it)->FinalizeSolutionStep(rCurrentProcessInfo); //we use this function to call the set initial contacts and the add continuum contacts.
              
                  //Rotate trihedron

                  if (trihedron_OPTION==1)
                  {
                      array_1d<double,3> dummy(3,0.0);
                      double dummy2 =0.0;

                      (it)->Calculate(PRESSURE, dummy2, rCurrentProcessInfo);
                  }
                  
                  
              } //loop over particles

          }// loop threads OpenMP
          
     
          KRATOS_CATCH("")
          
      }
        
      void FixVelocities()
	  {
         KRATOS_TRY

         KRATOS_WATCH("")
		 KRATOS_WATCH("FIXING VELOCITIES!")
         
          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements      = GetElements(r_model_part);
          
          #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
          #endif

          vector<unsigned int> element_partition;
          OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

          #pragma omp parallel for //private(index)

          for(int k=0; k<number_of_threads; k++)
          {
              typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
              typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

              for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
              {
				
				   if(  it->GetGeometry()(0)->GetSolutionStepValue(GROUP_ID) == 1) 
				   {
                        //(it)->GetGeometry()(0)->Free(VELOCITY_Y);
						(it)->GetGeometry()(0)->Fix(VELOCITY_Y);
											
                        (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y)   = rCurrentProcessInfo[FIXED_VEL_TOP];
                        
				   }
				   
				    if(  it->GetGeometry()(0)->GetSolutionStepValue(GROUP_ID) == 2   )   
				   {
                        //(it)->GetGeometry()(0)->Free(VELOCITY_Y);
						(it)->GetGeometry()(0)->Fix(VELOCITY_Y);
											
                        (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y)   = rCurrentProcessInfo[FIXED_VEL_BOT];
                        
				   }
				   
              } //loop over particles

          }// loop threads OpenMP
          
     
          KRATOS_CATCH("")

	  }
      
      void FreeVelocities()
	  {
         KRATOS_TRY
                  
          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements      = GetElements(r_model_part);
          
          #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
          #endif

          vector<unsigned int> element_partition;
          OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

          #pragma omp parallel for //private(index)

          for(int k=0; k<number_of_threads; k++)
          {
              typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
              typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

              for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
              {
		
				
				//KRATOS_WATCH(it->GetValue(GROUP_ID))
				   if (  it->GetGeometry()(0)->GetSolutionStepValue(GROUP_ID) == 1 ) 
				   {
					     
                        (it)->GetGeometry()(0)->Free(VELOCITY_Y);
						 rCurrentProcessInfo[FIXED_VEL_TOP] = (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y); //cutre way yeah!   
						 //I only store one value for every ball in the group ID
                        (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y)   = 0.0;

				   }
				   
				    if ( it->GetGeometry()(0)->GetSolutionStepValue(GROUP_ID) == 2 ) 
				   {
					     
                        (it)->GetGeometry()(0)->Free(VELOCITY_Y);
						 rCurrentProcessInfo[FIXED_VEL_BOT] = (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y); //cutre way yeah!   
                         //I only store one value for every ball in the group ID
						 (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y)   = 0.0;
						 
				   }
				   
              } //loop over particles

          }// loop threads OpenMP
          
     
          KRATOS_CATCH("")

	  }    
      
        
      void CalculateEnergies()
      {
      }
          
    protected:

    private:

    ModelPart& mcontacts_model_part;    
    const unsigned int    mdimension;
    Particle_Creator_Destructor mParticle_Creator_Destructor;

    unsigned int    minitial_conditions_size;
    unsigned int    mcontact_conditions_size;  
    bool   mInitialCalculations;
    bool   mElementsAreInitialized;
    bool   mConditionsAreInitialized;
    bool   mCalculateOldTime;
    bool   mSolutionStepIsInitialized;
    bool   mComputeTime;
    bool   mInitializeWasPerformed;
    bool   mInitialConditions;
    bool   mdelta_option;
    bool   mcontinuum_simulating_option;

    bool   mvirtual_mass;

    double mEnlargementFactor;
    double mdamping_ratio;
    double malpha_damp;
    double mbeta_damp; 
    double mfraction_delta_time;
    double mmax_delta_time;
    double molddelta_time;
    double mtimestep;
    int    mnstepsearch;
    double msafety_factor;

 
    typename IntegrationScheme::Pointer mpScheme;

     
           
    void InitializeElements()
    {
        KRATOS_TRY
        ModelPart& r_model_part             = BaseType::GetModelPart();
        ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements        = GetElements(r_model_part);
        
        int trihedron_OPTION = rCurrentProcessInfo[TRIHEDRON_OPTION];

        //Matrix MassMatrix;
        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif

        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

        #pragma omp parallel for //private(index, MassMatrix)  //M. proba de compilar sense mass matrix??
        for(int k=0; k<number_of_threads; k++)
        {
          typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
          typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
          for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
          {
            //  Element::GeometryType& geom = it->GetGeometry(); ///WARNING: COMMENTED AVOIDING WARNING COMPILATION
              (it)->Initialize();
              //(it)->MassMatrix(MassMatrix, rCurrentProcessInfo); //NELSON: fa falta????
              
              // 4. Set the Local Initial Axes for the trihedron Option
              if (trihedron_OPTION==1)
              {
                double dummy =0.0;

                (it)->Calculate(DUMMY_LOCAL_AXES, dummy, rCurrentProcessInfo);

              }
            
            }
        }

        //r_model_part.GetCommunicator().AssembleCurrentData(NODAL_MASS);
        mElementsAreInitialized   = true;
        KRATOS_CATCH("")
    }
      
       
    void Set_Initial_Contacts(const bool& delta_OPTION, const bool& continuum_simulating_OPTION)
    {       
        KRATOS_TRY

        ModelPart& r_model_part             = BaseType::GetModelPart();
        ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();  //M: ho necesitu aki per algoo?? per treure la tolerancia porser
        ElementsArrayType& pElements        = GetElements(r_model_part);
        
        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif

        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

        //unsigned int index = 0;

        #pragma omp parallel for //private(index)
        for(int k=0; k<number_of_threads; k++)

        {

            typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            
                   
            for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
            {   
                            
                (it)->InitializeSolutionStep(rCurrentProcessInfo); //we use this function to call the set initial contacts and the add continuum contacts.
                                          
            } //loop over particles

        }// loop threads OpenMP

        //modifying a switch
        rCurrentProcessInfo.SetValue(NEIGH_INITIALIZED,1);
     
  
        KRATOS_CATCH("")
    }  //Set_Initial_Contacts
    
    
  
   void CompressiveCheck()
      {
          KRATOS_TRY

          ModelPart& r_model_part             = BaseType::GetModelPart();
        ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
          
		  //CONTACT MODEL PART
                
        ElementsArrayType& pContactElements = GetAllElements(mcontacts_model_part);
		
		
	#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif
    
        vector<unsigned int> contact_element_partition;
        //OpenMPUtils::CreatePartition(number_of_threads, pContactElements.size(), contact_element_partition);
        OpenMPUtils::CreatePartition(number_of_threads, pContactElements.size(), contact_element_partition);

        #pragma omp parallel for //private(index)
        for(int k=0; k<number_of_threads; k++)

        {

          typename ElementsArrayType::iterator it_contact_begin=pContactElements.ptr_begin()+contact_element_partition[k];
          typename ElementsArrayType::iterator it_contact_end=pContactElements.ptr_begin()+contact_element_partition[k+1];
          
          for (ElementsArrayType::iterator it_contact= it_contact_begin; it_contact!=it_contact_end; ++it_contact)
          
        {
        
          (it_contact)->FinalizeSolutionStep(rCurrentProcessInfo); //calculates the sigma mean.
          
		  
		   
        } //loop over CONTACT ELEMENTS

        }// loop threads OpenMP

        KRATOS_CATCH("")
		
    } //CompressiveCheck
    
    
    
    
    /* All functions below need to be redefined in the mpi specialization */
    virtual void Synchronize(ModelPart& r_model_part,ModelPart& r_contact_model_part)
    {
        /* */
    }
    
    virtual void Repart(ModelPart& r_model_part, int CalculateBoundry)
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

    virtual void SearchIniNeighbours(ModelPart& r_model_part,bool extension_option)
    { 
        typedef Neighbours_Calculator<ParticleType> NeighboursCalculatorType;
              
        NeighboursCalculatorType neighbourCalc;
        neighbourCalc.Search_Ini_Neighbours(r_model_part, extension_option);
		
    }//SearchIniNeighbours


    virtual void SearchNeighbours(ModelPart& r_model_part,bool extension_option)
    {
        typedef Neighbours_Calculator<ParticleType> NeighboursCalculatorType;
              
        NeighboursCalculatorType neighbourCalc;
        neighbourCalc.Search_Neighbours(r_model_part, extension_option);
    }//SearchNeighbours
    
    virtual void PrepareContactModelPart(ModelPart& r_model_part, ModelPart& mcontacts_model_part)
    {  
        /* */
    }
    
    virtual bool ContactElementsParallelCondition(ElementsArrayType::ptr_iterator it, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator)
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

  
  }; // Class ExplicitSolverStrategy  


        
 /*
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
                    ExplicitSolverStrategy& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const ExplicitSolverStrategy& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
    */
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 




