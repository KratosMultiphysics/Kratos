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

#include "custom_elements/spheric_particle.h" //M: le afegit jo.. no hi era. cal que hi sigui oi???
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

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"

#include "custom_utilities/neighbours_calculator.h"
#include "custom_strategies/schemes/integration_scheme.h"

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
      typedef Neighbours_Calculator<ParticleType> NeighboursCalculatorType;
	  
      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverStrategy);
 
      
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ExplicitSolverStrategy(){}
      
      ExplicitSolverStrategy(ModelPart& model_part,   
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
      ) : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag),mdimension(dimension) //inicialitzacio de variables const. no poden inicialitzarse a l'esquerra d'un igual.
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
      virtual ~ExplicitSolverStrategy(){}
      
      void Initialized()
      {
          KRATOS_TRY

          //M: faig una primera búsqueda abans de inicialitzar elements pk allí guardaré veins inicials i altres coses.
          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();

          //1. Search Neighbours with tolerance (Not in mpi.)
          bool extension_option = false;

          if(rCurrentProcessInfo[CASE_OPTION] != 0)
          {
              extension_option = true;
          }

          //2. Initializing elements and perform the 1st repartition
          if(mElementsAreInitialized == false)
          {
//               #ifdef _OPENMPI
//               Neighbours_Calculator<2, DiscreteElement>::parallel_partitioning(r_model_part,true);
//               for(int i = 0; i < 30; i++)
//                   Neighbours_Calculator<2, DiscreteElement>::parallel_partitioning(r_model_part,true);
//               #endif

              InitializeElements();
          }
          
          //3. Search Neighbours with tolerance (afther first repartition process)
          SearchIniNeighbours(r_model_part,extension_option);  
          mInitializeWasPerformed   = true;

          // 3. Set Initial Contacts
          if(mdelta_option || mcontinuum_simulating_option)
          {
              Set_Initial_Contacts(mdelta_option, mcontinuum_simulating_option);  //delta option no fa falta i fer el continuu
          }

          //4.Final operations
          FinalizeSolutionStep();

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
          //0.0
          Synchronize(r_model_part);
          
          //1.0
          InitializeSolutionStep();

          //1. Get and Calculate the forces
          GetForce();
  
          //1.1. Calculate Local Dampings
          int rota_damp_id            = rCurrentProcessInfo[ROTA_DAMP_TYPE];
          int rotation_OPTION         = rCurrentProcessInfo[ROTATION_OPTION];
      
          if ( (rotation_OPTION == 1) && (rota_damp_id == 1) ) 
          {
              ApplyRotationalDampings();
          }
          
          //2. Motion Integration
          ComputeIntermedialVelocityAndNewDisplacement(); //llama al scheme, i aquesta ja fa el calcul dels despaçaments i tot
          
          //3. Neighbouring search. Every N times. +bounding box destruction
          if( time_step == 1)
          {
              mParticle_Creator_Destructor.CalculateSurroundingBoundingBox(r_model_part, mEnlargementFactor);
          }
          
          if ( (time_step + 1)%mnstepsearch == 0 )
          {
              if ( (time_step + 1)%(mnstepsearch*10) == 0 )
              {
                  Repart();
              }
              if(rCurrentProcessInfo[BOUNDING_BOX_OPTION]==1)
              {
                  BoundingBoxUtility(mEnlargementFactor);
              }
              
              SearchNeighbours(r_model_part,extension_option); //extension option false;
          }
          
          //4.Final operations
          FinalizeSolutionStep();
          
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

        //modifying a switch

        KRATOS_CATCH("")
      }

      void BoundingBoxUtility(double enlargement_factor)
      {
          KRATOS_TRY

          ModelPart& r_model_part = BaseType::GetModelPart();
          mParticle_Creator_Destructor.DestroyDistantParticles( r_model_part );

          KRATOS_CATCH("")
      }

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
                  //KRATOS_WATCH(trihedron_OPTION)
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
        
      void CalculateEnergies()
      {
      }
          
    protected:

    private:

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
    
    virtual void Synchronize(ModelPart& r_model_part)
    {
        /* */
    }
    
    virtual void Repart()
    {
        /* */
    }
    
    virtual ElementsArrayType& GetElements(ModelPart& r_model_part)
    {
        return r_model_part.Elements();
    }

    virtual void SearchIniNeighbours(ModelPart& r_model_part,bool extension_option)
    { 
        //WATCH: Aixo si que es pot fer static si vols, en plan:
        // Static NeighbourCalculatorType neighbourCalc;
        NeighboursCalculatorType neighbourCalc;
        neighbourCalc.Search_Ini_Neighbours(r_model_part, extension_option);
    }//SearchIniNeighbours


    virtual void SearchNeighbours(ModelPart& r_model_part,bool extension_option)
    {
        NeighboursCalculatorType neighbourCalc;
        neighbourCalc.Search_Neighbours(r_model_part, extension_option);
    }//SearchNeighbours

  
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




