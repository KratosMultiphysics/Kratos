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
#include "custom_utilities/spheric_particle_hertzian.h"
#include "custom_elements/spheric_particle.h" //M: le afegit jo.. no hi era. cal que hi sigui oi???
#include "custom_elements/spheric_particle.h" //M: le afegit jo.. no hi era. cal que hi sigui oi???
#include "includes/variables.h"

/* System includes */
#include <limits>
#include<iostream>
#include<iomanip>
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
  enum Constraint_Enforcement{Penalty_Methods, Lagrange_Multiplier_Methods};  
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
	  
      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverStrategy);
 
      
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
       ExplicitSolverStrategy(){}
      
       ExplicitSolverStrategy(ModelPart& model_part,   
			      const int        dimension,
			      const double     damping_ratio,       ///para calcular la matriz de amortiguamiento proporcional a la masa
                              const double     fraction_delta_time,
                              const double     max_delta_time,
		              const bool       MoveMeshFlag,
                              const bool       delta_option,
                              const bool       continuum_simulating_option,
		              typename         IntegrationScheme::Pointer pScheme
			)
			: SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag),
                        
                            mdimension(dimension) //inicialitzacio de variables const. no poden inicialitzarse a l'esquerra d'un igual.

                        {

                           mdelta_option                = delta_option;
                           mcontinuum_simulating_option = continuum_simulating_option;
                           mvirtual_mass                = false;  //M: it needs to be implemented.
                           mElementsAreInitialized      = false;
			   mConditionsAreInitialized    = false;
			   mCalculateOldTime            = false;
                           mSolutionStepIsInitialized   = false;
			   mInitializeWasPerformed      = false;
                           mComputeTime                 = false;
			   mInitialConditions           = false;
                           mdamping_ratio               = damping_ratio;
                           mfraction_delta_time         = fraction_delta_time;
                           mmax_delta_time              = max_delta_time;
                           molddelta_time               = 0.00;
                           mtimestep                    = 0.00;
                           mpScheme                     = pScheme;
			}
      /// Destructor.
      virtual ~ExplicitSolverStrategy(){}
       
      double Solve()
      {
	KRATOS_TRY

     //    KRATOS_WATCH("NO HA FALLAT0")
        std::cout<<std::fixed<<std::setw(15)<<std::scientific<<std::setprecision(5);
        ModelPart& r_model_part              = BaseType::GetModelPart();
	ProcessInfo& CurrentProcessInfo      = r_model_part.GetProcessInfo();
	
	/*
	#ifdef _OPENMP
	double start_prod = omp_get_wtime();   
	#endif
	*/
        std::cout<<"------------------------------------------------------------------------"<<std::endl;
        std::cout<<"                 KRATOS DEM APPLICATION. TIME STEPS = "           << CurrentProcessInfo[TIME_STEPS]     <<std::endl;                                                    
	std::cout<<"------------------------------------------------------------------------"<<std::endl;
	

        //STRATEGY:

        //0.PREVIOUS OPERATIONS

        if(mComputeTime==false){
	    ComputeCriticalTime();
	    mComputeTime = true;
	}

        
        //1. Get and Calculate the forces
        GetForce();
               
        //2. Motion Integration
        ComputeIntermedialVelocityAndNewDisplacement(); //llama al scheme, i aquesta ja fa el calcul dels despaçaments i tot

       //3. Búsqueda de veins
        SearchNeighbours(r_model_part,false); //extension option false;

	
        /*
        ///Initialize solution step
	InitializeSolutionStep();

	if(mInitialConditions==false){
	   ComputeInitialConditions();
	   GetForce();
	}

        /// Computa los valores del paso de las velocidades y aceleraciones
        ComputeOldVelocitiesAndAccelerations();

	/// Actualizacion de los desplazamientos
        MoveMesh();

        ///Computing The new internal and external forces  for n+1 step

        ///Finalize solution step
	FinalizeSolutionStep();
	
	///Computing energies
	CalculateEnergies();

        #ifdef _OPENMP
	double stop_prod = omp_get_wtime();
	std::cout << "  Time solving                                 = "<<   stop_prod - start_prod    << "  seconds" << std::endl; 
	#endif
        */

        std::cout <<"FINISHED SOLVE"<<std::endl;
	return 0.00;   
	KRATOS_CATCH("")
	
      }



       /*
         void CalculateVirtualMass()
        {
            KRATOS_TRY

            if(mvirtual_mass == true)
            {
              ModelPart& r_model_part          = BaseType::GetModelPart();
              ElementsArrayType& pElements     = r_model_part.Elements();

              ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();

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
                      if(CurrentProcessInfo[PARTICLE_IF_CAL_ROTATE] == 1)
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
	
	void Initialize()
        {
            KRATOS_TRY

       //M: faig una primera búsqueda abans de inicialitzar elements pk allí guardaré veins inicials i altres coses.

        ModelPart& r_model_part           = BaseType::GetModelPart();
        //ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();

        bool extension_option = true;
        SearchNeighbours(r_model_part,extension_option);
                       
        ///Initializing elements
	  if(mElementsAreInitialized == false)
	     InitializeElements();
	  mInitializeWasPerformed   = true;

          // if (self.delta_OPTION==True or self.continuum_simulating_OPTION==True):
          if(mdelta_option || mcontinuum_simulating_option){
              Set_Initial_Contacts(mdelta_option, mcontinuum_simulating_option);  //delta option no fa falta i fer el continuu

          }

          KRATOS_CATCH("")
        }
	
	void GetForce()
	{
          KRATOS_TRY
  
          //dummy variable
          //const Variable<double>& rDUMMY_FORCES = DUMMY_FORCES;
          Vector rhs_cond;
          //M: aixo es una xapuza

          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();  //M: ho necesitu aki per algoo?? per treure la tolerancia porser
          ElementsArrayType& pElements      = r_model_part.Elements();

          #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
           #endif

          vector<unsigned int> element_partition;
          OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

          unsigned int index = 0;

          #pragma omp parallel for private(index)
          for(int k=0; k<number_of_threads; k++)

          {

            typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
              {
                     
                        //(it)->Calculate(rDUMMY_FORCES, Output, rCurrentProcessInfo);
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

        
        void  ComputeCriticalTime()
	{

        KRATOS_TRY

        ModelPart& r_model_part          = BaseType::GetModelPart();
        ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();

        if(mvirtual_mass == true)
        {
            if(mtimestep > 0.9)
            {
                mtimestep = 0.9;
            }

            std::cout<<"******************Virtual Mass TimeStep is Used******************" <<std::endl;
        }
        else
        {
              double TimeStepTemp = 0.0;

              ElementsArrayType& pElements     = r_model_part.Elements();

              typename ElementsArrayType::iterator it_begin = pElements.ptr_begin();
              typename ElementsArrayType::iterator it_end   = pElements.ptr_end();

              for(ElementsArrayType::iterator it = it_begin; it!= it_end; it++)
              {
                  it->Calculate(DELTA_TIME, TimeStepTemp, rCurrentProcessInfo);

                  if(mtimestep > TimeStepTemp)
                  {
                      mtimestep = TimeStepTemp;
                  }
              }
              double safety_factor = 0.5;
              mtimestep = safety_factor * mtimestep;

              std::cout<<"******************Real Mass TimeStep is Used******************" <<std::endl;
        }


        rCurrentProcessInfo[DELTA_TIME] = mtimestep;

        std::cout<<"******************Calculating TimeStep Is "<<mtimestep<<  "******************" <<std::endl;

        KRATOS_CATCH("")





	}

        /*

	void ComputeInitialConditions()
	{
	}

	void InitializeSolutionStep()
	{
	}
	void ComputeOldVelocitiesAndAccelerations()
	{
	}
	  
	void MoveMesh()
	{
	}

	void FinalizeSolutionStep()
	{
	}
	  
	void CalculateEnergies()
	{
	}
	*/
      
      
    protected:



    private:
            

    const unsigned int    mdimension;
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
    
    double mdamping_ratio;
    double malpha_damp;
    double mbeta_damp; 
    double mfraction_delta_time;
    double mmax_delta_time;
    double molddelta_time;
    double mtimestep;  

    typename IntegrationScheme::Pointer mpScheme;

     
           
      void InitializeElements()
      {
          KRATOS_TRY
          ModelPart& r_model_part          = BaseType::GetModelPart();
          //ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements     = r_model_part.Elements();

          //Matrix MassMatrix;
          #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
           #endif

          vector<unsigned int> element_partition;
          OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);
          //unsigned int index = 0;

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
             
             }
          }

         //r_model_part.GetCommunicator().AssembleCurrentData(NODAL_MASS);
         mElementsAreInitialized   = true;
         KRATOS_CATCH("")
        }
      
       
      void Set_Initial_Contacts(const bool& delta_OPTION, const bool& continuum_simulating_OPTION)
      {
            
          KRATOS_TRY



          ModelPart& r_model_part          = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();  //M: ho necesitu aki per algoo?? per treure la tolerancia porser
          ElementsArrayType& pElements     = r_model_part.Elements();

          #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
           #endif

          vector<unsigned int> element_partition;
          OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

          unsigned int index = 0;

          #pragma omp parallel for private(index)
          for(int k=0; k<number_of_threads; k++)

          {

            typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
              {
                //Element::GeometryType& geom = it->GetGeometry();
                
                (it)->InitializeSolutionStep(rCurrentProcessInfo); //we use this function to call the set initial contacts and the add continuum contacts.
                                           
             } //loop over particles

          }// loop threads OpenMP

        //modifying a switch

        
         
        KRATOS_CATCH("")
      }  //Set_Initial_Contacts


      void SearchNeighbours(ModelPart r_model_part,bool extension_option)
      {

        typedef DiscreteElement                                                 ParticleType;
        typedef ParticleType::Pointer                                           ParticlePointerType;
        typedef ElementsContainerType                                           ParticleContainerType;
        typedef WeakPointerVector<Element>                                      ParticleWeakVectorType;  //R:hauria de ser Discrete_Element??
        typedef typename std::vector<ParticlePointerType>                       ParticlePointerVectorType;
        typedef WeakPointerVector<Element>::iterator                            ParticleWeakIteratorType;
        typedef typename std::vector<ParticleType>::iterator                    ParticleIteratorType;
        typedef typename std::vector<ParticlePointerType>::iterator             ParticlePointerIteratorType;
        typedef std::vector<double>                                             DistanceVectorType;
        typedef std::vector<double>::iterator                                   DistanceIteratorType;

        if (mdimension == 2)
            Neighbours_Calculator<2, ParticleType>::Search_Neighbours(r_model_part,extension_option);
        else if (mdimension == 3)
            Neighbours_Calculator<3, ParticleType>::Search_Neighbours(r_model_part,extension_option);

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




