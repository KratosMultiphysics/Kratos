//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:           MSantasusana $
//   Last modified by:    $Co-Author:          JMCarbonel $
//   Date:                $Date:               April 2014 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME)
#define  KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME


/* System includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"

namespace Kratos
{

  /**@name Kratos Globals */
  /*@{ */


  /*@} */
  /**@name Type Definitions */
  /*@{ */

  /*@} */


  /**@name  Enum's */
  /*@{ */


  /*@} */
  /**@name  Functions */
  /*@{ */



  /*@} */
  /**@name Kratos Classes */
  /*@{ */

  template<class TSparseSpace,
	   class TDenseSpace //= DenseSpace<double>
	   >
  class ExplicitCentralDifferencesScheme : public Scheme<TSparseSpace,TDenseSpace>
  {

  public:
    /**@name Type Definitions */
    /*@{ */

    //typedef boost::shared_ptr< ExplicitCentralDifferencesScheme<TSparseSpace,TDenseSpace> > Pointer;
    KRATOS_CLASS_POINTER_DEFINITION( ExplicitCentralDifferencesScheme );

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    
    typedef ModelPart::ElementsContainerType ElementsArrayType;

    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    
    typedef ModelPart::NodesContainerType NodesArrayType;


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ExplicitCentralDifferencesScheme(
				     const double  rMaximumDeltaTime,
				     const double  rDeltaTimeFraction,
				     const double  rDeltaTimePredictionLevel,
				     const bool    rRayleighDamping
				     )
      : Scheme<TSparseSpace,TDenseSpace>()
    {

      mDeltaTime.PredictionLevel  = rDeltaTimePredictionLevel;

      mDeltaTime.Maximum          = rMaximumDeltaTime;

      mDeltaTime.Fraction         = rDeltaTimeFraction;

      mRayleighDamping            = rRayleighDamping;
	
      //Allocate auxiliary memory
      int NumThreads = OpenMPUtils::GetNumThreads();

      //mMatrix.M.resize(NumThreads);
      mMatrix.D.resize(NumThreads);

      mVector.v.resize(NumThreads);
      //mVector.a.resize(NumThreads);
      //mVector.ap.resize(NumThreads);
      mSchemeIsInitialized = false;

    }

    
    /** Destructor.
     */
    virtual ~ExplicitCentralDifferencesScheme() {}


    /*@} */
    /**@name Operators
     */
    /*@{ */

    virtual int Check(ModelPart& r_model_part)
    {
      KRATOS_TRY

      BaseType::Check(r_model_part);

      if(r_model_part.GetBufferSize() < 2)
      {
        KRATOS_THROW_ERROR(std::logic_error, "Insufficient buffer size for Central Difference Scheme. It has to be 2", "")
	    }

      return 0;
      KRATOS_CATCH("");
    }


    virtual void Initialize(ModelPart& r_model_part)
    {
      KRATOS_TRY

      mModelPart = r_model_part;

      if( (mDeltaTime.PredictionLevel>0) && (mSchemeIsInitialized==false) )
      {
        CalculateDeltaTime(r_model_part);
      }

      ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

      //Preparing the time values for the first step (where time = initial_time + dt)
      mTime.Current         = rCurrentProcessInfo[TIME]+rCurrentProcessInfo[DELTA_TIME];
      mTime.Delta           = rCurrentProcessInfo[DELTA_TIME];
      mTime.Middle          = mTime.Current -  0.5*mTime.Delta;
      mTime.Previous        = mTime.Current -  mTime.Delta;
      mTime.PreviousMiddle  = mTime.Current - 1.5*mTime.Delta;

      if(mSchemeIsInitialized==false)
      {
        InitializeExplicitScheme(r_model_part);
      }
      else
      {
        SchemeCustomInitialization(r_model_part);
      }

      mSchemeIsInitialized = true;

      KRATOS_CATCH("")
    }

    //***************************************************************************

    virtual void InitializeSolutionStep(ModelPart& r_model_part,
					TSystemMatrixType& A,
					TSystemVectorType& Dx,
					TSystemVectorType& b
					)
    {
      KRATOS_TRY

      BaseType::InitializeSolutionStep(r_model_part,A,Dx,b);

      if(mDeltaTime.PredictionLevel>1)
      {
        CalculateDeltaTime(r_model_part);
      }

      InitializeResidual(r_model_part);       
	
      KRATOS_CATCH("")
    }

    //**************************************************************************
    
    void InitializeResidual( ModelPart& r_model_part )

    {
        KRATOS_TRY

        NodesArrayType& pNodes   = r_model_part.Nodes();

#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

        #pragma omp parallel for
        for(int k=0; k<number_of_threads; k++)
        {
            typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
            typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

            for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {
              array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(FORCE_RESIDUAL));  
              noalias(node_rhs)             = ZeroVector(3);
            }
        }

        KRATOS_CATCH("")
    }

    //***************************************************************************

    void CalculateDeltaTime(ModelPart& r_model_part)
    {

      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
      ElementsArrayType& pElements      = r_model_part.Elements();

#ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
#else
      int number_of_threads = 1;

#endif

      vector<unsigned int> element_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

      double safety_factor = 0.5;  //most autors recommend a value near 0.80 (Belytschko - Nonlinear FE.. 2000. chap 6. pag. 315)

      Vector delta_times(number_of_threads);

      double stable_delta_time = 0.00;

      for(int i = 0; i < number_of_threads; i++)
          delta_times[i] = mDeltaTime.Maximum/safety_factor;

#pragma omp parallel for private(stable_delta_time)
      for(int k=0; k<number_of_threads; k++)
      {
        typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
        typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
        
        for(ElementsArrayType::iterator it=it_begin; it!= it_end; it++)
        {

          //get geometric and material properties
          double length   = it->GetGeometry().Length();
          double alpha    = it->GetProperties()[RAYLEIGH_ALPHA];
          double beta     = it->GetProperties()[RAYLEIGH_BETA];
          double E        = it->GetProperties()[YOUNG_MODULUS];
          double v        = it->GetProperties()[POISSON_RATIO];
          double ro       = it->GetProperties()[DENSITY];

          //compute courant criterion
          double bulk       = E/(3.0*(1.0-2.0*v));               
          double wavespeed  = sqrt(bulk/ro);
          double w          = 2.0*wavespeed/length;   //frequency

          double psi        = 0.5*(alpha/w + beta*w); //critical ratio;
          stable_delta_time = (2.0/w)*(sqrt(1.0 + psi*psi)-psi);

          if(stable_delta_time > 0.00)
          {
            if(stable_delta_time < delta_times[k])
            {
              delta_times[k] = stable_delta_time;
            }
          }

        }
      }

      stable_delta_time  = *std::min_element(delta_times.begin(), delta_times.end());
      stable_delta_time *= safety_factor;// * 0.5; //extra factor added to get an stable delta time
        
      if(stable_delta_time < mDeltaTime.Maximum)
      {
          
        rCurrentProcessInfo[DELTA_TIME] = stable_delta_time;
        
      }
      std::cout<< "  [EXPLICIT PREDICTION LEVEL 1]:(computed stable time step = "<< stable_delta_time <<" s)"<< std::endl;
      std::cout<< "  Using  = "<< rCurrentProcessInfo[DELTA_TIME] <<" s as time step DELTA_TIME)"<< std::endl;
        
      KRATOS_CATCH("")
    }

    //***************************************************************************

    void InitializeExplicitScheme(ModelPart& r_model_part)
    {
      KRATOS_TRY

      //Set velocity to zero...

      NodesArrayType& pNodes        = r_model_part.Nodes();

#ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
#else
      int number_of_threads = 1;
#endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

#pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
      {
        typename NodesArrayType::iterator i_begin = pNodes.ptr_begin()+node_partition[k];
        typename NodesArrayType::iterator i_end   = pNodes.ptr_begin()+node_partition[k+1];

        for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
        {
          array_1d<double,3>& middle_velocity       = i->FastGetSolutionStepValue(MIDDLE_VELOCITY);
          array_1d<double,3>& current_velocity      = i->FastGetSolutionStepValue(VELOCITY);
          array_1d<double,3>& current_residual      = i->FastGetSolutionStepValue(FORCE_RESIDUAL);
          array_1d<double,3>& current_displacement  = i->FastGetSolutionStepValue(DISPLACEMENT);
          
          for (unsigned int j =0; j<3; j++)
          {
            
            middle_velocity[j]      = current_velocity[j] ;
            current_residual[j]     = 0.0;
            current_displacement[j] = 0.0;
          }
        }
      }

    KRATOS_CATCH("")
	}

  /**
      Performing the update of the solution.
  */
  //***************************************************************************
virtual void Update(ModelPart& r_model_part,
      DofsArrayType& rDofSet,
      TSystemMatrixType& A,
      TSystemVectorType& Dx,
      TSystemVectorType& b
      )
    {
      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes            = r_model_part.Nodes();

      //Step Update
      mTime.Current   = rCurrentProcessInfo[TIME];  //the first step is time = initial_time ( 0.0) + delta time
      mTime.Delta     = rCurrentProcessInfo[DELTA_TIME];

      mTime.Middle    = 0.5*(mTime.Previous + mTime.Current);


  #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
  #else
      int number_of_threads = 1;
  #endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

  #pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
      {
        typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
        typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

        for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
        {

          //Current step information "N+1" (before step update).

          const double& nodal_mass                    = i->FastGetSolutionStepValue(NODAL_MASS);
          array_1d<double,3>& current_residual        = i->FastGetSolutionStepValue(FORCE_RESIDUAL);

          array_1d<double,3>& current_velocity        = i->FastGetSolutionStepValue(VELOCITY);
          array_1d<double,3>& current_displacement    = i->FastGetSolutionStepValue(DISPLACEMENT);
          array_1d<double,3>& middle_velocity         = i->FastGetSolutionStepValue(MIDDLE_VELOCITY);

          array_1d<double,3>& current_acceleration    = i->FastGetSolutionStepValue(ACCELERATION);


          //Solution of the explicit equation:
          current_acceleration = current_residual/nodal_mass;

          int DoF = 2;
          bool Fix_displ[3] = {false, false, false};

          Fix_displ[0] = (i->pGetDof(DISPLACEMENT_X))->IsFixed();
          Fix_displ[1] = (i->pGetDof(DISPLACEMENT_Y))->IsFixed();

          if( i->HasDofFor(DISPLACEMENT_Z) )
          {
            DoF = 3;
            Fix_displ[2] = (i->pGetDof(DISPLACEMENT_Z))->IsFixed();
          }

          for (int j = 0; j < DoF; j++) 
          {
              
              if (Fix_displ[j] == true) 
              {
                  
                current_acceleration[j]  = 0.0;
                middle_velocity[j]       = 0.0; 
                
              }
              
              current_velocity[j]      = middle_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * current_acceleration[j]; //+ actual_velocity;
              middle_velocity[j]       = current_velocity[j] + (mTime.Middle - mTime.Previous) * current_acceleration[j] ; 
              current_displacement[j]  = current_displacement[j] + mTime.Delta * middle_velocity[j];      
              
          }//for DoF
          
        }//for Node 

      }//parallel

      mTime.Previous = mTime.Current;
      mTime.PreviousMiddle = mTime.Middle;


      KRATOS_CATCH("")
    }

    virtual void SchemeCustomInitialization(ModelPart& r_model_part)
    {
      KRATOS_TRY

      NodesArrayType& pNodes            = r_model_part.Nodes();


  #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
  #else
      int number_of_threads = 1;
  #endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

  #pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
      {
        typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
        typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

        for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
        {

          //Current step information "N+1" (before step update).

          const double& nodal_mass                    = i->FastGetSolutionStepValue(NODAL_MASS);
          array_1d<double,3>& current_residual        = i->FastGetSolutionStepValue(FORCE_RESIDUAL);

          array_1d<double,3>& current_velocity        = i->FastGetSolutionStepValue(VELOCITY);
          array_1d<double,3>& current_displacement    = i->FastGetSolutionStepValue(DISPLACEMENT);
          array_1d<double,3>& middle_velocity         = i->FastGetSolutionStepValue(MIDDLE_VELOCITY);

          array_1d<double,3>& current_acceleration    = i->FastGetSolutionStepValue(ACCELERATION);


          //Solution of the explicit equation:
          current_acceleration = current_residual/nodal_mass;

          int DoF = 2;
          bool Fix_displ[3] = {false, false, false};

          Fix_displ[0] = (i->pGetDof(DISPLACEMENT_X))->IsFixed();
          Fix_displ[1] = (i->pGetDof(DISPLACEMENT_Y))->IsFixed();

          if( i->HasDofFor(DISPLACEMENT_Z) )
          {
            DoF = 3;
            Fix_displ[2] = (i->pGetDof(DISPLACEMENT_Z))->IsFixed();
          }

          for (int j = 0; j < DoF; j++) 
          {
              
              if (Fix_displ[j] == true) 
              {
            
                current_acceleration[j]  = 0.0;
                middle_velocity[j]       = 0.0;
            
              }
              
            middle_velocity[j]       = 0.0 + (mTime.Middle - mTime.Previous) * current_acceleration[j] ;
            current_velocity[j]      = middle_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * current_acceleration[j]; //+ actual_velocity;
            current_displacement[j]  = 0.0;
              
          }//for DoF

        }//for node

      }//parallel

      mTime.Previous = mTime.Current;
      mTime.PreviousMiddle = mTime.Middle;


      KRATOS_CATCH("")
    }


  //***************************************************************************
  //***************************************************************************

  void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
          LocalSystemVectorType& RHS_Contribution,
          Element::EquationIdVectorType& EquationId,
          ProcessInfo& rCurrentProcessInfo)
  {

    KRATOS_TRY

    int thread = OpenMPUtils::ThisThread();

    //basic operations for the element considered
    (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);


    if(mRayleighDamping)
    {
      (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread],rCurrentProcessInfo);

      AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], rCurrentProcessInfo);

    }

    //add explicit contribution of the Element Residual (RHS) to nodal Force Residual (nodal RHS)
    (rCurrentElement) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
  }

  //Elements:
  //****************************************************************************


  void AddDynamicsToRHS(
      Element::Pointer rCurrentElement,
      LocalSystemVectorType& RHS_Contribution,
      LocalSystemMatrixType& D,
      /*LocalSystemMatrixType& M,*/
      ProcessInfo& CurrentProcessInfo)
  {
    int thread = OpenMPUtils::ThisThread();

    //adding damping contribution
    if (D.size1() != 0)
    {
      GetFirstDerivativesVector(rCurrentElement, mVector.v[thread]);

      noalias(RHS_Contribution) -= prod(D, mVector.v[thread]);
    }


  }

  void GetFirstDerivativesVector(Element::Pointer rCurrentElement, Vector& rValues ) //V at time n-1/2 old
  {

    const unsigned int number_of_nodes = rCurrentElement->GetGeometry().size();
    const unsigned int dimension       = rCurrentElement->GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
  unsigned int index = i * dimension;


  rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue(MIDDLE_VELOCITY);

  rValues[index]     = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[0];
  rValues[index + 1] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[1];

  if ( dimension == 3 )
    rValues[index + 2] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[2];

      }
  }


  //Conditions:
  //****************************************************************************

  void AddDynamicsToRHS(
      Condition::Pointer rCurrentCondition,
      LocalSystemVectorType& RHS_Contribution,
      LocalSystemMatrixType& D,
      /*LocalSystemMatrixType& M,*/
      ProcessInfo& CurrentProcessInfo)
  {
    int thread = OpenMPUtils::ThisThread();

    //adding damping contribution
    if (D.size1() != 0)
    {
      GetFirstDerivativesVector(rCurrentCondition, mVector.v[thread]);

      noalias(RHS_Contribution) -= prod(D, mVector.v[thread]);
    }


  }


  void GetFirstDerivativesVector(Condition::Pointer rCurrentCondition, Vector& rValues ) //V at time n-1/2 old
  {

    const unsigned int number_of_nodes = rCurrentCondition->GetGeometry().size();
    const unsigned int dimension       = rCurrentCondition->GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size    = number_of_nodes * dimension;

    if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      unsigned int index = i * dimension;

      rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue(MIDDLE_VELOCITY);

      rValues[index]     = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[0];
      rValues[index + 1] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[1];

      if ( dimension == 3 )
        rValues[index + 2] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[2];
    }
    
  }



  //***************************************************************************
  //***************************************************************************

    /** functions totally analogous to the precedent but applied to
  the "condition" objects
    */
  virtual void Condition_Calculate_RHS_Contribution(
                Condition::Pointer rCurrentCondition,
                LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId,
                ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    //int thread = OpenMPUtils::ThisThread();

    //basic operations for the element considered
    (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
          
    //if(mRayleighDamping)
    //    {
    //         (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread],rCurrentProcessInfo);

    //         AddDynamicsToRHS (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], rCurrentProcessInfo);
    //    }


    //add explicit contribution of the Condition Residual (RHS) to nodal Force Residual (nodal RHS)
    (rCurrentCondition) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);


    KRATOS_CATCH( "" )
  }


    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

  protected:
    /**@name Protected static Member Variables */
    /*@{ */


    struct GeneralMatrices
    {
      //std::vector< Matrix > M;     //first derivative matrix  (usually mass matrix)
      std::vector< Matrix > D;     //second derivative matrix (usually damping matrix)
    };

    struct GeneralVectors
    {
      std::vector< Vector > v;    //velocity
      //std::vector< Vector > a;    //acceleration
      //std::vector< Vector > ap;   //previous acceleration      
    };


    struct DeltaTimeParameters
    {
      double PredictionLevel; // 0, 1, 2

      double Maximum;  //maximum delta time
      double Fraction; //fraction of the delta time
    }; 


    struct TimeVariables
    {
      double PreviousMiddle; //n-1/2
      double Previous;       //n
      double Middle;         //n+1/2
      double Current;        //n+1

      double Delta;          //time step
    }; 


    GeneralMatrices     mMatrix;
    GeneralVectors      mVector;
             
    bool                mSchemeIsInitialized;

    ModelPart           mModelPart;


    TimeVariables       mTime;
    DeltaTimeParameters mDeltaTime;   
    bool                mRayleighDamping;


    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    /*@} */
    /**@name Protected Operations*/
    /*@{ */

    /*@} */
    /**@name Protected  Access */
    /*@{ */

    /*@} */
    /**@name Protected Inquiry */
    /*@{ */

    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */

    /*@} */

  private:
    /**@name Static Member Variables */
    /*@{ */

    /*@} */
    /**@name Member Variables */
    /*@{ */

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */

    /*@} */
    /**@name Private  Access */
    /*@{ */

    /*@} */
    /**@name Private Inquiry */
    /*@{ */

    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /*@} */

  }; /* Class Scheme */

  /*@} */

  /**@name Type Definitions */
  /*@{ */


  /*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME  defined */

