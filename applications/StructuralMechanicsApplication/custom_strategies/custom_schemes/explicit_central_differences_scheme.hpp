//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B Sautter (based on the work of MSantasusana)
//					 
//

#if !defined(KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME_HPP_INCLUDED)
#define  KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME_HPP_INCLUDED


/* System includes */


/* External includes */



/* Project includes */
#include "solving_strategies/schemes/scheme.h"



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

    ExplicitCentralDifferencesScheme(
				     const double  rMaximumDeltaTime,
				     const double  rDeltaTimeFraction,
				     const double  rDeltaTimePredictionLevel
				     )
      : Scheme<TSparseSpace,TDenseSpace>()
    {
      mDeltaTime.PredictionLevel  = rDeltaTimePredictionLevel;

      mDeltaTime.Maximum          = rMaximumDeltaTime;

      mDeltaTime.Fraction         = rDeltaTimeFraction;

      mSchemeIsInitialized = false;

    }

    
    /** Destructor.
     */
    virtual ~ExplicitCentralDifferencesScheme() {}


    /*@} */
    /**@name Operators
     */
    /*@{ */

    virtual int Check(ModelPart& rModelPart)
    {
      KRATOS_TRY

      BaseType::Check(rModelPart);

      KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2) 
        << "Insufficient buffer size for Central Difference Scheme. It has to be 2"
        << std::endl;

      return 0;
      KRATOS_CATCH("");
    }


    virtual void Initialize(ModelPart& rModelPart)
    {
      KRATOS_TRY

      if( (mDeltaTime.PredictionLevel>0) && (mSchemeIsInitialized==false) )
      {
        CalculateDeltaTime(rModelPart);
      }

      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

      //Preparing the time values for the first step (where time = initial_time + dt)
      mTime.Current         = rCurrentProcessInfo[TIME]+rCurrentProcessInfo[DELTA_TIME];
      mTime.Delta           = rCurrentProcessInfo[DELTA_TIME];
      mTime.Middle          = mTime.Current -  0.5*mTime.Delta;
      mTime.Previous        = mTime.Current -  mTime.Delta;
      mTime.PreviousMiddle  = mTime.Current - 1.5*mTime.Delta;

      if(mSchemeIsInitialized==false)
      {
        InitializeExplicitScheme(rModelPart);
      }
      else
      {
        SchemeCustomInitialization(rModelPart);
      }

      mSchemeIsInitialized = true;
      KRATOS_CATCH("")
    }

    //***************************************************************************

    virtual void InitializeSolutionStep(ModelPart& rModelPart,
					TSystemMatrixType& A,
					TSystemVectorType& Dx,
					TSystemVectorType& b)
    {
      KRATOS_TRY

      BaseType::InitializeSolutionStep(rModelPart,A,Dx,b);

      if(mDeltaTime.PredictionLevel>1)
      {
        CalculateDeltaTime(rModelPart);
      }

      InitializeResidual(rModelPart);       

	
      KRATOS_CATCH("")
    }

    //**************************************************************************
    
    void InitializeResidual( ModelPart& rModelPart )
    {
      KRATOS_TRY
       
      NodesArrayType& rNodes   = rModelPart.Nodes();

      #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
      #else
        int number_of_threads = 1;
      #endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, rNodes.size(), node_partition);

      #pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
      {
          typename NodesArrayType::iterator i_begin=rNodes.ptr_begin()+node_partition[k];
          typename NodesArrayType::iterator i_end=rNodes.ptr_begin()+node_partition[k+1];

          for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
          {
            array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(FORCE_RESIDUAL));  
            noalias(node_rhs)             = ZeroVector(3);

            if (i->HasDofFor(ROTATION_X))
            {
              array_1d<double,3>& node_rhs_moment  = (i->FastGetSolutionStepValue(MOMENT_RESIDUAL));  
              noalias(node_rhs_moment)             = ZeroVector(3);
            }
          }
      }
      KRATOS_CATCH("")
    }

    //***************************************************************************

    void CalculateDeltaTime(ModelPart& rModelPart)
    {

      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo  = rModelPart.GetProcessInfo();
      ElementsArrayType& rElements      = rModelPart.Elements();

      #ifdef _OPENMP
        const int number_of_threads = omp_get_max_threads();
      #else
        const int number_of_threads = 1;

      #endif

      vector<unsigned int> element_partition;
      OpenMPUtils::CreatePartition(number_of_threads, rElements.size(), element_partition);

      const double safety_factor = 0.5;  //most autors recommend a value near 0.80 (Belytschko - Nonlinear FE.. 2000. chap 6. pag. 315)

      std::vector<double> delta_times(number_of_threads);

      double stable_delta_time = 0.00;

      for(int i = 0; i < number_of_threads; i++)
          delta_times[i] = mDeltaTime.Maximum/safety_factor;

      #pragma omp parallel for private(stable_delta_time)

      for(int k=0; k<number_of_threads; k++)
      {
        typename ElementsArrayType::iterator it_begin=rElements.ptr_begin()+element_partition[k];
        typename ElementsArrayType::iterator it_end=rElements.ptr_begin()+element_partition[k+1];
        
        for(ElementsArrayType::iterator it=it_begin; it!= it_end; it++)
        {
          bool check_has_all_variables = true;
          double E(0.00), nu(0.00), roh(0.00), alpha(0.00), beta(0.00);
          //get geometric and material properties
          if (it->GetProperties().Has(RAYLEIGH_ALPHA))
          {
            alpha    = it->GetProperties()[RAYLEIGH_ALPHA];
          }
          if (it->GetProperties().Has(RAYLEIGH_BETA))
          {
            beta     = it->GetProperties()[RAYLEIGH_BETA];
          }
          if (it->GetProperties().Has(YOUNG_MODULUS))
          {
            E        = it->GetProperties()[YOUNG_MODULUS];
          }
          else check_has_all_variables = false;
          if (it->GetProperties().Has(POISSON_RATIO))
          {
            nu        = it->GetProperties()[POISSON_RATIO];
          }
          if (it->GetProperties().Has(DENSITY))
          {
            roh       = it->GetProperties()[DENSITY];
          }
          else check_has_all_variables = false;

          if (check_has_all_variables){
            const double length   = it->GetGeometry().Length();

            //compute courant criterion
            const double bulk_modulus       = E/(3.0*(1.0-2.0*nu));               
            const double wavespeed  = sqrt(bulk_modulus/roh);
            const double w          = 2.0*wavespeed/length;   //frequency

            const double psi        = 0.5*(alpha/w + beta*w); //critical ratio;
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
      }

      stable_delta_time  = *std::min_element(delta_times.begin(), delta_times.end());
      stable_delta_time *= safety_factor;// * 0.5; //extra factor added to get an stable delta time
        
      if(stable_delta_time < mDeltaTime.Maximum)
      {
          
        rCurrentProcessInfo[DELTA_TIME] = stable_delta_time;
        
      }

      std::cout<< "  [EXPLICIT PREDICTION LEVEL " << mDeltaTime.PredictionLevel 
        << " ] : (computed stable time step = "<< stable_delta_time <<" s)"<< std::endl;
      std::cout<< "  Using  = "<< rCurrentProcessInfo[DELTA_TIME] <<" s as time step DELTA_TIME)"<< std::endl;
        
      KRATOS_CATCH("")
    }

    //***************************************************************************

    void InitializeExplicitScheme(ModelPart& rModelPart)
    {
      KRATOS_TRY

      NodesArrayType& rNodes        = rModelPart.Nodes();

      #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
      #else
        int number_of_threads = 1;
      #endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, rNodes.size(), node_partition);

      #pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
      {
        typename NodesArrayType::iterator i_begin = rNodes.ptr_begin()+node_partition[k];
        typename NodesArrayType::iterator i_end   = rNodes.ptr_begin()+node_partition[k+1];

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

          if (i->HasDofFor(ROTATION_X))
          {
            array_1d<double,3>& middle_angular_velocity       = i->FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY);
            array_1d<double,3>& current_angular_velocity      = i->FastGetSolutionStepValue(ANGULAR_VELOCITY);
            array_1d<double,3>& current_residual_moment       = i->FastGetSolutionStepValue(MOMENT_RESIDUAL);
            array_1d<double,3>& current_rotation              = i->FastGetSolutionStepValue(ROTATION);    
            
            for (unsigned int j =0; j<3; j++)
            {
              
              middle_angular_velocity[j]      = current_angular_velocity[j] ;
              current_residual_moment[j]     = 0.0;
              current_rotation[j] = 0.0;
            }
          }



        }
      }
    KRATOS_CATCH("")
	  }

  /**
      Performing the update of the solution.
  */
  //***************************************************************************
    virtual void Update(ModelPart& rModelPart,
          DofsArrayType& rDofSet,
          TSystemMatrixType& A,
          TSystemVectorType& Dx,
          TSystemVectorType& b
          )
    {
      KRATOS_TRY
      ProcessInfo& rCurrentProcessInfo  = rModelPart.GetProcessInfo();
      NodesArrayType& rNodes            = rModelPart.Nodes();
      const double numerical_limit = std::numeric_limits<double>::epsilon();
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
      OpenMPUtils::CreatePartition(number_of_threads, rNodes.size(), node_partition);

      #pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
      {
        typename NodesArrayType::iterator i_begin=rNodes.ptr_begin()+node_partition[k];
        typename NodesArrayType::iterator i_end=rNodes.ptr_begin()+node_partition[k+1];

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
          if (nodal_mass > numerical_limit)  current_acceleration = current_residual/nodal_mass;
          else current_acceleration = ZeroVector(3);



          int DoF = 2;
          bool fix_displacements[3] = {false, false, false};

          fix_displacements[0] = (i->pGetDof(DISPLACEMENT_X))->IsFixed();
          fix_displacements[1] = (i->pGetDof(DISPLACEMENT_Y))->IsFixed();

          if( i->HasDofFor(DISPLACEMENT_Z) )
          {
            DoF = 3;
            fix_displacements[2] = (i->pGetDof(DISPLACEMENT_Z))->IsFixed();
          }

          for (int j = 0; j < DoF; j++) 
          {
              
              if (fix_displacements[j] == true) 
              {
                  
                current_acceleration[j]  = 0.0;
                middle_velocity[j]       = 0.0; 
                
              }
              
              current_velocity[j]      = middle_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * current_acceleration[j]; //+ actual_velocity;
              middle_velocity[j]       = current_velocity[j] + (mTime.Middle - mTime.Previous) * current_acceleration[j] ; 
              current_displacement[j]  = current_displacement[j] + mTime.Delta * middle_velocity[j];      
              
              
          }//for DoF


          ////// ROTATION DEGRESS OF FREEDOM
          if (i->HasDofFor(ROTATION_X))
          {
            array_1d<double,3> nodal_inertia     = i->FastGetSolutionStepValue(NODAL_INERTIA);  
            array_1d<double,3>& current_residual_moment          = i->FastGetSolutionStepValue(MOMENT_RESIDUAL);
            array_1d<double,3>& current_angular_velocity         = i->FastGetSolutionStepValue(ANGULAR_VELOCITY);
            array_1d<double,3>& current_rotation                 = i->FastGetSolutionStepValue(ROTATION);
            array_1d<double,3>& middle_angular_velocity          = i->FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY);
            array_1d<double,3>& current_angular_acceleration     = i->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

            for (int kk = 0; kk<3; ++kk)
            {         
              if (nodal_inertia[kk] > numerical_limit)  current_angular_acceleration[kk] = current_residual_moment[kk] / nodal_inertia[kk];
              else current_angular_acceleration[kk] = 0.00;
            }
            

            DoF = 2;
            bool fix_rotation[3] = {false, false, false};
            fix_rotation[0] = (i->pGetDof(ROTATION_X))->IsFixed();
            fix_rotation[1] = (i->pGetDof(ROTATION_Y))->IsFixed();   
            
            
            if (i->HasDofFor(ROTATION_Z))
            {
              DoF = 3;
              fix_rotation[1] = (i->pGetDof(ROTATION_Z))->IsFixed();  
            }

            for (int j = 0; j < DoF; j++)
            {
              if (fix_rotation[j])
              {
                current_angular_acceleration[j] = 0.00;
                middle_angular_velocity[j] = 0.00;
              }
              current_angular_velocity[j]  = middle_angular_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * current_angular_acceleration[j]; 
              middle_angular_velocity[j]   = current_angular_velocity[j] + (mTime.Middle - mTime.Previous) * current_angular_acceleration[j] ; 
              current_rotation[j]          = current_rotation[j] + mTime.Delta * middle_angular_velocity[j];
            }//for DoF
          }// Rot DoF

        }//for Node 
      }//parallel

      mTime.Previous = mTime.Current;
      mTime.PreviousMiddle = mTime.Middle;

      KRATOS_CATCH("")
    }

    virtual void SchemeCustomInitialization(ModelPart& rModelPart)
    {
      KRATOS_TRY
      NodesArrayType& rNodes            = rModelPart.Nodes();
      const double numerical_limit = std::numeric_limits<double>::epsilon();
      #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
      #else
        int number_of_threads = 1;
      #endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, rNodes.size(), node_partition);

      #pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
      {
        typename NodesArrayType::iterator i_begin=rNodes.ptr_begin()+node_partition[k];
        typename NodesArrayType::iterator i_end=rNodes.ptr_begin()+node_partition[k+1];

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
          if (nodal_mass > numerical_limit)  current_acceleration = current_residual/nodal_mass;
          else current_acceleration = ZeroVector(3);

          int DoF = 2;
          bool fix_displacements[3] = {false, false, false};

          fix_displacements[0] = (i->pGetDof(DISPLACEMENT_X))->IsFixed();
          fix_displacements[1] = (i->pGetDof(DISPLACEMENT_Y))->IsFixed();

          if( i->HasDofFor(DISPLACEMENT_Z) )
          {
            DoF = 3;
            fix_displacements[2] = (i->pGetDof(DISPLACEMENT_Z))->IsFixed();
          }

          for (int j = 0; j < DoF; j++) 
          {
              
              if (fix_displacements[j] == true) 
              {
            
                current_acceleration[j]  = 0.0;
                middle_velocity[j]       = 0.0;
            
              }
              
            middle_velocity[j]       = 0.0 + (mTime.Middle - mTime.Previous) * current_acceleration[j] ;
            current_velocity[j]      = middle_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * current_acceleration[j]; //+ actual_velocity;
            current_displacement[j]  = 0.0;

            
              
          }//for DoF
          ////// ROTATION DEGRESS OF FREEDOM
          if (i->HasDofFor(ROTATION_X))
          {
            
            array_1d<double,3> nodal_inertia     = i->FastGetSolutionStepValue(NODAL_INERTIA);  
            array_1d<double,3>& current_residual_moment          = i->FastGetSolutionStepValue(MOMENT_RESIDUAL);
            array_1d<double,3>& current_angular_velocity         = i->FastGetSolutionStepValue(ANGULAR_VELOCITY);
            array_1d<double,3>& current_rotation                 = i->FastGetSolutionStepValue(ROTATION);
            array_1d<double,3>& middle_angular_velocity          = i->FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY);
            array_1d<double,3>& current_angular_acceleration     = i->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

            
            for (int kk = 0; kk<3; ++kk)
            {         
              if (nodal_inertia[kk] > numerical_limit)  current_angular_acceleration[kk] = current_residual_moment[kk] / nodal_inertia[kk];
              else current_angular_acceleration[kk] = 0.00;
            }

            DoF = 2;
            bool fix_rotation[3] = {false, false, false};
            fix_rotation[0] = (i->pGetDof(ROTATION_X))->IsFixed();
            fix_rotation[1] = (i->pGetDof(ROTATION_Y))->IsFixed();   
            
            
            if (i->HasDofFor(ROTATION_Z))
            {
              DoF = 3;
              fix_rotation[1] = (i->pGetDof(ROTATION_Z))->IsFixed();  
            }

            for (int j = 0; j < DoF; j++)
            {
              if (fix_rotation[j])
              {
                current_angular_acceleration[j] = 0.00;
                middle_angular_velocity[j] = 0.00;
              }
               
              middle_angular_velocity[j]   = 0.00 + (mTime.Middle - mTime.Previous) * current_angular_acceleration[j] ; 
              current_angular_velocity[j]  = middle_angular_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * current_angular_acceleration[j];
              current_rotation[j]          = 0.00;
            }//for DoF
          }// Rot DoF




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

    //basic operations for the element considered
    Matrix DummyLHS;
    rCurrentElement -> CalculateLocalSystem(DummyLHS,RHS_Contribution,rCurrentProcessInfo);

    //add explicit contribution of the Element Residual (RHS) to nodal Force Residual (nodal RHS)
    (rCurrentElement) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);
    if (rCurrentElement->GetGeometry()[0].HasDofFor(ROTATION_X))
    {
      (rCurrentElement) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, MOMENT_RESIDUAL, rCurrentProcessInfo);
    }
    
    KRATOS_CATCH( "" )
  }

  //Elements:
  //****************************************************************************

  void GetFirstDerivativesVector(Element::Pointer rCurrentElement, Vector& rValues ) //V at time n-1/2 old
  {

    const unsigned int number_of_nodes = rCurrentElement->GetGeometry().size();
    const unsigned int dimension       = 3;
    unsigned int       element_size    = number_of_nodes * dimension;
    const bool check_has_rot_dof = rCurrentElement->GetGeometry()[0].HasDofFor(ROTATION_X);
    
    if (check_has_rot_dof) element_size *= 2;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int index = i * dimension;
        if (check_has_rot_dof) index *= 2;


        rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue(MIDDLE_VELOCITY);

        rValues[index]     = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[0];
        rValues[index + 1] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[1];
        rValues[index + 2] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[2];


        if (check_has_rot_dof)
        {
          rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue(MIDDLE_VELOCITY);
          
          rValues[index+dimension]     = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_ANGULAR_VELOCITY )[0];
          rValues[index+dimension+1] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_ANGULAR_VELOCITY )[1];
          rValues[index+dimension+2] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_ANGULAR_VELOCITY )[2];          
        }

      }
  }


  //Conditions:
  //****************************************************************************


  void GetFirstDerivativesVector(Condition::Pointer rCurrentCondition, Vector& rValues ) //V at time n-1/2 old
  {

    const unsigned int number_of_nodes = rCurrentCondition->GetGeometry().size();
    const unsigned int dimension       = 3;
    unsigned int       condition_size    = number_of_nodes * dimension;
    const bool check_has_rot_dof = rCurrentCondition->GetGeometry()[0].HasDofFor(ROTATION_X);

    if (check_has_rot_dof) condition_size *= 2;
    

    if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      unsigned int index = i * dimension;
      if (check_has_rot_dof) index *= 2;

      rValues[index]     = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[0];
      rValues[index + 1] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[1];
      rValues[index + 2] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[2];

        if (check_has_rot_dof)
        {          
          rValues[index+dimension]     = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_ANGULAR_VELOCITY )[0];
          rValues[index+dimension+1] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_ANGULAR_VELOCITY )[1];
          rValues[index+dimension+2] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_ANGULAR_VELOCITY )[2];          
        }
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
  
    (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

    (rCurrentCondition) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);
    if (rCurrentCondition->GetGeometry()[0].HasDofFor(ROTATION_X))
    {
    (rCurrentCondition) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, MOMENT_RESIDUAL, rCurrentProcessInfo);
    }
    KRATOS_CATCH( "" )
  }

  //***************************************************************************
  //***************************************************************************
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

             
    bool                mSchemeIsInitialized;


    TimeVariables       mTime;
    DeltaTimeParameters mDeltaTime;   


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

