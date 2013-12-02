//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_RIGID_WALL_CONTACT_CALCULATION_UTILITIES_H_INCLUDED )
#define  KRATOS_RIGID_WALL_CONTACT_CALCULATION_UTILITIES_H_INCLUDED


/* System includes */
#include <cmath>

/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/model_part.h"
#include "custom_conditions/wall_tip_condition.hpp"

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
  class RigidWallContactCalculationUtilities
  {

  private:

    enum ContactFace{ FreeSurface, RakeSurface, TipSurface, ClearanceSurface };

    typedef struct
    {
      double Radius;         //tool tip radius
      double RakeAngle;      //top angle, from vertical axis
      double ClearanceAngle; //bottom angle, from the horizontal axis
      
      double m_factor;  //tan(0.5*pi-RakeAngle) == cotangent(RakeAngle) == 1/tan(RakeAngle)
      double n_factor;  //tan(ClearanceAngle)

      array_1d<double,3>  OriginalCenter;  //tool tip arch center original position
      array_1d<double,3>  Center;    //tool tip arch center current position
      array_1d<double,3>  Velocity;  //tool velocity, expressed on the center velocity

    } RigidWallVariables;

  
    typedef struct
    {
      double gapN;              
      double penalty_factor;    //penalty to compute contact stiffness
      
      array_1d<double,3>  Point;      //point of the workpiece
      array_1d<double,3>  Projection; //point projection on tool
      array_1d<double,3>  Normal;     //contact normal on the tool

      bool active;

    } ContactVariables;


    typedef struct
    {

      double penalty_parameter; //factor to compute the penalty
      double young_modulus;     //factor to compute the penalty

    } MaterialVariables;

  public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( RigidWallContactCalculationUtilities );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

    typedef typename BaseType::TDataType                         TDataType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename Element::DofsVectorType                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
  
    typedef ModelPart::NodesContainerType                   NodesArrayType;

    typedef ModelPart::ElementsContainerType             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

    typedef Node<3>                                               NodeType;

    typedef Condition::EquationIdVectorType           EquationIdVectorType;       
    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    RigidWallContactCalculationUtilities(){mActiveWall=false;};
  
    /** Destructor.
     */
    ~RigidWallContactCalculationUtilities(){};

    /** Operators.
     */

    void SetEquationSystemSize(int EquationSystemSize)
    {
        mEquationSystemSize = EquationSystemSize;
    }

    //**************************************************************************
    //**************************************************************************
  
    void SetRigidWall (double Radius,
		       double RakeAngle,
		       double ClearanceAngle,
		       array_1d<double,3>  Center,
		       array_1d<double,3>  Velocity)
    {
      mActiveWall      = true;

      pi = 3.141592654;

      mRigidWall.Radius     =  Radius;
      mRigidWall.RakeAngle  =  RakeAngle * pi / 180;
      mRigidWall.ClearanceAngle = ClearanceAngle * pi / 180;

      mRigidWall.m_factor = tan(0.5*pi-mRigidWall.RakeAngle);
      mRigidWall.n_factor = tan(mRigidWall.ClearanceAngle);

      mRigidWall.OriginalCenter  = Center;
      mRigidWall.Center    = Center;
      mRigidWall.Velocity  = Velocity;

      std::cout<<" [TOOL:                        ] "<<std::endl;
      std::cout<<" [Radius:"<<mRigidWall.Radius<<"            ] "<<std::endl;
      std::cout<<" [Center:"<<mRigidWall.Center<<"  ] "<<std::endl;
      std::cout<<" [Velocity:"<<mRigidWall.Velocity<<"      ] "<<std::endl;
      std::cout<<" [Rake:"<<mRigidWall.RakeAngle<<"               ] "<<std::endl;
      std::cout<<" [Clearance:"<<mRigidWall.ClearanceAngle<<"          ] "<<std::endl;
      std::cout<<" [m_factor:"<<mRigidWall.m_factor<<"          ] "<<std::endl;
      std::cout<<" [n_factor:"<<mRigidWall.n_factor<<"          ] "<<std::endl;

    }


   //**************************************************************************
    //**************************************************************************
  
    array_1d<double,3> GetWallTipCenter ()
    {
      return mRigidWall.Center;
    }

   //**************************************************************************
    //**************************************************************************
  
    double GetWallTipRadius ()
    {
      return mRigidWall.Radius;
    }

    //**************************************************************************
    //**************************************************************************
  
    void SetContactParameters (double penalty_parameter,
			       double body_elastic_modulus)
    {
      mMaterial.penalty_parameter = penalty_parameter;
      mMaterial.young_modulus = body_elastic_modulus;
    }

    //**************************************************************************
    //**************************************************************************
    
    int Check(ModelPart& r_model_part,bool active)
    {
      KRATOS_TRY

	if(mActiveWall==active && mMaterial.young_modulus>0) 
	  return 0; //return ok
	else
	  return 1;

      KRATOS_CATCH( "" )
    }

    
    //**************************************************************************
    //**************************************************************************
    void Build(ModelPart& r_model_part,
	       TSystemMatrixType& A,
	       TSystemVectorType& b)
    
    { 
      KRATOS_TRY

      //getting the array of the conditions
      NodesArrayType& NodesArray = r_model_part.Nodes();

      //std::cout<<" Nodes "<<r_model_part.NumberOfNodes()<<std::endl;

      //contributions to the system
      LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
      LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

      //vector containing the localization in the system of the different
      //terms
      Condition::EquationIdVectorType EquationId;

      //double StartTime = GetTickCount();
	
      InitializeNonLinearIteration(r_model_part);

      // assemble all elements
#ifndef _OPENMP
      ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

      // assemble all conditions
      for ( NodesArrayType::ptr_iterator nd = NodesArray.ptr_begin(); nd != NodesArray.ptr_end(); ++nd)
	{
	  if(nd->Is(BOUNDARY)){
	  
	    //calculate elemental contribution
	    CalculateSystemContributions(*nd, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

	    //assemble the elemental contribution
	    AssembleLHS(A, LHS_Contribution, EquationId);
	    AssembleRHS(b, RHS_Contribution, EquationId);
	  }
	  
	}

#else
      //creating an array of lock variables of the size of the system matrix
      std::vector< omp_lock_t > lock_array(A.size1());

      int A_size = A.size1();
      for (int i = 0; i < A_size; i++)
	omp_init_lock(&lock_array[i]);

      //create a partition of the element array
      int number_of_threads = omp_get_max_threads();
      double start_prod     = omp_get_wtime();

      vector<unsigned int> nodes_partition;
      CreatePartition(number_of_threads, NodesArray.size(), nodes_partition);

      // if( r_model_part.GetCommunicator().MyPID() == 0)
      // 	{
      // 	  KRATOS_WATCH( number_of_threads )
      // 	  KRATOS_WATCH( nodes_partition )
      // 	}

#pragma omp parallel for
      for (int k = 0; k < number_of_threads; k++)
	{
	  //contributions to the system
	  LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
	  LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

	  Condition::EquationIdVectorType EquationId;

	  ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

	  typename NodesArrayType::ptr_iterator nd_begin = NodesArray.ptr_begin() + nodes_partition[k];
	  typename NodesArrayType::ptr_iterator nd_end   = NodesArray.ptr_begin() + nodes_partition[k + 1];

	  // assemble all elements
	  for (typename NodesArrayType::ptr_iterator nd = nd_begin; nd != nd_end; ++nd)
	    {
	      if((*nd)->Is(BOUNDARY)){

		//calculate elemental contribution
		CalculateSystemContributions(*nd, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
		
		//assemble the elemental contribution
		Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, lock_array);
		
		// #pragma omp critical
		// {
		//   //assemble the elemental contribution
		//   AssembleLHS(A,LHS_Contribution,EquationId);
		//   AssembleRHS(b,RHS_Contribution,EquationId);
		// }

	      }
	    }
	}



      double stop_prod = omp_get_wtime();
      if (r_model_part.GetCommunicator().MyPID() == 0)
	std::cout << "time: " << stop_prod - start_prod << std::endl;

      for (int i = 0; i < A_size; i++)
	omp_destroy_lock(&lock_array[i]);

      // if(r_model_part.GetCommunicator().MyPID() == 0)
      // 	{
      // 	  KRATOS_WATCH( "finished parallel building" )
      // 	}

      //                        //ensure that all the threads are syncronized here
      //                        #pragma omp barrier
#endif


      KRATOS_CATCH( "" )

	}


    //**************************************************************************
    //**************************************************************************
    void BuildRHS(ModelPart& r_model_part,
		  TSystemVectorType& b)
    
    { 
      KRATOS_TRY

      //getting the array of the conditions
      NodesArrayType& NodesArray = r_model_part.Nodes();

      //std::cout<<" Nodes "<<r_model_part.NumberOfNodes()<<std::endl;

      //contributions to the system
      LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

      //vector containing the localization in the system of the different
      //terms
      Condition::EquationIdVectorType EquationId;

      //double StartTime = GetTickCount();
	
      InitializeNonLinearIteration(r_model_part);

      // assemble all elements
#ifndef _OPENMP
      ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

      // assemble all conditions
      for (typename NodesArrayType::ptr_iterator nd = NodesArray.ptr_begin(); nd != NodesArray.ptr_end(); ++nd)
	{
	  if(nd->Is(BOUNDARY)){
	  
	    //calculate elemental contribution
	    CalculateSystemRHS(*nd, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

	    //assemble the elemental contribution
	    AssembleRHS(b, RHS_Contribution, EquationId);
	  }
	  
	}

#else
      //creating an array of lock variables of the size of the system matrix
      std::vector< omp_lock_t > lock_array(b.size());

      int b_size = b.size();
      for (int i = 0; i < b_size; i++)
	omp_init_lock(&lock_array[i]);

      //create a partition of the element array
      int number_of_threads = omp_get_max_threads();
      double start_prod     = omp_get_wtime();

      vector<unsigned int> nodes_partition;
      CreatePartition(number_of_threads, NodesArray.size(), nodes_partition);

      // if( r_model_part.GetCommunicator().MyPID() == 0)
      // 	{
      // 	  KRATOS_WATCH( number_of_threads )
      // 	  KRATOS_WATCH( nodes_partition )
      // 	}

#pragma omp parallel for
      for (int k = 0; k < number_of_threads; k++)
	{
	  //contributions to the system
	  LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
	  LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

	  Condition::EquationIdVectorType EquationId;

	  ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

	  typename NodesArrayType::ptr_iterator nd_begin = NodesArray.ptr_begin() + nodes_partition[k];
	  typename NodesArrayType::ptr_iterator nd_end   = NodesArray.ptr_begin() + nodes_partition[k + 1];

	  // assemble all elements
	  for (typename NodesArrayType::ptr_iterator nd = nd_begin; nd != nd_end; ++nd)
	    {
	      if((*nd)->Is(BOUNDARY)){

		//calculate elemental contribution
		CalculateSystemRHS(*nd, RHS_Contribution, EquationId, CurrentProcessInfo);
		
		//assemble the elemental contribution
		//AssembleRHS(b, RHS_Contribution, EquationId, lock_array);
		
                #pragma omp critical
		{
		  //assemble the elemental contribution
		//   AssembleLHS(A,LHS_Contribution,EquationId);
		  AssembleRHS(b,RHS_Contribution,EquationId);
		}

	      }
	    }
	}



      double stop_prod = omp_get_wtime();
      if (r_model_part.GetCommunicator().MyPID() == 0)
	std::cout << "time: " << stop_prod - start_prod << std::endl;

      for (int i = 0; i < b_size; i++)
	omp_destroy_lock(&lock_array[i]);

      // if(r_model_part.GetCommunicator().MyPID() == 0)
      // 	{
      // 	  KRATOS_WATCH( "finished parallel building" )
      // 	}

      //                        //ensure that all the threads are syncronized here
      //                        #pragma omp barrier
#endif


      KRATOS_CATCH( "" )

	}


    //**************************************************************************
    //**************************************************************************

    void InitializeSolutionStep( ModelPart& r_model_part,
				 TSystemMatrixType& A,
				 TSystemVectorType& Dx,
				 TSystemVectorType& b )
    {
      KRATOS_TRY

      ProcessInfo& CurrentProcessInfo= r_model_part.GetProcessInfo();	  
      double Time = CurrentProcessInfo[TIME];

      mRigidWall.Center = mRigidWall.OriginalCenter +  mRigidWall.Velocity * Time;

      if (Time == 0)
	KRATOS_ERROR( std::logic_error, "detected delta_time = 0 in the Solution Scheme ... check if the time step is created correctly for the current model part", "" )

      NodesArrayType& NodesArray = r_model_part.Nodes();
      for (typename NodesArrayType::ptr_iterator nd = NodesArray.ptr_begin(); nd != NodesArray.ptr_end(); ++nd)
	{
	  if((*nd)->FastGetSolutionStepValue(RIGID_WALL)==true){
    
	    (*nd)->Set(RIGID);
	    //(*nd)->Set(STRUCTURE);

	  }

	}


      KRATOS_CATCH( "" )
	}

    //**************************************************************************
    //**************************************************************************

    void FinalizeSolutionStep( ModelPart& r_model_part,
			       TSystemMatrixType& A,
			       TSystemVectorType& Dx,
			       TSystemVectorType& b )
    {
      KRATOS_TRY
	
      //getting the array of the conditions
      NodesArrayType& NodesArray = r_model_part.Nodes();
      ProcessInfo& CurrentProcessInfo= r_model_part.GetProcessInfo();	  
      double DeltaTime = CurrentProcessInfo[DELTA_TIME];
      
      for (typename NodesArrayType::ptr_iterator nd = NodesArray.ptr_begin(); nd != NodesArray.ptr_end(); ++nd)
	{
	  if((*nd)->FastGetSolutionStepValue(RIGID_WALL)==true){
	  
	    (*nd)->FastGetSolutionStepValue(DISPLACEMENT) += mRigidWall.Velocity * DeltaTime;
	    
	  }

	}
      
      //calculate elemental contribution
      KRATOS_CATCH( "" )      
     }

    //**************************************************************************
    //**************************************************************************

    void InitializeNonLinearIteration( ModelPart& r_model_part)
    {
      ProcessInfo& CurrentProcessInfo= r_model_part.GetProcessInfo();

      //getting the array of the conditions
      NodesArrayType& NodesArray = r_model_part.Nodes();

      for ( typename NodesArrayType::ptr_iterator nd = NodesArray.ptr_begin(); nd != NodesArray.ptr_end(); ++nd)
	{
	  ClearNodalForces(*nd);
	}

      CurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] = 0;
      CurrentProcessInfo[NUMBER_OF_STICK_CONTACTS]  = 0;
      CurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS]   = 0;

      std::cout<<" Set ACTIVE CONTACTS to zero "<<CurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS]<<std::endl;
    }

    //**************************************************************************
    //**************************************************************************
  
    void FinalizeNonLinearIteration( ModelPart& r_model_part)
    {

    }



    //************************************************************************************
    //************************************************************************************

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


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    //**************************************************************************

    void AssembleLHS(TSystemMatrixType& A,
		     LocalSystemMatrixType& LHS_Contribution,
		     Element::EquationIdVectorType& EquationId)
    {
        unsigned int local_size = LHS_Contribution.size1();
	
        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];
            if (i_global < mEquationSystemSize)
            {
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < mEquationSystemSize)
                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                }
            }
        }
    }
    


    //**************************************************************************

    void AssembleRHS(TSystemVectorType& b,
		     LocalSystemVectorType& RHS_Contribution,
		     Element::EquationIdVectorType& EquationId)
    {
      unsigned int local_size = RHS_Contribution.size();
      
      for (unsigned int i_local = 0; i_local < local_size; i_local++)
	{
	  unsigned int i_global = EquationId[i_local];
	  if (i_global < mEquationSystemSize) //on "free" DOFs
	    {
	      // ASSEMBLING THE SYSTEM VECTOR
	      b[i_global] += RHS_Contribution[i_local];
	    }
	  }
      
    }


    //************************************************************************************
    //************************************************************************************

#ifdef _OPENMP

    void Assemble(
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        std::vector< omp_lock_t >& lock_array
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

            if (i_global < mEquationSystemSize)
            {
                omp_set_lock(&lock_array[i_global]);

                b[i_global] += RHS_Contribution(i_local);
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < mEquationSystemSize)
                    {
                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                    }
                }

                omp_unset_lock(&lock_array[i_global]);


            }
            //note that computation of reactions is not performed here!
        }
    }
    
#endif


    //************************************************************************************
    //************************************************************************************

    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++)
            partitions[i] = partitions[i - 1] + partition_size;
    }


    //************************************************************************************
    //************************************************************************************

    void ClearNodalForces(NodeType::Pointer rCurrentNode)
    {

      KRATOS_TRY

	array_1d<double, 3 > & ContactForceNormal  = rCurrentNode->FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);
      ContactForceNormal.clear();

      array_1d<double, 3 > & ContactForceTangent = rCurrentNode->FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);
      ContactForceTangent.clear();
      
	

      KRATOS_CATCH( "" )
	}
  
    //************************************************************************************
    //************************************************************************************

    void EquationIdVector(NodeType::Pointer rCurrentNode, EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
      unsigned int number_of_nodes = 1;
      unsigned int dimension =2; 
      unsigned int size = number_of_nodes*dimension; //dimension
      rResult.resize(size);
      unsigned int index;
      for (unsigned int i=0; i<number_of_nodes; i++)
	{
	  index = i*dimension;
	  rResult[index]   = (rCurrentNode->GetDof(DISPLACEMENT_X).EquationId());
	  rResult[index+1] = (rCurrentNode->GetDof(DISPLACEMENT_Y).EquationId());
	
	}
    }
  
    //************************************************************************************
    //************************************************************************************
  
    void GetDofList(NodeType::Pointer rCurrentNode,DofsVectorType& NodeDofList,ProcessInfo& CurrentProcessInfo)
    {
      unsigned int number_of_nodes = 1;
      unsigned int dimension =2; 
      unsigned int size = number_of_nodes*dimension; //dimension
      NodeDofList.resize(size);
      unsigned int index;
      for (unsigned int i=0; i<number_of_nodes; i++)
	{
	  index = i*dimension;
	  NodeDofList[index]   = (rCurrentNode->pGetDof(DISPLACEMENT_X));
	  NodeDofList[index+1] = (rCurrentNode->pGetDof(DISPLACEMENT_Y));	
	}
    }


    //************************************************************************************
    //************************************************************************************
  
    void CalculateLocalRHS(NodeType::Pointer rCurrentNode, LocalSystemVectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      if(rRightHandSideVector.size() != 2)
	rRightHandSideVector.resize(2,false);


      noalias(rRightHandSideVector) = ZeroVector(2);
      
      ContactVariables Contact;
      Contact.active = false;

      switch( ContactSearch(rCurrentNode,Contact) )
	{
	  
	case FreeSurface:      //no contact in the current node
	  break;
	case RakeSurface:      CalculateRakeLocalRHS(rCurrentNode,rRightHandSideVector,rCurrentProcessInfo,Contact);
	  break;
	case TipSurface:       CalculateTipLocalRHS(rCurrentNode,rRightHandSideVector,rCurrentProcessInfo,Contact);
	  break;
	case ClearanceSurface: CalculateClearanceLocalRHS(rCurrentNode,rRightHandSideVector,rCurrentProcessInfo,Contact);
	  break;
	default:               //no contact in the current node
	  break;
	}


      if(Contact.active){
	rCurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] += 1;
      }

      // if(Contact.active)
      // 	std::cout<<" Contact.active "<<Contact.active<<" Fcont "<<rCurrentNode->FastGetSolutionStepValue(FORCE_CONTACT_NORMAL)<<" Center Position "<<mRigidWall.Center<<std::endl;


      KRATOS_CATCH( "" )
	}


    //************************************************************************************
    //************************************************************************************

    void CalculateSystemRHS ( NodeType::Pointer rCurrentNode,
			      LocalSystemVectorType& RHS_Contribution,
			      Element::EquationIdVectorType& EquationId,
			      ProcessInfo& CurrentProcessInfo)
    {
      KRATOS_TRY

	//int thread = OpenMPUtils::ThisThread();

      CalculateLocalRHS (rCurrentNode,RHS_Contribution,CurrentProcessInfo);
      EquationIdVector  (rCurrentNode,EquationId,CurrentProcessInfo);

      KRATOS_CATCH( "" )
	}
  


    //************************************************************************************
    //************************************************************************************
  
    void CalculateLocalSystem(NodeType::Pointer rCurrentNode, LocalSystemMatrixType& rLeftHandSideMatrix, LocalSystemVectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      if(rLeftHandSideMatrix.size1() != 2)
        rLeftHandSideMatrix.resize(2,2,false);

      if(rRightHandSideVector.size() != 2)
	rRightHandSideVector.resize(2,false);


      noalias(rLeftHandSideMatrix)  = ZeroMatrix(2,2);
      noalias(rRightHandSideVector) = ZeroVector(2);
      
      ContactVariables Contact;
      Contact.active = false;

      switch( ContactSearch(rCurrentNode,Contact) )
	{
	  
	case FreeSurface:      //no contact in the current node
	  break;
	case RakeSurface:      CalculateRakeLocalSystem(rCurrentNode,rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo,Contact);
	  break;
	case TipSurface:       CalculateTipLocalSystem(rCurrentNode,rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo,Contact);
	  break;
	case ClearanceSurface: CalculateClearanceLocalSystem(rCurrentNode,rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo,Contact);
	  break;
	default:               //no contact in the current node
	  break;
	}


      if(Contact.active){
	rCurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] += 1;
      }

      // if(Contact.active)
      // 	std::cout<<" Contact.active "<<Contact.active<<" Fcont "<<rCurrentNode->FastGetSolutionStepValue(FORCE_CONTACT_NORMAL)<<" Center Position "<<mRigidWall.Center<<std::endl;


      KRATOS_CATCH( "" )
	}


    //************************************************************************************
    //************************************************************************************

    void CalculateSystemContributions( NodeType::Pointer rCurrentNode,
				       LocalSystemMatrixType& LHS_Contribution,
				       LocalSystemVectorType& RHS_Contribution,
				       Element::EquationIdVectorType& EquationId,
				       ProcessInfo& CurrentProcessInfo)
    {
      KRATOS_TRY

	//int thread = OpenMPUtils::ThisThread();

      CalculateLocalSystem (rCurrentNode,LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
      EquationIdVector     (rCurrentNode,EquationId,CurrentProcessInfo);

      KRATOS_CATCH( "" )
	}
  


    //************************************************************************************
    //************************************************************************************


    void CalculateRakeLocalRHS(NodeType::Pointer rCurrentNode, LocalSystemVectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,ContactVariables &rContact)
    {
      KRATOS_TRY
      
      //1.-compute contact normal
      rContact.Normal[0] = -cos(mRigidWall.RakeAngle);
      rContact.Normal[1] =  sin(mRigidWall.RakeAngle);
      rContact.Normal[2] = 0;

      //2.-compute point projection
      
      rContact.Projection[0] = rContact.Point[0] + mRigidWall.m_factor * ( rContact.Point[1] - mRigidWall.m_factor * ( mRigidWall.Radius * cos(mRigidWall.RakeAngle) - mRigidWall.Center[0] ) - mRigidWall.Center[1] ) - mRigidWall.Radius * sin(mRigidWall.RakeAngle);

      rContact.Projection[1] = mRigidWall.m_factor * ( mRigidWall.m_factor * rContact.Point[1] + rContact.Point[0] + mRigidWall.Radius * cos(mRigidWall.RakeAngle) - mRigidWall.Center[0] ) + mRigidWall.Center[1] + mRigidWall.Radius * sin(mRigidWall.RakeAngle);
      
      rContact.Projection[2] = 0;
      
      rContact.Projection /= (1+mRigidWall.m_factor * mRigidWall.m_factor);


      if(mRigidWall.RakeAngle == 0){
	rContact.Projection[0] = mRigidWall.Center[0] - mRigidWall.Radius; 
	rContact.Projection[1] = rContact.Point[1]; 
      }

      //3.-compute gap

      rContact.gapN = inner_prod((rContact.Point - rContact.Projection),rContact.Normal);
      

      // if(rContact.gapN<0)
      // 	std::cout<<" RakeSystem [ gapN: "<<rContact.gapN<<" ] "<<std::endl;

      //4.-set active contact (if gap > 0 means overlapping -> contact active)
      
      if(rContact.gapN<0){
	
	rContact.active = true;
 
      
	//5.- compute penalty_factor
      
	rContact.penalty_factor =  CalculatePenaltyFactor (rCurrentNode,rContact);
      
      
	//6.-compute contact force

	array_1d<double,3>& force = rCurrentNode->GetSolutionStepValue(FORCE_CONTACT_NORMAL);
      
	force = (-1) * (rContact.penalty_factor *  rContact.gapN) * rContact.Normal;

	rRightHandSideVector[0] = force[0];
	rRightHandSideVector[1] = force[1];


      }
      else{
	
	rContact.active = false;
      }


      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    void CalculateTipLocalRHS(NodeType::Pointer rCurrentNode, LocalSystemVectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,ContactVariables &rContact)
    {
      KRATOS_TRY

 
      //1.-compute point projection
      
	rContact.Projection = mRigidWall.Radius * ( (rContact.Point-mRigidWall.Center)/ norm_2(rContact.Point-mRigidWall.Center) ) + mRigidWall.Center;

      
      //2.-compute contact normal
      rContact.Normal = (rContact.Projection-mRigidWall.Center)/mRigidWall.Radius;


      //3.-compute gap

      if( norm_2(mRigidWall.Center-rContact.Point) <= mRigidWall.Radius ){
	rContact.gapN = (-1) * norm_2(rContact.Point - rContact.Projection);
      }
      else{
	rContact.gapN = norm_2(rContact.Projection - rContact.Point);
      }
      

      // if(rContact.gapN<0)
      // 	std::cout<<" TipSystem [ gapN: "<<rContact.gapN<<" ] "<<std::endl;
      

      //4.-set active contact (if gap > 0 means overlapping -> contact active)
      
      if(rContact.gapN<0){
	
	rContact.active = true;
 
 
	//5.- compute penalty_factor
      
	rContact.penalty_factor =  CalculatePenaltyFactor (rCurrentNode,rContact);
      
      
	//6.-compute contact force

	array_1d<double,3>& force = rCurrentNode->GetSolutionStepValue(FORCE_CONTACT_NORMAL);
      
	force = (-1) * (rContact.penalty_factor *  rContact.gapN) * rContact.Normal;

	rRightHandSideVector[0] = force[0];
	rRightHandSideVector[1] = force[1];

      }
      else{
	
	rContact.active = false;
      }


    
      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    void CalculateClearanceLocalRHS(NodeType::Pointer rCurrentNode, LocalSystemVectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,ContactVariables &rContact)
    {
      KRATOS_TRY

      //1.-compute contact normal
      rContact.Normal[0] =  sin(mRigidWall.ClearanceAngle);
      rContact.Normal[1] = -cos(mRigidWall.ClearanceAngle);
      rContact.Normal[2] = 0;

      //2.-compute point projection
      
      rContact.Projection[0] = rContact.Point[0] + mRigidWall.n_factor * ( rContact.Point[1] + mRigidWall.n_factor * ( mRigidWall.Center[0] + mRigidWall.Radius * sin(mRigidWall.ClearanceAngle) ) - mRigidWall.Center[1] + mRigidWall.Radius * cos(mRigidWall.ClearanceAngle) );

      rContact.Projection[1] = mRigidWall.n_factor * ( mRigidWall.n_factor * rContact.Point[1] + rContact.Point[0] + mRigidWall.Center[0] - mRigidWall.Radius * sin(mRigidWall.ClearanceAngle) ) + mRigidWall.Center[1] - mRigidWall.Radius * cos(mRigidWall.ClearanceAngle);
      
      rContact.Projection[2] = 0;
      
      rContact.Projection /= (1+mRigidWall.n_factor * mRigidWall.n_factor);


      //3.-compute gap

      rContact.gapN = inner_prod((rContact.Point - rContact.Projection), rContact.Normal);
      
      // if(rContact.gapN<0)
      // 	std::cout<<" ClearanceSystem [ gapN: "<<rContact.gapN<<" ] "<<std::endl;

      //4.-set active contact (if gap > 0 means overlapping -> contact active)
      
      if(rContact.gapN<0){
	
	rContact.active = true;
 
      
	//5.- compute penalty_factor
      
	rContact.penalty_factor =  CalculatePenaltyFactor (rCurrentNode,rContact); 
      
	//6.-compute contact force

	array_1d<double,3>& force = rCurrentNode->GetSolutionStepValue(FORCE_CONTACT_NORMAL);
      
	force = (-1) * (rContact.penalty_factor *  rContact.gapN) * rContact.Normal;

	rRightHandSideVector[0] = force[0];
	rRightHandSideVector[1] = force[1];

      }
      else{
	
	rContact.active = false;
      }


    
      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    void CalculateRakeLocalSystem(NodeType::Pointer rCurrentNode, LocalSystemMatrixType& rLeftHandSideMatrix, LocalSystemVectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,ContactVariables &rContact)
    {
      KRATOS_TRY
      
      //1.-compute contact normal
      rContact.Normal[0] = -cos(mRigidWall.RakeAngle);
      rContact.Normal[1] =  sin(mRigidWall.RakeAngle);
      rContact.Normal[2] = 0;

      //2.-compute point projection
      
      rContact.Projection[0] = rContact.Point[0] + mRigidWall.m_factor * ( rContact.Point[1] - mRigidWall.m_factor * ( mRigidWall.Radius * cos(mRigidWall.RakeAngle) - mRigidWall.Center[0] ) - mRigidWall.Center[1] ) - mRigidWall.Radius * sin(mRigidWall.RakeAngle);

      rContact.Projection[1] = mRigidWall.m_factor * ( mRigidWall.m_factor * rContact.Point[1] + rContact.Point[0] + mRigidWall.Radius * cos(mRigidWall.RakeAngle) - mRigidWall.Center[0] ) + mRigidWall.Center[1] + mRigidWall.Radius * sin(mRigidWall.RakeAngle);
      
      rContact.Projection[2] = 0;
      
      rContact.Projection /= (1+mRigidWall.m_factor * mRigidWall.m_factor);


      if(mRigidWall.RakeAngle == 0){
	rContact.Projection[0] = mRigidWall.Center[0] - mRigidWall.Radius; 
	rContact.Projection[1] = rContact.Point[1]; 
      }

      //3.-compute gap

      rContact.gapN = inner_prod((rContact.Point - rContact.Projection),rContact.Normal);
      

      // if(rContact.gapN<0)
      // 	std::cout<<" RakeSystem [ gapN: "<<rContact.gapN<<" ] "<<std::endl;

      //4.-set active contact (if gap > 0 means overlapping -> contact active)
      
      if(rContact.gapN<0){
	
	rContact.active = true;
 
      
	//5.- compute penalty_factor
      
	rContact.penalty_factor =  CalculatePenaltyFactor (rCurrentNode,rContact);
      
      
	//6.-compute contact force

	array_1d<double,3>& force = rCurrentNode->GetSolutionStepValue(FORCE_CONTACT_NORMAL);
      
	force = (-1) * (rContact.penalty_factor *  rContact.gapN) * rContact.Normal;

	rRightHandSideVector[0] = force[0];
	rRightHandSideVector[1] = force[1];


	//7.-compute contact stiffness

	noalias(rLeftHandSideMatrix) =  rContact.penalty_factor * outer_prod_2(rContact.Normal,rContact.Normal);

      }
      else{
	
	rContact.active = false;
      }


      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    void CalculateTipLocalSystem(NodeType::Pointer rCurrentNode, LocalSystemMatrixType& rLeftHandSideMatrix, LocalSystemVectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,ContactVariables &rContact)
    {
      KRATOS_TRY

 
      //1.-compute point projection
      
	rContact.Projection = mRigidWall.Radius * ( (rContact.Point-mRigidWall.Center)/ norm_2(rContact.Point-mRigidWall.Center) ) + mRigidWall.Center;

      
      //2.-compute contact normal
      rContact.Normal = (rContact.Projection-mRigidWall.Center)/mRigidWall.Radius;


      //3.-compute gap

      if( norm_2(mRigidWall.Center-rContact.Point) <= mRigidWall.Radius ){
	rContact.gapN = (-1) * norm_2(rContact.Point - rContact.Projection);
      }
      else{
	rContact.gapN = norm_2(rContact.Projection - rContact.Point);
      }
      

      // if(rContact.gapN<0)
      // 	std::cout<<" TipSystem [ gapN: "<<rContact.gapN<<" ] "<<std::endl;
      

      //4.-set active contact (if gap > 0 means overlapping -> contact active)
      
      if(rContact.gapN<0){
	
	rContact.active = true;
 
 
	//5.- compute penalty_factor
      
	rContact.penalty_factor =  CalculatePenaltyFactor (rCurrentNode,rContact);
      
      
	//6.-compute contact force

	array_1d<double,3>& force = rCurrentNode->GetSolutionStepValue(FORCE_CONTACT_NORMAL);
      
	force = (-1) * (rContact.penalty_factor *  rContact.gapN) * rContact.Normal;

	rRightHandSideVector[0] = force[0];
	rRightHandSideVector[1] = force[1];

	//7.-compute contact stiffness

	noalias(rLeftHandSideMatrix) = rContact.penalty_factor * outer_prod_2(rContact.Normal,rContact.Normal);

      }
      else{
	
	rContact.active = false;
      }


    
      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    void CalculateClearanceLocalSystem(NodeType::Pointer rCurrentNode, LocalSystemMatrixType& rLeftHandSideMatrix, LocalSystemVectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,ContactVariables &rContact)
    {
      KRATOS_TRY

      //1.-compute contact normal
      rContact.Normal[0] =  sin(mRigidWall.ClearanceAngle);
      rContact.Normal[1] = -cos(mRigidWall.ClearanceAngle);
      rContact.Normal[2] = 0;

      //2.-compute point projection
      
      rContact.Projection[0] = rContact.Point[0] + mRigidWall.n_factor * ( rContact.Point[1] + mRigidWall.n_factor * ( mRigidWall.Center[0] + mRigidWall.Radius * sin(mRigidWall.ClearanceAngle) ) - mRigidWall.Center[1] + mRigidWall.Radius * cos(mRigidWall.ClearanceAngle) );

      rContact.Projection[1] = mRigidWall.n_factor * ( mRigidWall.n_factor * rContact.Point[1] + rContact.Point[0] + mRigidWall.Center[0] - mRigidWall.Radius * sin(mRigidWall.ClearanceAngle) ) + mRigidWall.Center[1] - mRigidWall.Radius * cos(mRigidWall.ClearanceAngle);
      
      rContact.Projection[2] = 0;
      
      rContact.Projection /= (1+mRigidWall.n_factor * mRigidWall.n_factor);


      //3.-compute gap

      rContact.gapN = inner_prod((rContact.Point - rContact.Projection), rContact.Normal);
      
      // if(rContact.gapN<0)
      // 	std::cout<<" ClearanceSystem [ gapN: "<<rContact.gapN<<" ] "<<std::endl;

      //4.-set active contact (if gap > 0 means overlapping -> contact active)
      
      if(rContact.gapN<0){
	
	rContact.active = true;
 
      
	//5.- compute penalty_factor
      
	rContact.penalty_factor =  CalculatePenaltyFactor (rCurrentNode,rContact); 
      
	//6.-compute contact force

	array_1d<double,3>& force = rCurrentNode->GetSolutionStepValue(FORCE_CONTACT_NORMAL);
      
	force = (-1) * (rContact.penalty_factor *  rContact.gapN) * rContact.Normal;

	rRightHandSideVector[0] = force[0];
	rRightHandSideVector[1] = force[1];


	//7.-compute contact stiffness

	noalias(rLeftHandSideMatrix) =  rContact.penalty_factor * outer_prod_2(rContact.Normal,rContact.Normal);

      }
      else{
	
	rContact.active = false;
      }


    
      KRATOS_CATCH( "" )
	}


    //**************************************************************************
    //**************************************************************************

    ContactFace ContactSearch(NodeType::Pointer rCurrentNode,ContactVariables &rContact)
    {

      KRATOS_TRY

      ContactFace Face = FreeSurface;
      
      // array_1d<double, 3 > & CurrentDisplacement  = rCurrentNode->FastGetSolutionStepValue(DISPLACEMENT);
      // array_1d<double, 3 > & PreviousDisplacement = rCurrentNode->FastGetSolutionStepValue(DISPLACEMENT,1);
      array_1d<double, 3 > & ReferencePosition    = rCurrentNode->Coordinates();
	    
      rContact.Point = ReferencePosition; // + (CurrentDisplacement-PreviousDisplacement);

      // rContact.Point[0] = rCurrentNode->X();
      // rContact.Point[1] = rCurrentNode->Y();
      // rContact.Point[2] = rCurrentNode->Z();
      
      double FaceR = CalculateRakeFace(FaceR,rContact);
      double FaceT = CalculateTipFace(FaceT,rContact);
      double FaceC = CalculateClearanceFace(FaceC,rContact);      
	
      double Face1=0,Face2=0,Face3=0;
      CalculateAuxiliarFaces(Face1,Face2,Face3,rContact);

      rCurrentNode->Reset(WallTipCondition::WALL_TIP);

      if(mRigidWall.RakeAngle>0){
	if(FaceR<=0 && Face3>=0 && Face1>=0){
	  Face = RakeSurface;
	}
	else if(FaceT<=0 && Face3<0 && Face2<=0){
	  Face = TipSurface;
	  rCurrentNode->Set(WallTipCondition::WALL_TIP);
	}
	else if(FaceC>=0 && Face2>=0 && Face1<0){
	  Face = ClearanceSurface;
	}
	else{
	  Face = FreeSurface;
	}
      }
      else if(mRigidWall.RakeAngle==0){

	if(FaceR<=0 && Face3>=0 && Face1>=0){
	  Face = RakeSurface;
	}
	else if(FaceT<=0 && Face3<=0 && Face2<=0){
	  Face = TipSurface;
	  rCurrentNode->Set(WallTipCondition::WALL_TIP);
	}
	else if(FaceC>=0 && Face2>=0 && Face1<=0){
	  Face = ClearanceSurface;
	}
	else{
	  Face = FreeSurface;
	}

      }
      // if(FaceR<=0 && Face3>=0 && Face1>=0){
      // 	Face = RakeSurface;
      // }
      // else if(FaceT<=0 && Face3<=0 && Face2>=0){
      // 	Face = TipSurface;
      // }
      // else if(FaceC>=0 && Face2<=0 && Face1<=0){
      // 	Face = ClearanceSurface;
      // }
      // else{
      // 	Face = FreeSurface;
      // }
      

      // if(Face!=FreeSurface)
      //  	std::cout<<"Node "<<rCurrentNode->Id()<<" Position :"<<rContact.Point<<" FACE "<<Face<<" [ Fr:"<<FaceR<<", Ft:"<<FaceT<<", Fc:"<<FaceC<<", F1:"<<Face1<<", F2:"<<Face2<<", F3:"<<Face3<<" ]"<<std::endl;

      return Face;


      KRATOS_CATCH( "" )

	}

    //************************************************************************************
    //************************************************************************************


    double& CalculateRakeFace(double& Face,ContactVariables &rContact)
    {
      KRATOS_TRY
	    
      Face = rContact.Point[1] - mRigidWall.m_factor * ( rContact.Point[0] + mRigidWall.Radius * cos(mRigidWall.RakeAngle) - mRigidWall.Center[0]) - mRigidWall.Center[1] - mRigidWall.Radius * sin(mRigidWall.RakeAngle); 
 
      if(mRigidWall.RakeAngle == 0)
	Face = rContact.Point[0] - mRigidWall.Center[0] - mRigidWall.Radius; 


      return Face;
    
      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    double& CalculateTipFace(double& Face,ContactVariables &rContact)
    {
      KRATOS_TRY
      
      Face = pow((rContact.Point[0] - mRigidWall.Center[0]),2) + pow((rContact.Point[1] - mRigidWall.Center[1]),2) - mRigidWall.Radius * mRigidWall.Radius; 
      
      return Face;	
    
      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    double& CalculateClearanceFace(double& Face,ContactVariables &rContact)
    {
      KRATOS_TRY
      
      Face = rContact.Point[1] - mRigidWall.n_factor * ( rContact.Point[0] - mRigidWall.Center[0] - mRigidWall.Radius * sin(mRigidWall.ClearanceAngle)) - mRigidWall.Center[1]  +  mRigidWall.Radius * cos(mRigidWall.ClearanceAngle); 
 
      return Face;
 

    
      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    void CalculateAuxiliarFaces(double & Face1, double& Face2,  double& Face3, ContactVariables &rContact  )
    {
      KRATOS_TRY
	     
      // double newx =   mRigidWall.m_factor * ( mRigidWall.Center[0] - mRigidWall.Radius * cos(mRigidWall.RakeAngle) ) - mRigidWall.Center[1] - mRigidWall.Radius * sin(mRigidWall.RakeAngle);
      // newx += mRigidWall.Center[1] - mRigidWall.Radius * cos(mRigidWall.ClearanceAngle) - mRigidWall.n_factor * ( mRigidWall.Center[0]  + mRigidWall.Radius * sin(mRigidWall.ClearanceAngle) );
      // newx /= (mRigidWall.m_factor - mRigidWall.n_factor);

      // double newy =    mRigidWall.n_factor * ( mRigidWall.m_factor * ( mRigidWall.Center[0] - mRigidWall.Radius * cos(mRigidWall.RakeAngle) ) -  mRigidWall.Center[1] - mRigidWall.Radius * sin(mRigidWall.RakeAngle) );
      // newy += mRigidWall.m_factor * ( mRigidWall.Center[1] - mRigidWall.n_factor * ( mRigidWall.Center[0] + mRigidWall.Radius * sin(mRigidWall.ClearanceAngle) ) - mRigidWall.Radius * cos(mRigidWall.ClearanceAngle) );
      // newy /= (mRigidWall.m_factor - mRigidWall.n_factor);

      // double s = (newy - mRigidWall.Center[1])/(newx - mRigidWall.Center[0]);


      // Face1 = rContact.Point[1] - s * ( rContact.Point[0] - mRigidWall.Center[0] ) -  mRigidWall.Center[1];
 
	Face1 = rContact.Point[1] - tan(pi*0.25) * ( rContact.Point[0] - mRigidWall.Center[0] ) - mRigidWall.Center[1];

      // Face2 = rContact.Point[1] + mRigidWall.n_factor * ( rContact.Point[0] - mRigidWall.Center[0] ) - mRigidWall.Center[1];
	Face2 = rContact.Point[0]  + mRigidWall.n_factor * ( rContact.Point[1] - mRigidWall.Center[1] ) - mRigidWall.Center[0];

	Face3 = rContact.Point[1] + tan(mRigidWall.RakeAngle) * ( rContact.Point[0] - mRigidWall.Center[0] ) - mRigidWall.Center[1];

      if(mRigidWall.RakeAngle==0)
      	Face3 = rContact.Point[1]-mRigidWall.Center[1];


      // Face2 = mRigidWall.Center[0] - rContact.Point[0] + mRigidWall.n_factor * ( rContact.Point[1] - mRigidWall.Center[1] );

      // Face3 = mRigidWall.Center[0] - rContact.Point[0] + mRigidWall.m_factor * ( rContact.Point[1] - mRigidWall.Center[1] );

      // if(mRigidWall.RakeAngle==0)
      // 	Face3 = rContact.Point[1]-mRigidWall.Center[1];
     
      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    double CalculatePenaltyFactor(NodeType::Pointer rCurrentNode,ContactVariables &rContact)
    {

      KRATOS_TRY

      WeakPointerVector<Node<3> >& rN = rCurrentNode->GetValue(NEIGHBOUR_NODES);
      array_1d<double,3> Neighb_Point;
      double distance = 0;
      double counter = 0;
      for(unsigned int i = 0; i < rN.size(); i++)
	{
	  if(rN[i].Is(BOUNDARY)){
	    
	    Neighb_Point[0] = rN[i].X();
	    Neighb_Point[1] = rN[i].Y();
	    Neighb_Point[2] = rN[i].Z();
	    
	    distance += norm_2(rContact.Point-Neighb_Point);

	    counter ++;
	  }
	}

      distance /= counter;

      rContact.penalty_factor = distance * 20 * mMaterial.penalty_parameter * mMaterial.young_modulus;
      
      return rContact.penalty_factor;
    
      KRATOS_CATCH( "" )
    }



    static inline LocalSystemMatrixType outer_prod_2(const array_1d<double, 3>& a, const array_1d<double, 3>& b)
    {
        LocalSystemMatrixType A(2,2);
        A(0,0)=a[0]*b[0];
        A(0,1)=a[0]*b[1];
        A(1,0)=a[1]*b[0];
        A(1,1)=a[1]*b[1];
	//A(0,2)=a[0]*b[2];
	//A(1,2)=a[1]*b[2];
        //A(2,0)=a[2]*b[0];
        //A(2,1)=a[2]*b[1];
        //A(2,2)=a[2]*b[2];
        return A;
    }


   static inline double inner_prod(const array_1d<double, 3>& a, const array_1d<double, 3>& b)
    {
        double temp =a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
        return temp;
    }

    static inline double norm_2(const array_1d<double, 3>& a)
    {
        double temp = pow(a[0],2) + pow(a[1],2) + pow(a[2],2);
        temp = sqrt(temp);
        return temp;
    }

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
    
    double pi;

      /*@} */
      /**@name Member Variables */
      /*@{ */

    RigidWallVariables mRigidWall;
    MaterialVariables   mMaterial;

    bool         mActiveWall;
    unsigned int mEquationSystemSize;

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

  }; /* Class ComputeLineSearch */

  /*@} */

  /**@name Type Definitions */
  /*@{ */


  /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RIGID_WALL_CONTACT_CALCULATION_UTILITIES_H_INCLUDED defined */

