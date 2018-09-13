//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:            December 2015 $
//   Revision:            $Revision:                  0.1 $
//
//

#if !defined(KRATOS_EXPLICIT_HAMILTON_BUILDER_AND_SOLVER )
#define  KRATOS_EXPLICIT_HAMILTON_BUILDER_AND_SOLVER


/* System includes */
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
//#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"
#include "utilities/beam_math_utilities.hpp"

#include "solid_mechanics_application_variables.h"

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
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ExplicitHamiltonBuilderAndSolver : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( ExplicitHamiltonBuilderAndSolver );

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::ElementsArrayType ElementsArrayType;

    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    typedef BeamMathUtils<double> BeamMathUtilsType;

    typedef Quaternion<double> QuaternionType;
    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ExplicitHamiltonBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
      //std::cout<<" EXPLICIT HAMILTON BUILDER AND SOLVER "<<std::endl;
    }

    /** Destructor.
     */
    virtual ~ExplicitHamiltonBuilderAndSolver()
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */


    //**************************************************************************
    //**************************************************************************

    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A)
    {
        KRATOS_TRY

         //Set Nodal Mass to zero
         NodesArrayType& pNodes             = r_model_part.Nodes();
         ElementsArrayType& pElements       = r_model_part.Elements();
         ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();

         #ifdef _OPENMP
                 int number_of_threads = omp_get_max_threads();
         #else
                 int number_of_threads = 1;
         #endif

         vector<unsigned int> node_partition;
         OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

         vector<unsigned int> element_partition;
         OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);


         #pragma omp parallel
         {

           #pragma omp for

	   for(int k=0; k<number_of_threads; k++)
	     {
	       typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	       typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	       for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
		 {
		   double& nodal_mass           =  i->FastGetSolutionStepValue(NODAL_MASS);
		   Matrix& nodal_inertia_dyadic =  i->FastGetSolutionStepValue(INERTIA_DYADIC);

		   nodal_mass = 0.0;
		   nodal_inertia_dyadic = ZeroMatrix(3,3);
		 }
	     }

         }

         //Calculate and assemble Mass Matrix on nodes

         unsigned int index = 0;

         #pragma omp parallel
         {
	   int k = OpenMPUtils::ThisThread();
	   typename ElementsArrayType::iterator ElemBegin = pElements.begin() + element_partition[k];
	   typename ElementsArrayType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];

	   for (typename ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)
             {
	       Matrix MassMatrix;

	       Element::GeometryType& geometry = itElem->GetGeometry(); //element nodes

	       (itElem)->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

	       const unsigned int dimension   = geometry.WorkingSpaceDimension();

	       index = 0;
	       for (unsigned int i = 0; i <geometry.size(); i++)
                 {
		   index = i * ( dimension * 2 );

		   double& mass                 = geometry(i)->FastGetSolutionStepValue(NODAL_MASS);
		   Matrix& nodal_inertia_dyadic = geometry(i)->FastGetSolutionStepValue(INERTIA_DYADIC);

		   geometry(i)->SetLock();

		   mass += MassMatrix(index,index);
		   for(unsigned int k=0; k<dimension; k++)
		     {
		       for(unsigned int l=0; l<dimension; l++)
			 {
			   nodal_inertia_dyadic(k,l) += MassMatrix(k+index+dimension, l+index+dimension);
			 }
		     }
		   geometry(i)->UnSetLock();
                 }
             }
         }



         #pragma omp parallel
         {

           #pragma omp for

	   for(int k=0; k<number_of_threads; k++)
	     {
	       typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	       typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	       for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
		 {
		   double& nodal_mass           =  i->FastGetSolutionStepValue(NODAL_MASS);
		   Matrix& nodal_inertia_dyadic =  i->FastGetSolutionStepValue(INERTIA_DYADIC);

		   //inertia_dyadic is multiplied by the mass and divided by the total mass
		   //this is done to not increase the total inertia of the section
		   //the end is an average inertia for the nodal section
		   nodal_inertia_dyadic /= nodal_mass;

		   //std::cout<<" Node ["<<i->Id()<<"] (mass:"<<nodal_mass<<", inertia:"<<nodal_inertia_dyadic<<") "<<std::endl;

		 }
	     }

         }

        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b)
    {
        KRATOS_TRY


        //Set Nodal Liapunov variables to zero

 	//getting nodes from the model
        NodesArrayType& pNodes           = r_model_part.Nodes();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

	vector<unsigned int> node_partition;
	OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

        #pragma omp parallel
	{

          #pragma omp for

	  for(int k=0; k<number_of_threads; k++)
	    {
	      typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	      typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	      for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
		{
		  (i)->FastGetSolutionStepValue(FORCE_RESIDUAL).clear();
		  (i)->FastGetSolutionStepValue(MOMENT_RESIDUAL).clear();
		}
	    }

        }

        // Compute condition contributions to RHS.
        CalculateAndAddConditionsRHS(pScheme, r_model_part);

        // Compute element contributions to RHS.
        CalculateAndAddElementsRHS(pScheme, r_model_part);

	KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        if (!pScheme)
            KRATOS_THROW_ERROR( std::runtime_error, "No scheme provided!", "" )

	NodesArrayType& pNodes           = r_model_part.Nodes();
        ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

	double DeltaTime  = rCurrentProcessInfo[DELTA_TIME];
	double Alpha      = rCurrentProcessInfo[ALPHA_TRAPEZOIDAL_RULE];

        #ifdef _OPENMP
	int number_of_threads = omp_get_max_threads();
        #else
	int number_of_threads = 1;
        #endif

        vector<unsigned int> node_partition;
	OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

	//initial build to initialize all
	this->Build(pScheme, r_model_part, A, b);

         #pragma omp parallel
         {

           #pragma omp for

	   for(int k=0; k<number_of_threads; k++)
	     {
	       typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	       typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	       for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
		 {
		   //is active by default
		   bool node_is_active = true;
		   if( (i)->IsDefined(ACTIVE) )
		     node_is_active = (i)->Is(ACTIVE);

		   if(node_is_active)
		     {

		       //Rotation
		       array_1d<double, 3 >& CurrentStepRotation  = (i)->FastGetSolutionStepValue(STEP_ROTATION);

		       //Moment and momentum: Rotation Part
		       array_1d<double,3>& moment_residual        = (i)->FastGetSolutionStepValue(MOMENT_RESIDUAL);
		       array_1d<double,3>& rotation_momentum      = (i)->FastGetSolutionStepValue(ROTATION_MOMENTUM);

		       //moment_residual is (Mext-Mint) then the sign is changed respect the reference formulation
		       array_1d<double,3> Y = DeltaTime * ( rotation_momentum + Alpha * DeltaTime * moment_residual );

		       // (0) Solving Lyapunov - type  equation:  WRONG it must be done at integration point!!

		       int iter = 0;
		       int max_iters = 15;

		       double residual  = 1;
		       double tolerance = 1e-12;

		       Vector Residual  = ZeroVector(3); //current iteration residual vector
		       Vector Rotation  = ZeroVector(3); //current step rotation vector
		       Vector ReferenceResidual = ZeroVector(3);

		       for(unsigned int j=0; j<3; j++)
			 {
			   Rotation[j]            = 1e-5; //CurrentStepRotation[j];
			   ReferenceResidual[j]   = Y[j];
			 }


		       //std::cout<<" Node ["<<i->Id()<<"] (moment:"<<moment_residual<<",momentum:"<<rotation_momentum<<") norm residual "<<norm_2(ReferenceResidual)<<std::endl;


		       // if( norm_2(ReferenceResidual) <= tolerance ){

		       // 	 Rotation = ZeroVector(3);
		       // 	 CurrentStepRotation.clear();

		       // 	 //std::cout<<" Node ["<<i->Id()<<"] (STEP_ROTATION:"<<Rotation<<") "<<iter<<" reference_residual "<<norm_2(ReferenceResidual)<<std::endl;

		       // }
		       // else {

			 Matrix Tangent        = ZeroMatrix(3,3);
			 Matrix InverseTangent = ZeroMatrix(3,3);
			 double determinant    = 0;

			 for(unsigned int j=0; j<3; j++)
			   {
			     CurrentStepRotation[j] = Rotation[j];
			   }


			 while ( residual > tolerance && iter < max_iters )
			   {
			     // (1) Build the residual :

			     // double BuildTimeStart = Timer::GetTime();
			     this->Build(*i.base(), pScheme, r_model_part, A, b);
			     // BuildTime = Timer::GetTime() - BuildTime;

			     // if (this->GetEchoLevel() > 1)
			     // {
			     //   std::cout<<"Time Building :"<<BuildTime<<std::endl;
			     // }

			     array_1d<double,3>& LyapunovResidual = i->FastGetSolutionStepValue(RESIDUAL_LYAPUNOV);
			     Matrix& LyapunovTangent = i->FastGetSolutionStepValue(TANGENT_LYAPUNOV);


			     //std::cout<<"(id:"<<i->Id()<<") LyapunovResidual "<<LyapunovResidual<<" LyapunovTangent "<<LyapunovTangent<<std::endl;

			     for(unsigned int j=0; j<3; j++)
			       {
				 Residual[j] = ReferenceResidual[j] - LyapunovResidual[j];
			       }

			     // (3) Solve the system:
			     MathUtils<double>::InvertMatrix( LyapunovTangent, InverseTangent, determinant );

			     Rotation += prod( InverseTangent, Residual );


			     // (4) Update Step Rotations:
			     if( (i->pGetDof(ROTATION_X))->IsFree() )
			       CurrentStepRotation[0]  = Rotation[0];
			     if( (i->pGetDof(ROTATION_Y))->IsFree() )
			       CurrentStepRotation[1]  = Rotation[1];
			     if( (i->pGetDof(ROTATION_Z))->IsFree() )
			       CurrentStepRotation[2]  = Rotation[2];


			     // (5) Update Rotations:
			     bool update_rotations = false;

			     if( update_rotations ){

			       array_1d<double, 3 >& PreviousRotation     = (i)->FastGetSolutionStepValue(ROTATION,1);
			       array_1d<double, 3 >& CurrentRotation      = (i)->FastGetSolutionStepValue(ROTATION);

			       Vector CurrentRotationVector = ZeroVector(3);
			       for( unsigned int j=0; j<3; j++)
				 {
				   CurrentRotationVector[j]  = PreviousRotation[j];    //previous iteration total rotation
				 }

			       QuaternionType ReferenceRotationQuaternion = QuaternionType::FromRotationVector( CurrentRotationVector );
			       QuaternionType StepRotationQuaternion      = QuaternionType::FromRotationVector( Rotation );
			       QuaternionType CurrentRotationQuaternion   = StepRotationQuaternion * ReferenceRotationQuaternion;

			       CurrentRotationQuaternion.ToRotationVector( CurrentRotationVector );



			       for( unsigned int j=0; j<3; j++)
			     	 {
			     	   CurrentRotation[j] = CurrentRotationVector[j];
			     	 }

			     }

			     //(6) Check Residual:
			     residual = norm_2( Residual );

			     if (this->GetEchoLevel() > 1){
			       std::cout<<" ("<<iter<<") : Rotation "<<CurrentStepRotation<<" Residual norm: "<<residual<<std::endl;
			     }

			     iter++;

			   }



			 if( iter >= max_iters ){
			     std::cout<<" Node ["<<i->Id()<<"] : LYAPUNOV ROTATION EQUATION NOT CONVERGED, iters:"<<iter<<", residual: "<<residual<<" STEP_ROTATION:"<<Rotation<<std::endl;
			 }
			 else{
			   if (this->GetEchoLevel() > 1){
			     std::cout<<" Node ["<<i->Id()<<"] (STEP_ROTATION:"<<Rotation<<") iters: "<<iter<<" residual "<<residual<<std::endl;
			   }
			 }

		       }
		     }
		   //}

	     }
	 }

	KRATOS_CATCH("")

    }

    //***************************************************************************
    //***************************************************************************

    void Build(
	       Node<3>::Pointer pNode,
	       typename TSchemeType::Pointer pScheme,
	       ModelPart& r_model_part,
	       TSystemMatrixType& A,
	       TSystemVectorType& b)
    {

        KRATOS_TRY

	ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();

        //Set Nodal Liapunov variables to zero
	array_1d<double,3>& LyapunovResidual = (pNode)->FastGetSolutionStepValue(RESIDUAL_LYAPUNOV);
	Matrix& LyapunovTangent = (pNode)->FastGetSolutionStepValue(TANGENT_LYAPUNOV);

	LyapunovResidual.clear();
	LyapunovTangent = ZeroMatrix(3,3);

	//IT HAS TO BE MODIFIED TO ONLY COMPUTE THE NODE NEIGHBOUR ELEMENTS AND CONDITIONS:

	WeakPointerVector<Element >& rE = pNode->GetValue(NEIGHBOUR_ELEMENTS);

	//std::cout<<" node ("<<(pNode)->Id()<<"): "<<rE.size()<<std::endl;

	//vector containing the localization in the system of the different terms
        Element::EquationIdVectorType EquationId;

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

	for(WeakPointerVector< Element >::iterator ie = rE.begin(); ie!=rE.end(); ++ie)
	  {
            //calculate elemental contribution
            pScheme->CalculateSystemContributions(Element::Pointer( *(ie.base()) ), LHS_Contribution, RHS_Contribution, EquationId, rCurrentProcessInfo);

            // clean local elemental memory
            pScheme->CleanMemory(Element::Pointer( *(ie.base()) ));
	  }



        KRATOS_CATCH( "" )

    }


    //***************************************************************************
    //***************************************************************************

    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& b)
    {

        KRATOS_TRY

        //Set Nodal Liapunov variables to zero

 	//getting nodes from the model
        NodesArrayType& pNodes           = r_model_part.Nodes();
	//getting elements from the model
	ElementsArrayType& pElements     = r_model_part.Elements();
        //getting conditions from the model
        ConditionsArrayType& pConditions = r_model_part.Conditions();


	ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

	vector<unsigned int> node_partition;
	OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

        #pragma omp parallel
	{

          #pragma omp for

	  for(int k=0; k<number_of_threads; k++)
	    {
	      typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	      typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	      for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
		{
		   array_1d<double,3>& LyapunovResidual = i->FastGetSolutionStepValue(RESIDUAL_LYAPUNOV);
		   Matrix& LyapunovTangent = i->FastGetSolutionStepValue(TANGENT_LYAPUNOV);

		   LyapunovResidual.clear();
		   LyapunovTangent = ZeroMatrix(3,3);
		}
	    }

        }

	//IT HAS TO BE MODIFIED TO ONLY COMPUTE THE NODE NEIGHBOUR ELEMENTS AND CONDITIONS:
	bool compute_everything = false;
	if( compute_everything ){

	  //vector containing the localization in the system of the different terms
	  Element::EquationIdVectorType EquationId;

	  //contributions to the system
	  LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
	  LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);


	  for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
	    {
	      //calculate elemental contribution
	      pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, rCurrentProcessInfo);

	      // clean local elemental memory
	      pScheme->CleanMemory(*it);
	    }

	  LHS_Contribution.resize(0, 0, false);
	  RHS_Contribution.resize(0, false);

	  // assemble all conditions
	  for (typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it)
	    {
	      //calculate condition contribution
	      pScheme->Condition_CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, rCurrentProcessInfo);

	      //clean local condition memory
	      pScheme->CleanMemory(*it);

	    }
	}


        KRATOS_CATCH( "" )

    }

    //***************************************************************************
    //***************************************************************************


    void CalculateAndAddConditionsRHS(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part )
    {

    KRATOS_TRY

    ProcessInfo& rCurrentProcessInfo      = r_model_part.GetProcessInfo();
    ConditionsArrayType& pConditions      = r_model_part.Conditions();

#ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
#else
    int number_of_threads = 1;
#endif

    vector<unsigned int> condition_partition;
    OpenMPUtils::CreatePartition(number_of_threads, pConditions.size(), condition_partition);


    #pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
    {
       typename ConditionsArrayType::ptr_iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
       typename ConditionsArrayType::ptr_iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];

       for (typename ConditionsArrayType::ptr_iterator it= it_begin; it!=it_end; ++it)
       {

           LocalSystemVectorType RHS_Condition_Contribution = LocalSystemVectorType(0);

           Element::EquationIdVectorType EquationId; //Dummy

	   //is active by default
	   bool condition_is_active = true;
	   if( (*it)->IsDefined(ACTIVE) )
	     condition_is_active = (*it)->Is(ACTIVE);

	   if(condition_is_active)
	     {
	       pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Condition_Contribution, EquationId, rCurrentProcessInfo);
	     }


       }
    }

    KRATOS_CATCH("")
    }


    //***************************************************************************
    //***************************************************************************


    void CalculateAndAddElementsRHS(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part )
    {

        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements        = r_model_part.Elements();

#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

        #pragma omp parallel for
        for(int k=0; k<number_of_threads; k++)
        {
            typename ElementsArrayType::ptr_iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (typename ElementsArrayType::ptr_iterator it= it_begin; it!=it_end; ++it)
            {

                LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
                Element::EquationIdVectorType EquationId; //Dummy

		//is active by default
		bool element_is_active = true;
		if( (*it)->IsDefined(ACTIVE) )
		  element_is_active = (*it)->Is(ACTIVE);

		if(element_is_active)
		  {
		    pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, rCurrentProcessInfo);
		  }

            }
        }

        KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************


    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
	KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
    }

    //**************************************************************************
    //**************************************************************************

    void ApplyPointLoads(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b)
    {
    }

    /**
    this function is intended to be called at the end of the solution step to clean up memory
    storage not needed
     */
    void Clear()
    {
        this->mDofSet = DofsArrayType();

        if (this->mpReactionsVector != NULL)
            TSparseSpace::Clear((this->mpReactionsVector));
        //      this->mReactionsVector = TSystemVectorType();

        this->mpLinearSystemSolver->Clear();

        if (this->GetEchoLevel() > 1)
        {
            std::cout << "ExplicitHamiltonBuilderAndSolver Clear Function called" << std::endl;
        }
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    virtual int Check(ModelPart& r_model_part)
    {
        KRATOS_TRY

        return 0;

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
    /*@} */
    /**@name Protected member Variables */
    /*@{ */
    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    //******************************************************************************************
    //******************************************************************************************


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

}; /* Class ExplicitHamiltonBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */
/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_HAMILTON_BUILDER_AND_SOLVER  defined */

