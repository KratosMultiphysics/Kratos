//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_BOSSAK_SCHEME )
#define  KRATOS_RESIDUAL_BASED_BOSSAK_SCHEME

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
#include "pfem_solid_mechanics_application_variables.h"

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
/*@} */

template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedBossakScheme: public Scheme<TSparseSpace,TDenseSpace>
{
protected:

    struct GeneralAlphaMethod
    {

        double f;  //alpha Hilbert
        double m;  //alpha Bosssak

    };

    struct NewmarkMethod
    {

        double beta;
        double gamma;

        //system constants
        double c0;
        double c1;
        double c2;
        double c3;
        double c4;
        double c5;
        double c6;

        //static-dynamic parameter
        double static_dynamic;

    };


    struct  GeneralMatrices
    {

        std::vector< Matrix > M;     //first derivative matrix  (usually mass matrix)
        std::vector< Matrix > D;     //second derivative matrix (usually damping matrix)

    };

    struct GeneralVectors
    {

        std::vector< Vector > v;    //velocity
        std::vector< Vector > a;    //acceleration
        std::vector< Vector > ap;   //previous acceleration

    };



public:


    /**@name Type Definitions */

    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedBossakScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

    typedef typename BaseType::TDataType                         TDataType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename Element::DofsVectorType                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;
    
    typedef typename BaseType::Pointer                     BaseTypePointer;
  
    /*@} */

    /**
     * Constructor.
     * The bossak method
     */
    ResidualBasedBossakScheme(double rAlpham=0,double rDynamic=1)
        :Scheme<TSparseSpace,TDenseSpace>()
    {
        //For pure Newmark Scheme
        mAlpha.f= 0;
        mAlpha.m= rAlpham;

        mNewmark.beta= (1.0+mAlpha.f-mAlpha.m)*(1.0+mAlpha.f-mAlpha.m)*0.25;
        mNewmark.gamma= 0.5+mAlpha.f-mAlpha.m;

        mNewmark.static_dynamic= rDynamic;

	//std::cout << " MECHANICAL SCHEME: The Bossak Time Integration Scheme [alpha_m= "<<mAlpha.m<<" beta= "<<mNewmark.beta<<" gamma= "<<mNewmark.gamma<<"]"<<std::endl;


        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();

        mMatrix.M.resize(NumThreads);
        mMatrix.D.resize(NumThreads);

        mVector.v.resize(NumThreads);
        mVector.a.resize(NumThreads);
        mVector.ap.resize(NumThreads);


    }


    /** Copy Constructor.
     */
    ResidualBasedBossakScheme(ResidualBasedBossakScheme& rOther)
      :BaseType(rOther)
      ,mAlpha(rOther.mAlpha)
      ,mNewmark(rOther.mNewmark)
      ,mMatrix(rOther.mMatrix)
      ,mVector(rOther.mVector)
    {
    }


    /** Destructor.
     */
    virtual ~ResidualBasedBossakScheme
    () {}

   /*@} */
    /**@name Operators
     */
    /*@{ */


    /**
     * Clone 
     */
    BaseTypePointer Clone() override
    {
      return BaseTypePointer( new ResidualBasedBossakScheme(*this) );
    }



    //***************************************************************************
    //***************************************************************************

    /**
     * Performing the update of the solution
     * Incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param r_model_part
     * @param rDofSet set of all primary variables
     * @param A	LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b ) override
    {
        KRATOS_TRY

        //std::cout << " Update " << std::endl;
        //update of displacement (by DOF)
        for (typename DofsArrayType::iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
        {
            if (i_dof->IsFree() )
            {
                i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
            }
        }

        //updating time derivatives (nodally for efficiency)
        array_1d<double, 3 > DeltaDisplacement;
        for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                i != r_model_part.NodesEnd(); ++i)
        {

            noalias(DeltaDisplacement) = (i)->FastGetSolutionStepValue(DISPLACEMENT) - (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);


            array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY, 0);
            array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
            array_1d<double, 3 > & PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

            UpdateVelocity     (CurrentVelocity, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

            UpdateAcceleration (CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

        }

        KRATOS_CATCH( "" )
    }


    //***************************************************************************
    //***************************************************************************

    //predicts the solution for the current step:
    // x = xold + vold * Dt

    void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        //std::cout << " Prediction " << std::endl;
        array_1d<double, 3 > DeltaDisplacement;

        //double DeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];


        for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                i != r_model_part.NodesEnd(); ++i)
        {
            //predicting displacement = PreviousDisplacement + PreviousVelocity * DeltaTime;
            //ATTENTION::: the prediction is performed only on free nodes

            array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);
            array_1d<double, 3 > & PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3 > & CurrentDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & ImposedDisplacement  = (i)->FastGetSolutionStepValue(IMPOSED_DISPLACEMENT);


            if ((i->pGetDof(DISPLACEMENT_X))->IsFixed() == false)
            {
                CurrentDisplacement[0] = PreviousDisplacement[0];// + DeltaTime * PreviousVelocity[0];
            }
            else
            {
                CurrentDisplacement[0]  = PreviousDisplacement[0] + ImposedDisplacement[0];//to impose fixed displacements;
                //PreviousDisplacement[0] = 0;
            }

            if (i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false)
            {
                CurrentDisplacement[1] = PreviousDisplacement[1]; //+ DeltaTime * PreviousVelocity[1];
            }
            else
            {
                CurrentDisplacement[1]  = PreviousDisplacement[1] + ImposedDisplacement[1];//to impose fixed displacements;
                //PreviousDisplacement[1] = 0;
            }


            if (i->HasDofFor(DISPLACEMENT_Z))
            {
                if (i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false)
                {
                    CurrentDisplacement[2] = PreviousDisplacement[2]; // + DeltaTime * PreviousVelocity[2];
                }
                else
                {
                    CurrentDisplacement[2]  = PreviousDisplacement[2] + ImposedDisplacement[2];//to impose fixed displacements;
                    //PreviousDisplacement[2] = 0;
                }
            }


            // std::cout<<" DispPre "<<PreviousDisplacement<<" ID "<<i->Id()<<std::endl;
            // std::cout<<" DispCur "<<CurrentDisplacement<<" ID "<<i->Id()<<std::endl;

            if (i->HasDofFor(PRESSURE))
            {
                double& PreviousPressure    = (i)->FastGetSolutionStepValue(PRESSURE, 1);
                double& CurrentPressure     = (i)->FastGetSolutionStepValue(PRESSURE);

                if ((i->pGetDof(PRESSURE))->IsFixed() == false)
                    CurrentPressure = PreviousPressure;

                //std::cout<<" PressureCur [1] "<<CurrentPressure<<" PressurePre [1] "<<PreviousPressure<<" ID "<<i->Id()<<std::endl;
            }



            //updating time derivatives ::: please note that displacements and its time derivatives
            //can not be consistently fixed separately

            noalias(DeltaDisplacement) = CurrentDisplacement - PreviousDisplacement;

            array_1d<double, 3 > & PreviousAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 1);
            array_1d<double, 3 > & CurrentVelocity       = (i)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & CurrentAcceleration   = (i)->FastGetSolutionStepValue(ACCELERATION);

            UpdateVelocity     (CurrentVelocity, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

            UpdateAcceleration (CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

        }

    }

    //***************************************************************************
    //***************************************************************************

    /**
    this is the place to initialize the elements.
    This is intended to be called just once when the strategy is initialized
     */
    void InitializeElements(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin() + ElementPartition[k];
            ElementsArrayType::iterator ElemEnd = rModelPart.Elements().begin() + ElementPartition[k + 1];

            for (ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)
            {
                itElem->Initialize(); //function to initialize the element
            }

        }

        this->mElementsAreInitialized = true;

        //std::cout<<" mechanical elements are initialized "<<std::endl;

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    /**
    this is the place to initialize the conditions.
    This is intended to be called just once when the strategy is initialized
    */
    void InitializeConditions(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        if(this->mElementsAreInitialized==false)
            KRATOS_THROW_ERROR( std::logic_error, "Before initilizing Conditions, initialize Elements FIRST", "" )

            int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Conditions().size(), NumThreads, ConditionPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin() + ConditionPartition[k];
            ConditionsArrayType::iterator CondEnd = rModelPart.Conditions().begin() + ConditionPartition[k + 1];

            for (ConditionsArrayType::iterator itCond = CondBegin; itCond != CondEnd; itCond++)
            {
                itCond->Initialize(); //function to initialize the condition
            }

        }

        this->mConditionsAreInitialized = true;
        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    /**
     * initializes time step solution
     * only for reasons if the time step solution is restarted
     * @param r_model_part
     * @param A	LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        ProcessInfo CurrentProcessInfo= r_model_part.GetProcessInfo();


        Scheme<TSparseSpace,TDenseSpace>::InitializeSolutionStep(r_model_part,A,Dx,b);


        double DeltaTime = CurrentProcessInfo[DELTA_TIME];

        //if (DeltaTime == 0)
	  //KRATOS_THROW_ERROR( std::logic_error, "detected delta_time = 0 in the Solution Scheme ... check if the time step is created correctly for the current model part", "" )

	if (DeltaTime == 0){
	  std::cout<<" WARNING: detected delta_time = 0 in the Solution Scheme "<<std::endl;
	  std::cout<<" DELTA_TIME set to 1 considering a Quasistatic step with one step only "<<std::endl;
	  std::cout<<" PLEASE : check if the time step is created correctly for the current model part "<<std::endl;
	  
	  CurrentProcessInfo[DELTA_TIME] = 1;
	  DeltaTime = CurrentProcessInfo[DELTA_TIME];
	}


        //initializing Newmark constants
        mNewmark.c0 = ( 1.0 / (mNewmark.beta * DeltaTime * DeltaTime) );
        mNewmark.c1 = ( mNewmark.gamma / (mNewmark.beta * DeltaTime) );
        mNewmark.c2 = ( 1.0 / (mNewmark.beta * DeltaTime) );
        mNewmark.c3 = ( 0.5 / (mNewmark.beta) - 1.0 );
        mNewmark.c4 = ( (mNewmark.gamma / mNewmark.beta) - 1.0  );
        mNewmark.c5 = ( DeltaTime * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2 ) );


        //std::cout<<" Newmark Variables "<<mNewmark.c0<<" "<<mNewmark.c1<<" "<<mNewmark.c2<<" "<<mNewmark.c3<<" "<<mNewmark.c4<<" "<<mNewmark.c5<<std::endl;

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    /**
    function called once at the end of a solution step, after convergence is reached if
    an iterative process is needed
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY
        //finalizes solution step for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ElementsArrayType::iterator ElementsBegin = rElements.begin() + ElementPartition[k];
            ElementsArrayType::iterator ElementsEnd = rElements.begin() + ElementPartition[k + 1];

            for (ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                itElem->FinalizeSolutionStep(CurrentProcessInfo);
            }
        }

        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ConditionsArrayType::iterator ConditionsBegin = rConditions.begin() + ConditionPartition[k];
            ConditionsArrayType::iterator ConditionsEnd = rConditions.begin() + ConditionPartition[k + 1];

            for (ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; itCond++)
            {
                itCond->FinalizeSolutionStep(CurrentProcessInfo);
            }
        }
        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void InitializeNonLinIteration(ModelPart& r_model_part,
                                   TSystemMatrixType& A,
                                   TSystemVectorType& Dx,
                                   TSystemVectorType& b) override
    {
        KRATOS_TRY
        ElementsArrayType& pElements = r_model_part.Elements();
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for (ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            (it) -> InitializeNonLinearIteration(CurrentProcessInfo);
        }

        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            (it) -> InitializeNonLinearIteration(CurrentProcessInfo);
        }
        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void InitializeNonLinearIteration(Condition::Pointer rCurrentCondition,
                                      ProcessInfo& CurrentProcessInfo) override
    {
        (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }


    //***************************************************************************
    //***************************************************************************

    void InitializeNonLinearIteration(Element::Pointer rCurrentElement,
                                      ProcessInfo& CurrentProcessInfo) override
    {
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }

    //***************************************************************************
    //***************************************************************************

    /** this function is designed to be called in the builder and solver to introduce*/

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread],CurrentProcessInfo);

            (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToLHS (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

            AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        }

        //AssembleTimeSpaceLHS(rCurrentElement, LHS_Contribution, DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {

        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

            (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);
        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);
        }

        KRATOS_CATCH( "" )

    }

    //***************************************************************************
    //***************************************************************************

    /** functions totally analogous to the precedent but applied to
          the "condition" objects
    */
    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {


        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

	
        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

            (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);
	
        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToLHS  (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

            AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);
        }

        //AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix,CurrentProcessInfo);


        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current condition
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

            (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread], CurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

        //adding the dynamic contributions (static is already included)

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        }

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    /** Function that returns the list of Degrees of freedom to be
    assembled in the system for a Given Element
     */
    void GetElementalDofList(
        Element::Pointer rCurrentElement,
        Element::DofsVectorType& ElementalDofList,
        ProcessInfo& CurrentProcessInfo) override
    {
        rCurrentElement->GetDofList(ElementalDofList, CurrentProcessInfo);
    }

    //***************************************************************************
    //***************************************************************************

    /** Function that returns the list of Degrees of freedom to be
    assembled in the system for a Given Element
     */
    void GetConditionDofList(
        Condition::Pointer rCurrentCondition,
        Element::DofsVectorType& ConditionDofList,
        ProcessInfo& CurrentProcessInfo) override
    {
        rCurrentCondition->GetDofList(ConditionDofList, CurrentProcessInfo);
    }

    //***************************************************************************
    //***************************************************************************

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    int Check(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        int err = Scheme<TSparseSpace, TDenseSpace>::Check(r_model_part);
        if(err!=0) return err;

        //check for variables keys
        //verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )

        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin();
                it!=r_model_part.NodesEnd(); it++)
        {
            if (it->SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if (it->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if (it->SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
        }

        //check that dofs exist
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin();
                it!=r_model_part.NodesEnd(); it++)
        {
            if(it->HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_X dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Y dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Z dof on node ",it->Id() )
        }


        //check for admissible value of the AlphaBossak
        if(mAlpha.m > 0.0 || mAlpha.m < -0.3)
            KRATOS_THROW_ERROR( std::logic_error,"Value not admissible for AlphaBossak. Admissible values should be between 0.0 and -0.3. Current value is ", mAlpha.m )

            //check for minimum value of the buffer index
            //verify buffer size
            if (r_model_part.GetBufferSize() < 2)
                KRATOS_THROW_ERROR( std::logic_error, "insufficient buffer size. Buffer size should be greater than 2. Current size is", r_model_part.GetBufferSize() )


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

protected:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */
    
    GeneralAlphaMethod  mAlpha;
    NewmarkMethod       mNewmark;

    GeneralMatrices     mMatrix;
    GeneralVectors      mVector;

    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

    inline void UpdateVelocity(array_1d<double, 3 > & CurrentVelocity,
                               const array_1d<double, 3 > & DeltaDisplacement,
                               const array_1d<double, 3 > & PreviousVelocity,
                               const array_1d<double, 3 > & PreviousAcceleration)
    {

        noalias(CurrentVelocity) =  (mNewmark.c1 * DeltaDisplacement - mNewmark.c4 * PreviousVelocity
                                     - mNewmark.c5 * PreviousAcceleration) * mNewmark.static_dynamic;

    }


    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************

    inline void UpdateAcceleration(array_1d<double, 3 > & CurrentAcceleration,
                                   const array_1d<double, 3 > & DeltaDisplacement,
                                   const array_1d<double, 3 > & PreviousVelocity,
                                   const array_1d<double, 3 > & PreviousAcceleration)
    {

        noalias(CurrentAcceleration) =  (mNewmark.c0 * DeltaDisplacement - mNewmark.c2 * PreviousVelocity
                                         -  mNewmark.c3 * PreviousAcceleration) * mNewmark.static_dynamic;


    }


    //Elements:
    //****************************************************************************

    /**
    Atangent = M*c0 + D*c1 + K

     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {

        // adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += M * (1-mAlpha.m) * mNewmark.c0 * mNewmark.static_dynamic;

            //std::cout<<" Mass Matrix "<<M<<" coeficient "<<(1-mAlpha.m)*mNewmark.c0<<std::endl;
        }

        //adding  damping contribution
        if (D.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += D * (1-mAlpha.f) * mNewmark.c1 * mNewmark.static_dynamic;

        }

    }

    //Elements:
    //****************************************************************************

    /**
    bdyn = b - M*a - D*v

     */
    void AddDynamicsToRHS(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= (1.00 - mAlpha.m) * mNewmark.static_dynamic ;

            rCurrentElement->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mAlpha.m * mVector.ap[thread] * mNewmark.static_dynamic;

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
            //KRATOS_WATCH( prod(M, macc[thread] ) )

        }

        //adding damping contribution
        if (D.size1() != 0)
        {
            rCurrentElement->GetFirstDerivativesVector(mVector.v[thread], 0);

            (mVector.v[thread]) *= mNewmark.static_dynamic ;

            noalias(RHS_Contribution) -= prod(D, mVector.v[thread]);
        }


    }


    //Conditions:
    //****************************************************************************

    void AddDynamicsToRHS(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentCondition->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= (1.00 - mAlpha.m) * mNewmark.static_dynamic;

            rCurrentCondition->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mAlpha.m * mVector.ap[thread] * mNewmark.static_dynamic;

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

        //adding damping contribution
        //damping contribution
        if (D.size1() != 0)
        {
            rCurrentCondition->GetFirstDerivativesVector(mVector.v[thread], 0);

            (mVector.v[thread]) *= mNewmark.static_dynamic ;

            noalias(RHS_Contribution) -= prod(D, mVector.v [thread]);
        }

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
    /**@name Unaccessible methods */
    /*@{ */
}; /* Class ResidualBasedBossakScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BOSSAK_SCHEME defined */


