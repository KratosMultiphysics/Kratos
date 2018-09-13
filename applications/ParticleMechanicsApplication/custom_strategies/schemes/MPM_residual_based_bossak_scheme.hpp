//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//

#if !defined(KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME )
#define  KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "containers/array_1d.h"
#include "solving_strategies/schemes/scheme.h"

#include "particle_mechanics_application.h"
#include "custom_utilities/mpm_boundary_rotation_utility.h"


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
class MPMResidualBasedBossakScheme: public Scheme<TSparseSpace,TDenseSpace>
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

        // System constants
        double c0;
        double c1;
        double c2;
        double c3;
        double c4;
        double c5;
        double c6;

        // Static-dynamic parameter
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
    KRATOS_CLASS_POINTER_DEFINITION( MPMResidualBasedBossakScheme );

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
    MPMResidualBasedBossakScheme(ModelPart& grid_model_part, unsigned int DomainSize, unsigned int BlockSize, double rAlpham=0,double rDynamic=1)
        :Scheme<TSparseSpace,TDenseSpace>(), mr_grid_model_part(grid_model_part), mRotationTool(DomainSize,BlockSize,IS_STRUCTURE)
    {
        // For pure Newmark Scheme
        mAlpha.f= 0;
        mAlpha.m= rAlpham;

        mNewmark.beta= 0.25;//(1.0+mAlpha.f-mAlpha.m)*(1.0+mAlpha.f-mAlpha.m)*0.25;
        mNewmark.gamma= 0.5;//0.5+mAlpha.f-mAlpha.m;

        mNewmark.static_dynamic= rDynamic;

        mDomainSize = DomainSize;
        mBlockSize  = BlockSize;

        // Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();

        mMatrix.M.resize(NumThreads);
        mMatrix.D.resize(NumThreads);

        mVector.v.resize(NumThreads);
        mVector.a.resize(NumThreads);
        mVector.ap.resize(NumThreads);
    }

    /** Copy Constructor.
     */
     MPMResidualBasedBossakScheme(MPMResidualBasedBossakScheme& rOther)
        :BaseType(rOther)
        ,mAlpha(rOther.mAlpha)
        ,mNewmark(rOther.mNewmark)
        ,mMatrix(rOther.mMatrix)
        ,mVector(rOther.mVector)
        ,mr_grid_model_part(rOther.mr_grid_model_part)
        ,mDomainSize(rOther.mDomainSize)
        ,mRotationTool(rOther.mDomainSize,rOther.mBlockSize,IS_STRUCTURE)
    {
    }

    /** Destructor.
     */
    virtual ~MPMResidualBasedBossakScheme
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
        return BaseTypePointer( new MPMResidualBasedBossakScheme(*this) );
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

        // Rotate displacement to the local coordinate system since Dx is in the local coordinate system
        // Do not confuse with the name RotateVelocities, what the function really do is to RotateDisplacements
        mRotationTool.RotateVelocities(r_model_part);

        // Update of displacement (by DOF)
        for (typename DofsArrayType::iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
        {
            if (i_dof->IsFree() )
            {
                i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
            }
        }

        // Rotate the displacement back to the original coordinate system to calculate the velocity and acceleration
        // Do not confuse with the name RecoverVelocities, what the function really do is to RecoverDisplacements
        mRotationTool.RecoverVelocities(r_model_part);
        
		#pragma omp parallel for
		for(int iter = 0; iter < static_cast<int>(r_model_part.Nodes().size()); ++iter)
		{
			auto i = r_model_part.NodesBegin() + iter;
			array_1d<double, 3 > & DeltaDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT);


            array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY, 0);
            const array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
            const array_1d<double, 3 > & PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

            UpdateVelocity(CurrentVelocity, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);
            UpdateAcceleration(CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);
		}

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    // Predicts the solution for the current step:
    void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
		#pragma omp parallel for
		for(int iter = 0; iter < static_cast<int>(r_model_part.Nodes().size()); ++iter)
		{
			
			auto i = r_model_part.NodesBegin() + iter;
			const array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);
            const array_1d<double, 3 > & PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3 > & CurrentDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);
            //array_1d<double, 3 > & ImposedDisplacement  = (i)->FastGetSolutionStepValue(IMPOSED_DISPLACEMENT);
            const array_1d<double, 3 > & PreviousAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 1);
            
            if ((i->pGetDof(DISPLACEMENT_X))->IsFixed() == false)
            {
                CurrentDisplacement[0] = 0.0;
            }
            else
            {
                CurrentDisplacement[0]  = PreviousDisplacement[0];// + ImposedDisplacement[0];
            }

            if (i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false)
            {
                CurrentDisplacement[1] = 0.0;
            }
            else
            {
                CurrentDisplacement[1]  = PreviousDisplacement[1];// + ImposedDisplacement[1];
            }

            if (i->HasDofFor(DISPLACEMENT_Z))
            {
                if (i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false)
                {
                    CurrentDisplacement[2] = 0.0;
                }
                else
                {
                    CurrentDisplacement[2]  = PreviousDisplacement[2];// + ImposedDisplacement[2];
                }
            }

            if (i->HasDofFor(PRESSURE))
            {
                const double& PreviousPressure    = (i)->FastGetSolutionStepValue(PRESSURE, 1);
                double& CurrentPressure     = (i)->FastGetSolutionStepValue(PRESSURE);

                if ((i->pGetDof(PRESSURE))->IsFixed() == false)
                    CurrentPressure = PreviousPressure;
            }

            // Updating time derivatives
            array_1d<double, 3 > & CurrentVelocity       = (i)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & CurrentAcceleration   = (i)->FastGetSolutionStepValue(ACCELERATION);

            UpdateVelocity(CurrentVelocity, CurrentDisplacement, PreviousVelocity, PreviousAcceleration);
            UpdateAcceleration (CurrentAcceleration, CurrentDisplacement, PreviousVelocity, PreviousAcceleration);
			
		}
    }

    //***************************************************************************
    //***************************************************************************

    /**
    This is the place to initialize the elements.
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

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    /**
    This is the place to initialize the conditions.
    This is intended to be called just once when the strategy is initialized
    */
    void InitializeConditions(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(this->mElementsAreInitialized==false) << "Before initilizing Conditions, initialize Elements FIRST" << std::endl;

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
                itCond->Initialize(); // Function to initialize the condition
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

        // LOOP OVER THE GRID NODES PERFORMED FOR CLEAR ALL NODAL INFORMATION
		#pragma omp parallel for
		for(int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
		{
			
			auto i = mr_grid_model_part.NodesBegin() + iter;
			if( (i)->SolutionStepsDataHas(NODAL_MOMENTUM) && (i)->SolutionStepsDataHas(NODAL_MASS) && (i)->SolutionStepsDataHas(NODAL_INERTIA))//&& (i)->SolutionStepsDataHas(NODAL_INTERNAL_FORCE) )
            {
                array_1d<double, 3 > & NodalMomentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                array_1d<double, 3 > & NodalInertia = (i)->FastGetSolutionStepValue(NODAL_INERTIA);
                double & NodalMass = (i)->FastGetSolutionStepValue(NODAL_MASS);
                double & NodalPressure = (i)->FastGetSolutionStepValue(PRESSURE,1);

                double & NodalDensity = (i)->FastGetSolutionStepValue(DENSITY);
                double & NodalAuxR = (i)->FastGetSolutionStepValue(AUX_R);
                array_1d<double, 3 > & NodalAuxRVel = (i)->FastGetSolutionStepValue(AUX_R_VEL);
                array_1d<double, 3 > & NodalAuxRAcc = (i)->FastGetSolutionStepValue(AUX_R_ACC);

                NodalMomentum.clear();
                NodalInertia.clear();
                NodalMass= 0.0;
                NodalPressure = 0.0;

                NodalDensity = 0.0;
                NodalAuxR = 0.0;
                NodalAuxRVel.clear();
                NodalAuxRAcc.clear();
            }

            if((i)->SolutionStepsDataHas(DISPLACEMENT) && (i)->SolutionStepsDataHas(VELOCITY) && (i)->SolutionStepsDataHas(ACCELERATION) )
            {
                array_1d<double, 3 > & NodalDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT);
                double & NodalPressure = (i)->FastGetSolutionStepValue(PRESSURE);
                array_1d<double, 3 > & NodalVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);
                array_1d<double, 3 > & NodalAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);
                array_1d<double, 3 > & DeltaNodalVelocity = (i)->FastGetSolutionStepValue(AUX_VELOCITY);
                array_1d<double, 3 > & DeltaNodalAcceleration = (i)->FastGetSolutionStepValue(AUX_ACCELERATION);

                NodalDisplacement.clear();
                NodalPressure = 0.0;
                NodalVelocity.clear();
                NodalAcceleration.clear();
                DeltaNodalVelocity.clear();
                DeltaNodalAcceleration.clear();
            }
		}

        double NormVel = 1.0;
        double NormAcc = 1.0;
        double NormPres = 1.0;
        double NormDeltaVel = 1.0;
        double NormDeltaAcc = 1.0;
        double NormDeltaPres = 1.0;

        int ItNum = 1;

        while (ItNum <2)
        {
            #pragma omp parallel for
			for( int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
			{
                auto i = mr_grid_model_part.NodesBegin() + iter;
                if( (i)->SolutionStepsDataHas(NODAL_MOMENTUM) && (i)->SolutionStepsDataHas(NODAL_MASS) && (i)->SolutionStepsDataHas(NODAL_INERTIA))//&& (i)->SolutionStepsDataHas(NODAL_INTERNAL_FORCE) )
                {

                    array_1d<double, 3 > & NodalMomentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                    array_1d<double, 3 > & NodalInertia = (i)->FastGetSolutionStepValue(NODAL_INERTIA);
                    array_1d<double, 3 > & DeltaNodalVelocity = (i)->FastGetSolutionStepValue(AUX_VELOCITY,1);
                    array_1d<double, 3 > & DeltaNodalAcceleration = (i)->FastGetSolutionStepValue(AUX_ACCELERATION,1);

                    double & NodalMass = (i)->FastGetSolutionStepValue(NODAL_MASS);
                    NodalMomentum.clear();
                    NodalInertia.clear();
                    DeltaNodalVelocity.clear();
                    DeltaNodalAcceleration.clear();

                    NodalMass = 0.0;

                    if(i->SolutionStepsDataHas(NODAL_MPRESSURE)) {
                        double & NodalMPressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);
                        NodalMPressure = 0.0;
                    }
                }
			}
                     
            NormVel = 0.0;
            NormAcc = 0.0;
            NormPres = 0.0;
            NormDeltaVel = 0.0;
            NormDeltaAcc = 0.0;
            NormDeltaPres = 0.0;

            int nodes_counter = 0;
            Scheme<TSparseSpace,TDenseSpace>::InitializeSolutionStep(r_model_part,A,Dx,b);

			#pragma omp parallel for reduction(+:NormVel,NormAcc,NormPres,NormDeltaVel,NormDeltaAcc,NormDeltaPres)
			for(int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
			{
			
			    auto i = mr_grid_model_part.NodesBegin() + iter;
			    double & NodalMass     = (i)->FastGetSolutionStepValue(NODAL_MASS);
                
                double DeltaNodalPressure = 0.0;

                if (NodalMass > 1.0e-16 )
                {
                    array_1d<double, 3 > & DeltaNodalVelocity = (i)->FastGetSolutionStepValue(AUX_VELOCITY,1);
                    array_1d<double, 3 > & DeltaNodalAcceleration = (i)->FastGetSolutionStepValue(AUX_ACCELERATION,1);

                    array_1d<double, 3 > & NodalMomentum     = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                    array_1d<double, 3 > & NodalInertia    = (i)->FastGetSolutionStepValue(NODAL_INERTIA);

                    array_1d<double, 3 > & NodalVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);
                    array_1d<double, 3 > & NodalAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);
                    double & NodalPressure = (i)->FastGetSolutionStepValue(PRESSURE,1);
                    
                    if (i->HasDofFor(PRESSURE) && i->SolutionStepsDataHas(NODAL_MPRESSURE))
                    {
                        double & NodalMPressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);
                        DeltaNodalPressure = NodalMPressure/NodalMass;
                    }
            
                    DeltaNodalVelocity = NodalMomentum/NodalMass;
                    DeltaNodalAcceleration = NodalInertia/NodalMass;

                    NodalVelocity += DeltaNodalVelocity;
                    NodalAcceleration += DeltaNodalAcceleration;

                    NodalPressure += DeltaNodalPressure;

                    NormDeltaVel += (DeltaNodalVelocity[0]*DeltaNodalVelocity[0]+DeltaNodalVelocity[1]*DeltaNodalVelocity[1]+DeltaNodalVelocity[2]*DeltaNodalVelocity[2]);
                    NormDeltaAcc += (DeltaNodalAcceleration[0]*DeltaNodalAcceleration[0]+DeltaNodalAcceleration[1]*DeltaNodalAcceleration[1]+DeltaNodalAcceleration[2]*DeltaNodalAcceleration[2]);
                    NormDeltaPres += (DeltaNodalPressure * DeltaNodalPressure);

                    NormVel += (NodalVelocity[0]*NodalVelocity[0]+NodalVelocity[1]*NodalVelocity[1]+NodalVelocity[2]*NodalVelocity[2]);
                    NormAcc += (NodalAcceleration[0]*NodalAcceleration[0]+NodalAcceleration[1]*NodalAcceleration[1]+NodalAcceleration[2]*NodalAcceleration[2]);
                    NormPres += (NodalPressure * NodalPressure);
                    
                    ++nodes_counter;
                }
            }

            NormVel = sqrt(NormVel);
            NormAcc = sqrt(NormAcc);
            NormPres = sqrt(NormPres);

            NormDeltaVel = sqrt(NormDeltaVel);
            NormDeltaAcc = sqrt(NormDeltaAcc);
            NormDeltaPres = sqrt(NormDeltaPres);
            
            ++ItNum;
        }

        double DeltaTime = CurrentProcessInfo[DELTA_TIME];
        KRATOS_ERROR_IF(DeltaTime == 0) << "Detected delta_time = 0 in the Solution Scheme ... check if the time step is created correctly for the current model part" << std::endl;

        //initializing Newmark constants
        mNewmark.c0 = ( 1.0 / (mNewmark.beta * DeltaTime * DeltaTime) );
        mNewmark.c1 = ( mNewmark.gamma / (mNewmark.beta * DeltaTime) );
        mNewmark.c2 = ( 1.0 / (mNewmark.beta * DeltaTime) );
        mNewmark.c3 = ( 0.5 / (mNewmark.beta) - 1.0 );
        mNewmark.c4 = ( (mNewmark.gamma / mNewmark.beta) - 1.0  );
        mNewmark.c5 = ( DeltaTime * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2 ) );

        KRATOS_CATCH( "" )
    }


    //***************************************************************************
    //***************************************************************************
    /**
    Function called once at the end of a solution step, after convergence is reached if
    an iterative process is needed
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

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
        for (ModelPart::NodeIterator i = mr_grid_model_part.NodesBegin();
                i != mr_grid_model_part.NodesEnd(); ++i)
        {
            if( (i)->SolutionStepsDataHas(EXTERNAL_FORCE) && (i)->SolutionStepsDataHas(INTERNAL_FORCE) )
            {
                array_1d<double, 3 > & ExternalForce = (i)->FastGetSolutionStepValue(EXTERNAL_FORCE);
                array_1d<double, 3 > & InternalForce = (i)->FastGetSolutionStepValue(INTERNAL_FORCE);
                ExternalForce.clear();
                InternalForce.clear();
            }
        }
    }

    //***************************************************************************
    //***************************************************************************

    /** This function is designed to be called in the builder and solver to introduce*/

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

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

        // If there is a slip condition, apply it on a rotated system of coordinates
        mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());
        mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());
       
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

        // If there is a slip condition, apply it on a rotated system of coordinates
        mRotationTool.RotateRHS(RHS_Contribution,rCurrentElement->GetGeometry());
        mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentElement->GetGeometry());
        
        KRATOS_CATCH( "" )

    }

    //***************************************************************************
    //***************************************************************************

    /** Functions totally analogous to the precedent but applied to
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

        // Basic operations for the element considered
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

        // Rotate contributions (to match coordinates for slip conditions)
        mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());
        mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());        
        

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

        // Basic operations for the element considered
        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {
            (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

            (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread], CurrentProcessInfo);
        }

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {
            AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);
        }

        // Rotate contributions (to match coordinates for slip conditions)
        mRotationTool.Rotate(RHS_Contribution,rCurrentCondition->GetGeometry());
        mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentCondition->GetGeometry());        
        
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

        //check that the variables are correctly initialized
        KRATOS_ERROR_IF(DISPLACEMENT.Key() == 0) <<"DISPLACEMENT has Key zero! (check if the application is correctly registered"<<std::endl;
        KRATOS_ERROR_IF(VELOCITY.Key() == 0) <<"VELOCITY has Key zero! (check if the application is correctly registered"<<std::endl;
        KRATOS_ERROR_IF(ACCELERATION.Key() == 0) <<"ACCELERATION has Key zero! (check if the application is correctly registered"<<std::endl;
        
        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin();
                it!=r_model_part.NodesEnd(); it++)
        {
            KRATOS_ERROR_IF(it->SolutionStepsDataHas(DISPLACEMENT) == false) << "DISPLACEMENT variable is not allocated for node "<< it->Id() <<std::endl;
            KRATOS_ERROR_IF(it->SolutionStepsDataHas(VELOCITY) == false) << "VELOCITY variable is not allocated for node "<< it->Id() <<std::endl;
            KRATOS_ERROR_IF(it->SolutionStepsDataHas(ACCELERATION) == false) << "ACCELERATION variable is not allocated for node " << it->Id() <<std::endl;
        }

        //check that dofs exist
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin();
                it!=r_model_part.NodesEnd(); it++)
        {
            KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_X) == false) <<"Missing DISPLACEMENT_X dof on node "<<it->Id() <<std::endl;
            KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_Y) == false) <<"Missing DISPLACEMENT_Y dof on node "<<it->Id() <<std::endl;
            KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_Z) == false) <<"Missing DISPLACEMENT_Z dof on node "<<it->Id() <<std::endl;
        }

        //check for admissible value of the AlphaBossak
        KRATOS_ERROR_IF(mAlpha.m > 0.0 || mAlpha.m < -0.3) << "Value not admissible for AlphaBossak. Admissible values should be between 0.0 and -0.3. Current value is "<< mAlpha.m << std::endl;

        //check for minimum value of the buffer index
        KRATOS_ERROR_IF(r_model_part.GetBufferSize() < 2) << "Insufficient buffer size. Buffer size should be greater than 2. Current size is" << r_model_part.GetBufferSize() <<std::endl;

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

    ModelPart& mr_grid_model_part;

    unsigned int    mDomainSize;
    unsigned int    mBlockSize;

    MPMBoundaryRotationUtility<LocalSystemMatrixType,LocalSystemVectorType> mRotationTool;

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

        // Adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += M * (1-mAlpha.m) * mNewmark.c0 * mNewmark.static_dynamic;

        }

        // Adding  damping contribution
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

        // Adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= (1.00 - mAlpha.m) * mNewmark.static_dynamic ;

            rCurrentElement->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mAlpha.m * mVector.ap[thread] * mNewmark.static_dynamic;

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

        // Adding damping contribution
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

        // Adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentCondition->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= (1.00 - mAlpha.m) * mNewmark.static_dynamic;

            rCurrentCondition->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mAlpha.m * mVector.ap[thread] * mNewmark.static_dynamic;

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

        // Adding damping contribution
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
}; /* Class MPMResidualBasedBossakScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME defined */


