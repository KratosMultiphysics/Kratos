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
        int num_threads = OpenMPUtils::GetNumThreads();

        mMatrix.M.resize(num_threads);
        mMatrix.D.resize(num_threads);

        mVector.v.resize(num_threads);
        mVector.a.resize(num_threads);
        mVector.ap.resize(num_threads);
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
			const array_1d<double, 3 > & delta_displacement = (i)->FastGetSolutionStepValue(DISPLACEMENT);

            array_1d<double, 3 > & current_velocity     = (i)->FastGetSolutionStepValue(VELOCITY, 0);
            array_1d<double, 3 > & current_acceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 0);

            const array_1d<double, 3 > & previous_velocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);
            const array_1d<double, 3 > & previous_acceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);


            UpdateVelocity(current_velocity, delta_displacement, previous_velocity, previous_acceleration);
            UpdateAcceleration(current_acceleration, delta_displacement, previous_velocity, previous_acceleration);
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
            const array_1d<double, 3 > & previous_displacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
			const array_1d<double, 3 > & previous_velocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);
            const array_1d<double, 3 > & previous_acceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

            array_1d<double, 3 > & current_displacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);
            //array_1d<double, 3 > & ImposedDisplacement  = (i)->FastGetSolutionStepValue(IMPOSED_DISPLACEMENT);

            if ((i->pGetDof(DISPLACEMENT_X))->IsFixed() == false)
            {
                current_displacement[0] = 0.0;
            }
            else
            {
                current_displacement[0]  = previous_displacement[0];// + ImposedDisplacement[0];
            }

            if (i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false)
            {
                current_displacement[1] = 0.0;
            }
            else
            {
                current_displacement[1]  = previous_displacement[1];// + ImposedDisplacement[1];
            }

            if (i->HasDofFor(DISPLACEMENT_Z))
            {
                if (i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false)
                {
                    current_displacement[2] = 0.0;
                }
                else
                {
                    current_displacement[2]  = previous_displacement[2];// + ImposedDisplacement[2];
                }
            }

            if (i->HasDofFor(PRESSURE))
            {
                double& current_pressure        = (i)->FastGetSolutionStepValue(PRESSURE);
                const double& previous_pressure = (i)->FastGetSolutionStepValue(PRESSURE, 1);

                if ((i->pGetDof(PRESSURE))->IsFixed() == false)
                    current_pressure = previous_pressure;
            }

            // Updating time derivatives
            array_1d<double, 3 > & current_velocity       = (i)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & current_acceleration   = (i)->FastGetSolutionStepValue(ACCELERATION);

            UpdateVelocity(current_velocity, current_displacement, previous_velocity, previous_acceleration);
            UpdateAcceleration (current_acceleration, current_displacement, previous_velocity, previous_acceleration);

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

        int num_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector element_partition;
        OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(), num_threads, element_partition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ElementsArrayType::iterator element_begin = rModelPart.Elements().begin() + element_partition[k];
            ElementsArrayType::iterator element_end   = rModelPart.Elements().begin() + element_partition[k + 1];

            for (ElementsArrayType::iterator itElem = element_begin; itElem != element_end; itElem++)
            {
                itElem->Initialize(); // function to initialize the element
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

        int num_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector condition_partition;
        OpenMPUtils::DivideInPartitions(rModelPart.Conditions().size(), num_threads, condition_partition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ConditionsArrayType::iterator condition_begin = rModelPart.Conditions().begin() + condition_partition[k];
            ConditionsArrayType::iterator condition_end   = rModelPart.Conditions().begin() + condition_partition[k + 1];

            for (ConditionsArrayType::iterator itCond = condition_begin; itCond != condition_end; itCond++)
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
                array_1d<double, 3 > & nodal_momentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                array_1d<double, 3 > & nodal_inertia  = (i)->FastGetSolutionStepValue(NODAL_INERTIA);
                double & nodal_mass     = (i)->FastGetSolutionStepValue(NODAL_MASS);
                double & nodal_pressure = (i)->FastGetSolutionStepValue(PRESSURE,1);

                double & nodal_density = (i)->FastGetSolutionStepValue(DENSITY);
                double & nodal_aux_R   = (i)->FastGetSolutionStepValue(AUX_R);
                array_1d<double, 3 > & nodal_aux_R_vel = (i)->FastGetSolutionStepValue(AUX_R_VEL);
                array_1d<double, 3 > & nodal_aux_R_acc = (i)->FastGetSolutionStepValue(AUX_R_ACC);

                nodal_momentum.clear();
                nodal_inertia.clear();
                nodal_mass= 0.0;
                nodal_pressure = 0.0;

                nodal_density = 0.0;
                nodal_aux_R = 0.0;
                nodal_aux_R_vel.clear();
                nodal_aux_R_acc.clear();
            }

            if((i)->SolutionStepsDataHas(DISPLACEMENT) && (i)->SolutionStepsDataHas(VELOCITY) && (i)->SolutionStepsDataHas(ACCELERATION) )
            {
                array_1d<double, 3 > & nodal_displacement = (i)->FastGetSolutionStepValue(DISPLACEMENT);
                double & nodal_pressure = (i)->FastGetSolutionStepValue(PRESSURE);
                array_1d<double, 3 > & nodal_velocity     = (i)->FastGetSolutionStepValue(VELOCITY,1);
                array_1d<double, 3 > & nodal_acceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);
                array_1d<double, 3 > & delta_nodal_velocity     = (i)->FastGetSolutionStepValue(AUX_VELOCITY);
                array_1d<double, 3 > & delta_nodal_acceleration = (i)->FastGetSolutionStepValue(AUX_ACCELERATION);

                nodal_displacement.clear();
                nodal_pressure = 0.0;
                nodal_velocity.clear();
                nodal_acceleration.clear();
                delta_nodal_velocity.clear();
                delta_nodal_acceleration.clear();
            }
		}

        double norm_velocity = 1.0;
        double norm_acceleration = 1.0;
        double norm_pressure = 1.0;
        double norm_delta_velocity = 1.0;
        double norm_delta_acceleration = 1.0;
        double norm_delta_pressure = 1.0;

        int it_num = 1;

        while (it_num <2)
        {
            #pragma omp parallel for
			for( int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
			{
                auto i = mr_grid_model_part.NodesBegin() + iter;
                if( (i)->SolutionStepsDataHas(NODAL_MOMENTUM) && (i)->SolutionStepsDataHas(NODAL_MASS) && (i)->SolutionStepsDataHas(NODAL_INERTIA))//&& (i)->SolutionStepsDataHas(NODAL_INTERNAL_FORCE) )
                {
                    array_1d<double, 3 > & nodal_momentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                    array_1d<double, 3 > & nodal_inertia = (i)->FastGetSolutionStepValue(NODAL_INERTIA);
                    array_1d<double, 3 > & delta_nodal_velocity = (i)->FastGetSolutionStepValue(AUX_VELOCITY,1);
                    array_1d<double, 3 > & delta_nodal_acceleration = (i)->FastGetSolutionStepValue(AUX_ACCELERATION,1);

                    double & nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
                    nodal_momentum.clear();
                    nodal_inertia.clear();
                    delta_nodal_velocity.clear();
                    delta_nodal_acceleration.clear();

                    nodal_mass = 0.0;

                    if(i->SolutionStepsDataHas(NODAL_MPRESSURE)) {
                        double & nodal_mpressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);
                        nodal_mpressure = 0.0;
                    }
                }
			}

            norm_velocity = 0.0;
            norm_acceleration = 0.0;
            norm_pressure = 0.0;
            norm_delta_velocity = 0.0;
            norm_delta_acceleration = 0.0;
            norm_delta_pressure = 0.0;

            Scheme<TSparseSpace,TDenseSpace>::InitializeSolutionStep(r_model_part,A,Dx,b);

			#pragma omp parallel for reduction(+:norm_velocity,norm_acceleration,norm_pressure,norm_delta_velocity,norm_delta_acceleration,norm_delta_pressure)
			for(int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
			{

			    auto i = mr_grid_model_part.NodesBegin() + iter;
			    const double & nodal_mass     = (i)->FastGetSolutionStepValue(NODAL_MASS);

                if (nodal_mass > 1.0e-16 )
                {
                    array_1d<double, 3 > & delta_nodal_velocity     = (i)->FastGetSolutionStepValue(AUX_VELOCITY,1);
                    array_1d<double, 3 > & delta_nodal_acceleration = (i)->FastGetSolutionStepValue(AUX_ACCELERATION,1);

                    const array_1d<double, 3 > & nodal_momentum   = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                    const array_1d<double, 3 > & nodal_inertia    = (i)->FastGetSolutionStepValue(NODAL_INERTIA);

                    array_1d<double, 3 > & nodal_velocity     = (i)->FastGetSolutionStepValue(VELOCITY,1);
                    array_1d<double, 3 > & nodal_acceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);
                    double & nodal_pressure = (i)->FastGetSolutionStepValue(PRESSURE,1);

                    double delta_nodal_pressure = 0.0;

                    if (i->HasDofFor(PRESSURE) && i->SolutionStepsDataHas(NODAL_MPRESSURE))
                    {
                        double & nodal_mpressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);
                        delta_nodal_pressure = nodal_mpressure/nodal_mass;
                    }

                    delta_nodal_velocity = nodal_momentum/nodal_mass;
                    delta_nodal_acceleration = nodal_inertia/nodal_mass;

                    nodal_velocity += delta_nodal_velocity;
                    nodal_acceleration += delta_nodal_acceleration;

                    nodal_pressure += delta_nodal_pressure;

                    norm_delta_velocity += (delta_nodal_velocity[0]*delta_nodal_velocity[0]+delta_nodal_velocity[1]*delta_nodal_velocity[1]+delta_nodal_velocity[2]*delta_nodal_velocity[2]);
                    norm_delta_acceleration += (delta_nodal_acceleration[0]*delta_nodal_acceleration[0]+delta_nodal_acceleration[1]*delta_nodal_acceleration[1]+delta_nodal_acceleration[2]*delta_nodal_acceleration[2]);
                    norm_delta_pressure += (delta_nodal_pressure * delta_nodal_pressure);

                    norm_velocity += (nodal_velocity[0]*nodal_velocity[0]+nodal_velocity[1]*nodal_velocity[1]+nodal_velocity[2]*nodal_velocity[2]);
                    norm_acceleration += (nodal_acceleration[0]*nodal_acceleration[0]+nodal_acceleration[1]*nodal_acceleration[1]+nodal_acceleration[2]*nodal_acceleration[2]);
                    norm_pressure += (nodal_pressure * nodal_pressure);
                }
            }

            norm_velocity     = std::sqrt(norm_velocity);
            norm_acceleration = std::sqrt(norm_acceleration);
            norm_pressure     = std::sqrt(norm_pressure);

            norm_delta_velocity     = std::sqrt(norm_delta_velocity);
            norm_delta_acceleration = std::sqrt(norm_delta_acceleration);
            norm_delta_pressure     = std::sqrt(norm_delta_pressure);

            ++it_num;
        }

        double delta_time = CurrentProcessInfo[DELTA_TIME];
        KRATOS_ERROR_IF(delta_time == 0) << "Detected delta_time = 0 in the Solution Scheme ... check if the time step is created correctly for the current model part" << std::endl;

        //initializing Newmark constants
        mNewmark.c0 = ( 1.0 / (mNewmark.beta * delta_time * delta_time) );
        mNewmark.c1 = ( mNewmark.gamma / (mNewmark.beta * delta_time) );
        mNewmark.c2 = ( 1.0 / (mNewmark.beta * delta_time) );
        mNewmark.c3 = ( 0.5 / (mNewmark.beta) - 1.0 );
        mNewmark.c4 = ( (mNewmark.gamma / mNewmark.beta) - 1.0  );
        mNewmark.c5 = ( delta_time * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2 ) );

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

        int num_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector element_partition;
        OpenMPUtils::DivideInPartitions(rElements.size(), num_threads, element_partition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ElementsArrayType::iterator element_begin = rElements.begin() + element_partition[k];
            ElementsArrayType::iterator element_end   = rElements.begin() + element_partition[k + 1];

            for (ElementsArrayType::iterator itElem = element_begin; itElem != element_end; itElem++)
            {
                itElem->FinalizeSolutionStep(CurrentProcessInfo);

            }
        }

        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector condition_partition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), num_threads, condition_partition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ConditionsArrayType::iterator condition_begin = rConditions.begin() + condition_partition[k];
            ConditionsArrayType::iterator condition_end   = rConditions.begin() + condition_partition[k + 1];

            for (ConditionsArrayType::iterator itCond = condition_begin; itCond != condition_end; itCond++)
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
                array_1d<double, 3 > & external_force = (i)->FastGetSolutionStepValue(EXTERNAL_FORCE);
                array_1d<double, 3 > & internal_force = (i)->FastGetSolutionStepValue(INTERNAL_FORCE);
                external_force.clear();
                internal_force.clear();
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
    inline void UpdateVelocity(array_1d<double, 3 > & rCurrentVelocity,
                                      const array_1d<double, 3 > & rDeltaDisplacement,
                                      const array_1d<double, 3 > & rPreviousVelocity,
                                      const array_1d<double, 3 > & rPreviousAcceleration)
    {
        noalias(rCurrentVelocity) =  (mNewmark.c1 * rDeltaDisplacement - mNewmark.c4 * rPreviousVelocity
                                     - mNewmark.c5 * rPreviousAcceleration) * mNewmark.static_dynamic;
    }

    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************

    inline void UpdateAcceleration(array_1d<double, 3 > & rCurrentAcceleration,
                                   const array_1d<double, 3 > & rDeltaDisplacement,
                                   const array_1d<double, 3 > & rPreviousVelocity,
                                   const array_1d<double, 3 > & rPreviousAcceleration)
    {
        noalias(rCurrentAcceleration) =  (mNewmark.c0 * rDeltaDisplacement - mNewmark.c2 * rPreviousVelocity
                                         -  mNewmark.c3 * rPreviousAcceleration) * mNewmark.static_dynamic;
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


