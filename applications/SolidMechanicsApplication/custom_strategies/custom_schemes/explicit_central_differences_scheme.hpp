//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_EXPLICIT_CENTRAL_DIFFERENCES)
#define  KRATOS_EXPLICIT_CENTRAL_DIFFERENCES


/* System includes */


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
        const double  max_delta_time,
        const double  fraction_delta_time,
        const double  time_step_prediction_level,
        const bool    mRayleighDamping
    )
        : Scheme<TSparseSpace,TDenseSpace>()
    {

        mtimestep                   = 0.00;
        mCalculateOldTime           = false;
        mTimeStepPredictionLevel    = time_step_prediction_level;
        mmax_delta_time             = max_delta_time;
        mfraction_delta_time        = fraction_delta_time;
        mdamping_coeficients        = false;
        mInitialConditions          = true;

        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();

        //mMatrix.M.resize(NumThreads);
        mMatrix.D.resize(NumThreads);

        mVector.v.resize(NumThreads);
        //mVector.a.resize(NumThreads);
        //mVector.ap.resize(NumThreads);


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
             KRATOS_ERROR(std::logic_error, "Insufficient buffer size for Central Difference Scheme. It has to be 2", "")
        }

        return 0;
        KRATOS_CATCH("");
    }

    virtual void Initialize(
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY

        mModelPart = r_model_part;

        if(mTimeStepPredictionLevel>0)
        {

           CalculateDeltaTime(r_model_part);
        }


        ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

        molddelta_time = rCurrentProcessInfo[DELTA_TIME]; //initializing the variable.

        mSchemeIsInitialized = true;

        mold_time = 0.0;
        mold_mid_step_time = 0.0;

        KRATOS_CATCH("")
    }

    //***************************************************************************

    virtual void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY
        BaseType::InitializeSolutionStep(r_model_part,A,Dx,b);

       if(mTimeStepPredictionLevel>1)
       {
            CalculateDeltaTime(r_model_part);
       }

       KRATOS_CATCH("")
    }

    //***************************************************************************


    void CalculateDeltaTime(ModelPart& r_model_part)
    {

        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements     = r_model_part.Elements();

#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;

#endif

        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

        double safety_factor = 0.65;  //most autors recommend a value near 0.80 (Belytschko - Nonlinear FE.. 2000. chap 6. pag. 315)

        Vector dts(number_of_threads);
        double delta_time_a = 0.00;
        for(int i = 0; i < number_of_threads; i++)
            dts[i] = mmax_delta_time/safety_factor;

        #pragma omp parallel for private(delta_time_a)
        for(int k=0; k<number_of_threads; k++)
        {
            typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for(ElementsArrayType::iterator it=it_begin; it!= it_end; it++)
            {

                //The following procedure is similar to the one followed in Abaqus for central difference scheme with damping.

                double length   = it->GetGeometry().Length();
                double alpha    = it->GetProperties()[RAYLEIGH_ALPHA];
                double beta     = it->GetProperties()[RAYLEIGH_BETA];
                double E        = it->GetProperties()[YOUNG_MODULUS];
                double v        = it->GetProperties()[POISSON_RATIO];
                double ro       = it->GetProperties()[DENSITY];
                double bulk     = E/(3*(1-2*v));
                
                double wavespeed = sqrt(bulk/ro);
                double w = 2*wavespeed/length; //freq

                double psi = 0.5*(alpha/w + beta*w); //critical ratio;

                delta_time_a = (2/w)*(sqrt(1+psi*psi)-psi);

                if(delta_time_a>0.00)
                    if(delta_time_a < dts[k])
                        dts[k] = delta_time_a;
            }
        }

        double new_min_time_step = *std::min_element(dts.begin(), dts.end());

        new_min_time_step *= safety_factor;
        
        if(new_min_time_step < mmax_delta_time)
        {
        
        rCurrentProcessInfo[DELTA_TIME] = new_min_time_step;

        std::cout<< "  New computed stable Time step                                    = "<< new_min_time_step << "         " << std::endl;
        std::cout<< "   " << std::endl;
        }
        
        KRATOS_CATCH("")
    }

    //***************************************************************************

    void ComputeInitialConditions(ModelPart& r_model_part)
    {
        KRATOS_TRY

        //Set Vel to zero

        NodesArrayType& pNodes           = r_model_part.Nodes();

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
                    array_1d<double,3>& mid_pos_velocity      = i->FastGetSolutionStepValue(MID_POS_VELOCITY);
                    array_1d<double,3>& current_velocity      = i->FastGetSolutionStepValue(VELOCITY);
                    array_1d<double,3>& current_Rhs           = i->FastGetSolutionStepValue(RHS);
                    array_1d<double,3>& current_displacement  = i->FastGetSolutionStepValue(DISPLACEMENT);

                    for (unsigned int j =0; j<3; j++)
                    {

                        mid_pos_velocity[j] = 0.0;
                        current_velocity[j] = 0.0;
                        current_Rhs[j] = 0.0;
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
    virtual void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY


        if(mInitialConditions)
        {
            ComputeInitialConditions(r_model_part);
            mInitialConditions  = false;
        }

        else
        {
            Central_differences(r_model_part);
        }

        KRATOS_CATCH( "" )
    }


    /** this function is designed to be called in the builder and solver to introduce
    the selected time integration scheme. It "asks" the matrix needed to the element and
    performs the operations needed to introduce the seected time integration scheme.

    this function calculates at the same time the contribution to the LHS and to the RHS
    of the system
    */
    virtual void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
    )
    {
        KRATOS_TRY
        //Initializing the non linear iteration for the current element
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentElement)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        (rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }


    //***************************************************************************
    //***************************************************************************

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
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

            AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread],/* mMatrix.M[thread],*/ rCurrentProcessInfo);

        }


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


            rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue(MID_POS_VELOCITY);

            rValues[index]     = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MID_POS_VELOCITY )[0];
            rValues[index + 1] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MID_POS_VELOCITY )[1];

            if ( dimension == 3 )
                rValues[index + 2] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MID_POS_VELOCITY )[2];

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

            rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue(MID_POS_VELOCITY);

            rValues[index]     = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MID_POS_VELOCITY )[0];
            rValues[index + 1] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MID_POS_VELOCITY )[1];

            if ( dimension == 3 )
                rValues[index + 2] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MID_POS_VELOCITY )[2];

        }
    }



    //***************************************************************************
    //***************************************************************************

    /** functions totally analogous to the precedent but applied to
    the "condition" objects
    */
    virtual void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
        (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);
        KRATOS_CATCH( "" )
    }

    virtual void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

 //       int thread = OpenMPUtils::ThisThread();

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);


//        if(mRayleighDamping)
//        {
//            (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread],rCurrentProcessInfo);

//            AddDynamicsToRHS (rCurrentCondition, RHS_Contribution, mMatrix.D[thread],/* mMatrix.M[thread],*/ rCurrentProcessInfo);
//        }


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


    struct  GeneralMatrices
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


    GeneralMatrices     mMatrix;
    GeneralVectors      mVector;
    
    double mmax_delta_time;
    double mfraction_delta_time;
    
    double mtimestep;
    
    bool mdamping_coeficients;
    
    bool mCalculateOldTime;
    
    double molddelta_time;

    bool mSchemeIsInitialized;

    bool mInitialConditions;

    double mTimeStepPredictionLevel;

    ModelPart mModelPart;

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    /*@} */
    /**@name Protected Operations*/
    /*@{ */

 
    void Central_differences( ModelPart& r_model_part ) //(Belytschko pg.313 Explicit methods)
    {

        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
        NodesArrayType& pNodes           = r_model_part.Nodes();

        //Old Step available variables;
        /*
        - mold_time;                       //old t N+1
        - mold_mid_step_time;              //old t N+1/2
        - mid_pos_velocity;                //old v N+1/2 which is used and then updated
        */
        //Step Update

        const double current_time = rCurrentProcessInfo[TIME];  //new N+1 //t N + current delta_time
        const double mid_step_time    = 0.5*( mold_time + current_time);
        const double current_delta_time = rCurrentProcessInfo[DELTA_TIME]; //new Delta_t N+1/2

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
                array_1d<double,3>& current_Rhs             = i->FastGetSolutionStepValue(RHS);
                array_1d<double,3>& current_velocity        = i->FastGetSolutionStepValue(VELOCITY);
                array_1d<double,3>& current_displacement    = i->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double,3>& mid_pos_velocity        = i->FastGetSolutionStepValue(MID_POS_VELOCITY);

                array_1d<double,3>& current_acceleration    = i->FastGetSolutionStepValue(ACCELERATION);

                current_acceleration[0] = 0.0;
                current_acceleration[1] = 0.0;
                current_acceleration[2] = 0.0;

                array_1d<double,3> acceleration = current_Rhs/nodal_mass;
                
                if( (i->pGetDof(DISPLACEMENT_X))->IsFree() )
                {
                                        //old_mid_pos_vel
                    current_velocity[0] = mid_pos_velocity[0] + (mold_time - mold_mid_step_time) * acceleration[0]; //+ actual_velocity;

                    //STEP UPDATE N <-- N+1

                    mid_pos_velocity[0] = current_velocity[0] + ( mid_step_time - mold_time) * acceleration[0] ;

                    current_displacement[0]  = current_displacement[0] + current_delta_time * mid_pos_velocity[0];

                    current_acceleration[0]  = acceleration[0]; //export acceleration on free nodes


                }


                if( (i->pGetDof(DISPLACEMENT_Y))->IsFree() )
                {

                    current_velocity[1] = mid_pos_velocity[1] + (mold_time - mold_mid_step_time) * acceleration[1]; //+ actual_velocity;

                    //STEP UPDATE N <-- N+1

                    mid_pos_velocity[1] = current_velocity[1] + ( mid_step_time - mold_time) * acceleration[1] ;

                    current_displacement[1]  = current_displacement[1] + current_delta_time * mid_pos_velocity[1];

                    current_acceleration[1]  = acceleration[1]; //export acceleration on free nodes


                }


                if( i->HasDofFor(DISPLACEMENT_Z))
                {

                    if( (i->pGetDof(DISPLACEMENT_Z))->IsFree() )
                    {

                        current_velocity[2] = mid_pos_velocity[2] + (mold_time - mold_mid_step_time) * acceleration[2]; //+ actual_velocity;

                        //STEP UPDATE N <-- N+1

                        mid_pos_velocity[2] = current_velocity[2] + ( mid_step_time - mold_time) * acceleration[2] ;

                        current_displacement[2]  = current_displacement[2] + current_delta_time * mid_pos_velocity[2];

                        current_acceleration[2]  = acceleration[2]; //export acceleration on free nodes

                    }

                }
            }
        }

         mold_time = current_time;
         mold_mid_step_time = mid_step_time;

        KRATOS_CATCH("")
    }


    void ForwardEuler(ModelPart& r_model_part )
    {

        KRATOS_TRY

        ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
        NodesArrayType& pNodes           = r_model_part.Nodes();

        const double current_delta_time  = CurrentProcessInfo[DELTA_TIME];

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



                array_1d<double, 3 > & vel             = i->FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3 > & displ           = i->FastGetSolutionStepValue(DISPLACEMENT);


                array_1d<double,3>& current_Rhs           = i->FastGetSolutionStepValue(RHS);              /// Fext(n) - Fint(n)

                const double& nodal_mass    =  i->FastGetSolutionStepValue(NODAL_MASS);


                double aux = current_delta_time / nodal_mass;


                if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == false )
                {
                    vel[0] += aux * current_Rhs[0];

                    displ[0] +=  current_delta_time * vel[0];

                }


                if( (i->pGetDof(DISPLACEMENT_Y))->IsFixed() == false )
                {
                    vel[1] += aux * current_Rhs[1];

                    displ[1] += current_delta_time * vel[1];

                }


                if( i->HasDofFor(DISPLACEMENT_Z))
                {
                    if( (i->pGetDof(DISPLACEMENT_Z))->IsFixed() == false )
                    {
                        vel[2] += aux * current_Rhs[2];

                        displ[2] +=   current_delta_time * vel[2];

                    }

                }
            }


        }


        KRATOS_CATCH("")
    }



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

     double mold_time;
     double mold_mid_step_time;
     bool mRayleighDamping;



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

#endif /* KRATOS_EXPLICIT_CENTRAL_DIFFERENCES  defined */

