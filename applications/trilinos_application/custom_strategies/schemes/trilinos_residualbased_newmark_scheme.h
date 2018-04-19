//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_TRILINOS_RESIDUALBASED_NEWMARK_SCHEME )
#define  KRATOS_TRILINOS_RESIDUALBASED_NEWMARK_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
#include "includes/legacy_structural_app_vars.h"

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
/** Short class definition.

This class provides the implementation of the basic tasks that are needed by the solution strategy.
It is intended to be the place for tailoring the solution strategies to problem specific tasks.

Detail class definition.

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}*/

template<class TSparseSpace,  class TDenseSpace >

class TrilinosResidualBasedNewmarkScheme: public Scheme<TSparseSpace, TDenseSpace>
{

public:
    /**
     * Type Definitions
     */
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosResidualBasedNewmarkScheme);

    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::ElementsArrayType ElementsArrayType;

    typedef typename Element::DofsVectorType DofsVectorType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    /**
     * Constructor.
     * @param mDissipationRadius if the scheme is numerically energy conserving or not
     *                              == 1.0 energy conserving
     *                              < 1.0 numerically discipating
     * @ref Chung&Hulbert: "A time integration algorithm for structural dynamics with improved
     *  numerical dissipation: The generalized alpha method" Journal of applied mechanics, 60, 371-375
     */
    TrilinosResidualBasedNewmarkScheme( double mDissipationRadius ): Scheme<TSparseSpace, TDenseSpace>()
    {
        mElementalDofList.resize( 5 );

//                      mAlpha= mDissipationRadius/(1.0+mDissipationRadius);

//          mAlpha_m= (2.0*mDissipationRadius-1.0)/(mDissipationRadius+1.0);
        //For pure Newmark Scheme
        mAlpha = 0.0;
        mAlpha_m = 0.0;

        mBeta = ( 1.0 + mAlpha - mAlpha_m ) * ( 1.0 + mAlpha - mAlpha_m ) / 4.0;

        mGamma = 0.5 + mAlpha - mAlpha_m;

        std::cout << "using the Generalized alpha Time Integration Scheme with radius= " << mDissipationRadius << " alpha_f= " << mAlpha << " alpha_m= " << mAlpha_m << " beta= " << mBeta << " gamma= " << mGamma << std::endl;
        std::cout << "PURE Newmark Time !!!!!!!!!!!!!!!!!!!!!" << std::endl;
        //(...)_NULL DOF at the begin of time step, (...)_EINS DOF at the end of time step, (...)
        // DOF at the midpoint of time step, Please recognize that the name of the DOF is (...)
        // while the iteration is done towards the (...)_EINS DOF value at end of time step
    }

    /** Destructor.*/
    virtual ~TrilinosResidualBasedNewmarkScheme
    () {}

    /**@name Operators */

    /** Performing the update of the solution.*/
    //***************************************************************************
    /**
     * incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param r_model_part
     * @param rDofSet set of all primary variables
     * @param A LHS matrix
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

        for ( ModelPart::NodeIterator i = r_model_part.NodesBegin() ;
                i != r_model_part.NodesEnd() ; i++ )
        {
            if ( i->HasDofFor( DISPLACEMENT_X ) )
            {
                if ( i->GetDof( DISPLACEMENT_X ).IsFree() )
                {
                    i->GetSolutionStepValue( DISPLACEMENT_EINS_X ) += *( Dx[i->GetDof( DISPLACEMENT_X ).EquationId()] );
                }
            }

            if ( i->HasDofFor( DISPLACEMENT_Y ) )
            {
                if ( i->GetDof( DISPLACEMENT_Y ).IsFree() )
                {
                    i->GetSolutionStepValue( DISPLACEMENT_EINS_Y ) += *( Dx[i->GetDof( DISPLACEMENT_Y ).EquationId()] );
                }
            }

            if ( i->HasDofFor( DISPLACEMENT_Z ) )
            {
                if ( i->GetDof( DISPLACEMENT_Z ).IsFree() )
                {
                    i->GetSolutionStepValue( DISPLACEMENT_EINS_Z ) += *( Dx[i->GetDof( DISPLACEMENT_Z ).EquationId()] );
                }
            }

            if ( i->HasDofFor( WATER_PRESSURE ) )
            {
                if ( i->GetDof( WATER_PRESSURE ).IsFree() )
                {
                    i->GetSolutionStepValue( WATER_PRESSURE_EINS ) += *( Dx[i->GetDof( WATER_PRESSURE ).EquationId()] );
                }
            }

            if ( i->HasDofFor( AIR_PRESSURE ) )
            {
                if ( i->GetDof( AIR_PRESSURE ).IsFree() )
                {
                    i->GetSolutionStepValue( AIR_PRESSURE_EINS ) += *( Dx[i->GetDof( AIR_PRESSURE ).EquationId()] );
                }
            }
        }

        KRATOS_CATCH( "" )
    }

    /**
     * initializes next newton step by calculating
     * u_{n+1-alpha}= u_n*alpha+U_n+1^k+1*(1-alpha)
     * @param r_model_part
     * @param A LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeNonLinIteration( ModelPart& r_model_part,
                                    TSystemMatrixType& A,
                                    TSystemVectorType& Dx,
                                    TSystemVectorType& b ) override
    {
        KRATOS_TRY

        ProcessInfo CurrentProcessInfo = r_model_part.GetProcessInfo();
        //Update nodal values and nodal velocities at mAlpha

        for ( ModelPart::NodeIterator i = r_model_part.NodesBegin() ;
                i != r_model_part.NodesEnd() ; i++ )
        {
            if ( i->HasDofFor( DISPLACEMENT_X ) )
            {
                i->GetSolutionStepValue( ACCELERATION_EINS_X )
                = 1.0 / ( mBeta * CurrentProcessInfo[DELTA_TIME]
                          * CurrentProcessInfo[DELTA_TIME] )
                  * ( i->GetSolutionStepValue( DISPLACEMENT_EINS_X )
                      - i->GetSolutionStepValue( DISPLACEMENT_NULL_X ) )
                  - 1.0 / ( mBeta * CurrentProcessInfo[DELTA_TIME] )
                  * i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_X )
                  - ( 1.0 - 2.0 * mBeta ) / ( 2.0 * mBeta ) * i->GetSolutionStepValue( ACCELERATION_NULL_X );

                i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_X )
                = ( i->GetSolutionStepValue( DISPLACEMENT_EINS_X )
                    - i->GetSolutionStepValue( DISPLACEMENT_NULL_X ) )
                  * mGamma / ( mBeta * CurrentProcessInfo[DELTA_TIME] )
                  - ( mGamma - mBeta ) / mBeta * ( i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_X ) )
                  - ( mGamma - 2.0 * mBeta ) / ( 2.0 * mBeta ) * CurrentProcessInfo[DELTA_TIME]
                  * ( i->GetSolutionStepValue( ACCELERATION_NULL_X ) );

                i->GetSolutionStepValue( ACCELERATION_X )
                = mAlpha_m * i->GetSolutionStepValue( ACCELERATION_NULL_X )
                  + ( 1.0 - mAlpha_m ) * i->GetSolutionStepValue( ACCELERATION_EINS_X );

                i->GetSolutionStepValue( DISPLACEMENT_DT_X )
                = mAlpha * i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_X )
                  + ( 1.0 - mAlpha ) * i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_X );

                i->GetSolutionStepValue( DISPLACEMENT_X )
                = mAlpha * i->GetSolutionStepValue( DISPLACEMENT_NULL_X )
                  + ( 1.0 - mAlpha ) * i->GetSolutionStepValue( DISPLACEMENT_EINS_X );

            }

            if ( i->HasDofFor( DISPLACEMENT_Y ) )
            {
                i->GetSolutionStepValue( ACCELERATION_EINS_Y )
                = 1.0 / ( mBeta * CurrentProcessInfo[DELTA_TIME] * CurrentProcessInfo[DELTA_TIME] )
                  * ( i->GetSolutionStepValue( DISPLACEMENT_EINS_Y )
                      - i->GetSolutionStepValue( DISPLACEMENT_NULL_Y ) )
                  - 1.0 / ( mBeta * CurrentProcessInfo[DELTA_TIME] )
                  * i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Y )
                  - ( 1.0 - 2.0 * mBeta ) / ( 2.0 * mBeta ) *
                  i->GetSolutionStepValue( ACCELERATION_NULL_Y );

                i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_Y )
                = ( i->GetSolutionStepValue( DISPLACEMENT_EINS_Y )
                    - i->GetSolutionStepValue( DISPLACEMENT_NULL_Y ) )
                  * mGamma / ( mBeta * CurrentProcessInfo[DELTA_TIME] )
                  - ( mGamma - mBeta ) / mBeta * ( i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Y ) )
                  - ( mGamma - 2.0 * mBeta ) / ( 2.0 * mBeta ) * CurrentProcessInfo[DELTA_TIME]
                  * ( i->GetSolutionStepValue( ACCELERATION_NULL_Y ) );

                i->GetSolutionStepValue( ACCELERATION_Y )
                = mAlpha_m * i->GetSolutionStepValue( ACCELERATION_NULL_Y )
                  + ( 1.0 - mAlpha_m ) * i->GetSolutionStepValue( ACCELERATION_EINS_Y );

                i->GetSolutionStepValue( DISPLACEMENT_DT_Y )
                = mAlpha * i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Y )
                  + ( 1.0 - mAlpha ) * i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_Y );

                i->GetSolutionStepValue( DISPLACEMENT_Y )
                = mAlpha * i->GetSolutionStepValue( DISPLACEMENT_NULL_Y ) + ( 1.0 - mAlpha ) *
                  i->GetSolutionStepValue( DISPLACEMENT_EINS_Y );
            }

            if ( i->HasDofFor( DISPLACEMENT_Z ) )
            {
                i->GetSolutionStepValue( ACCELERATION_EINS_Z )
                = 1.0 / ( mBeta * CurrentProcessInfo[DELTA_TIME] * CurrentProcessInfo[DELTA_TIME] )
                  * ( i->GetSolutionStepValue( DISPLACEMENT_EINS_Z )
                      - i->GetSolutionStepValue( DISPLACEMENT_NULL_Z ) )
                  - 1.0 / ( mBeta * CurrentProcessInfo[DELTA_TIME] )
                  * i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Z )
                  - ( 1.0 - 2.0 * mBeta ) / ( 2.0 * mBeta ) * i->GetSolutionStepValue( ACCELERATION_NULL_Z );

                i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_Z )
                = ( i->GetSolutionStepValue( DISPLACEMENT_EINS_Z )
                    - i->GetSolutionStepValue( DISPLACEMENT_NULL_Z ) )
                  * mGamma / ( mBeta * CurrentProcessInfo[DELTA_TIME] ) - ( mGamma - mBeta ) / mBeta *
                  ( i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Z ) )
                  - ( mGamma - 2.0 * mBeta ) / ( 2.0 * mBeta ) * CurrentProcessInfo[DELTA_TIME]
                  * ( i->GetSolutionStepValue( ACCELERATION_NULL_Z ) );

                i->GetSolutionStepValue( ACCELERATION_Z )
                = mAlpha_m * i->GetSolutionStepValue( ACCELERATION_NULL_Z )
                  + ( 1.0 - mAlpha_m ) * i->GetSolutionStepValue( ACCELERATION_EINS_Z );

                i->GetSolutionStepValue( DISPLACEMENT_DT_Z )
                = mAlpha * i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Z )
                  + ( 1.0 - mAlpha ) * i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_Z );

                i->GetSolutionStepValue( DISPLACEMENT_Z )
                = mAlpha * i->GetSolutionStepValue( DISPLACEMENT_NULL_Z )
                  + ( 1.0 - mAlpha ) * i->GetSolutionStepValue( DISPLACEMENT_EINS_Z );
            }

            if ( i->HasDofFor( WATER_PRESSURE ) )
            {
                i->GetSolutionStepValue( WATER_PRESSURE_EINS_ACCELERATION )
                = 1.0 / ( mBeta * CurrentProcessInfo[DELTA_TIME] * CurrentProcessInfo[DELTA_TIME] )
                  * ( i->GetSolutionStepValue( WATER_PRESSURE_EINS )
                      - i->GetSolutionStepValue( WATER_PRESSURE_NULL ) )
                  - 1.0 / ( mBeta * CurrentProcessInfo[DELTA_TIME] )
                  * i->GetSolutionStepValue( WATER_PRESSURE_NULL_DT )
                  - ( 1.0 - 2.0 * mBeta ) / ( 2.0 * mBeta ) * i->GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION );

                i->GetSolutionStepValue( WATER_PRESSURE_EINS_DT )
                = ( i->GetSolutionStepValue( WATER_PRESSURE_EINS )
                    - i->GetSolutionStepValue( WATER_PRESSURE_NULL ) )
                  * mGamma / ( mBeta * CurrentProcessInfo[DELTA_TIME] )
                  - ( mGamma - mBeta ) / mBeta * ( i->GetSolutionStepValue( WATER_PRESSURE_NULL_DT ) )
                  - ( mGamma - 2.0 * mBeta ) / ( 2.0 * mBeta ) * CurrentProcessInfo[DELTA_TIME]
                  * ( i->GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION ) );

                i->GetSolutionStepValue( WATER_PRESSURE_ACCELERATION )
                = mAlpha_m * i->GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION )
                  + ( 1.0 - mAlpha_m ) * i->GetSolutionStepValue( WATER_PRESSURE_EINS_ACCELERATION );

                i->GetSolutionStepValue( WATER_PRESSURE_DT )
                = mAlpha * i->GetSolutionStepValue( WATER_PRESSURE_NULL_DT )
                  + ( 1.0 - mAlpha ) * i->GetSolutionStepValue( WATER_PRESSURE_EINS_DT );

                i->GetSolutionStepValue( WATER_PRESSURE )
                = mAlpha * i->GetSolutionStepValue( WATER_PRESSURE_NULL )
                  + ( 1.0 - mAlpha ) * i->GetSolutionStepValue( WATER_PRESSURE_EINS );
            }

            if ( i->HasDofFor( AIR_PRESSURE ) )
            {
                i->GetSolutionStepValue( AIR_PRESSURE_EINS_ACCELERATION )
                = 1.0 / ( mBeta * CurrentProcessInfo[DELTA_TIME] * CurrentProcessInfo[DELTA_TIME] )
                  * ( i->GetSolutionStepValue( AIR_PRESSURE_EINS )
                      - i->GetSolutionStepValue( AIR_PRESSURE_NULL ) )
                  - 1.0 / ( mBeta * CurrentProcessInfo[DELTA_TIME] )
                  * i->GetSolutionStepValue( AIR_PRESSURE_NULL_DT )
                  - ( 1.0 - 2.0 * mBeta ) / ( 2.0 * mBeta ) * i->GetSolutionStepValue( AIR_PRESSURE_NULL_ACCELERATION );

                i->GetSolutionStepValue( AIR_PRESSURE_EINS_DT )
                = ( i->GetSolutionStepValue( AIR_PRESSURE_EINS )
                    - i->GetSolutionStepValue( AIR_PRESSURE_NULL ) )
                  * mGamma / ( mBeta * CurrentProcessInfo[DELTA_TIME] )
                  - ( mGamma - mBeta ) / mBeta *
                  ( i->GetSolutionStepValue( AIR_PRESSURE_NULL_DT ) )
                  - ( mGamma - 2.0 * mBeta ) / ( 2.0 * mBeta ) * CurrentProcessInfo[DELTA_TIME]
                  * ( i->GetSolutionStepValue( AIR_PRESSURE_NULL_ACCELERATION ) );

                i->GetSolutionStepValue( AIR_PRESSURE_ACCELERATION )
                = mAlpha_m * i->GetSolutionStepValue( AIR_PRESSURE_NULL_ACCELERATION )
                  + ( 1.0 - mAlpha_m ) * i->GetSolutionStepValue( AIR_PRESSURE_EINS_ACCELERATION );

                i->GetSolutionStepValue( AIR_PRESSURE_DT )
                = mAlpha * i->GetSolutionStepValue( AIR_PRESSURE_NULL_DT )
                  + ( 1.0 - mAlpha ) * i->GetSolutionStepValue( AIR_PRESSURE_EINS_DT );

                i->GetSolutionStepValue( AIR_PRESSURE )
                = mAlpha * i->GetSolutionStepValue( AIR_PRESSURE_NULL )
                  + ( 1.0 - mAlpha ) * i->GetSolutionStepValue( AIR_PRESSURE_EINS );
            }
        }

        //For total Lagrangian
        for ( ModelPart::NodeIterator i = r_model_part.NodesBegin() ;
                i != r_model_part.NodesEnd() ; ++i )
        {
            ( i )->X() = ( i )->X0() + i->GetSolutionStepValue( DISPLACEMENT_X );
            ( i )->Y() = ( i )->Y0() + i->GetSolutionStepValue( DISPLACEMENT_Y );
            ( i )->Z() = ( i )->Z0() + i->GetSolutionStepValue( DISPLACEMENT_Z );
        }

        KRATOS_CATCH( "" )
    }

    /**
     */
    void FinalizeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b ) override
    {
        KRATOS_TRY

        KRATOS_CATCH( "" )
    }

    /**
     * initializes time step solution
     * only for reasons if the time step solution is restarted
     * @param r_model_part
     * @param A LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b ) override
    {
        KRATOS_TRY

        ProcessInfo CurrentProcessInfo = r_model_part.GetProcessInfo();

        Scheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep( r_model_part, A, Dx, b );
        //Update nodal values and nodal velocities at mAlpha

        for ( ModelPart::NodeIterator i = r_model_part.NodesBegin() ;
                i != r_model_part.NodesEnd() ; i++ )
        {
            if ( i->HasDofFor( DISPLACEMENT_X ) &&  i->GetDof( DISPLACEMENT_X ).IsFree() )
            {
                i->GetSolutionStepValue( ACCELERATION_EINS_X ) =
                    i->GetSolutionStepValue( ACCELERATION_NULL_X );
                i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_X ) =
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_X );
                i->GetSolutionStepValue( DISPLACEMENT_EINS_X ) =
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_X );
            }

            if ( i->HasDofFor( DISPLACEMENT_Y ) &&  i->GetDof( DISPLACEMENT_Y ).IsFree() )
            {
                i->GetSolutionStepValue( ACCELERATION_EINS_Y ) =
                    i->GetSolutionStepValue( ACCELERATION_NULL_Y );
                i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_Y ) =
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Y );
                i->GetSolutionStepValue( DISPLACEMENT_EINS_Y ) =
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_Y );
            }

            if ( i->HasDofFor( DISPLACEMENT_Z ) && i->GetDof( DISPLACEMENT_Z ).IsFree() )
            {
                i->GetSolutionStepValue( ACCELERATION_EINS_Y ) =
                    i->GetSolutionStepValue( ACCELERATION_NULL_Z );
                i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_Z ) =
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Z );
                i->GetSolutionStepValue( DISPLACEMENT_EINS_Z ) =
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_Z );
            }

            if ( i->HasDofFor( WATER_PRESSURE ) && i->GetDof( WATER_PRESSURE ).IsFree() )
            {
                i->GetSolutionStepValue( WATER_PRESSURE_EINS_DT ) =
                    i->GetSolutionStepValue( WATER_PRESSURE_NULL_DT );
                i->GetSolutionStepValue( WATER_PRESSURE_EINS_ACCELERATION ) =
                    i->GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION );
                i->GetSolutionStepValue( WATER_PRESSURE_EINS ) =
                    i->GetSolutionStepValue( WATER_PRESSURE_NULL );
            }

            if ( i->HasDofFor( AIR_PRESSURE ) && i->GetDof( AIR_PRESSURE ).IsFree() )
            {
                i->GetSolutionStepValue( AIR_PRESSURE_EINS_DT ) =
                    i->GetSolutionStepValue( AIR_PRESSURE_NULL_DT );
                i->GetSolutionStepValue( AIR_PRESSURE_EINS_ACCELERATION ) =
                    i->GetSolutionStepValue( AIR_PRESSURE_NULL_ACCELERATION );
                i->GetSolutionStepValue( AIR_PRESSURE_EINS ) =
                    i->GetSolutionStepValue( AIR_PRESSURE_NULL );
            }
        }

        KRATOS_CATCH( "" )
    }

    /**
     * finalizes time step solution
     * by setting u_n= u_n+1^k etc.
     * u_{n+1-alpha}= u_n*alpha+U_n+1^k+1*(1-alpha)
     * @param r_model_part
     * @param A LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b ) override
    {

        InitializeNonLinIteration( r_model_part, A, Dx, b );

        ProcessInfo CurrentProcessInfo = r_model_part.GetProcessInfo();

        ElementsArrayType& pElements = r_model_part.Elements();

        for ( typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin();
                it != pElements.ptr_end(); ++it )
        {
            //calculate elemental contribution
            ( *it )->FinalizeSolutionStep( r_model_part.GetProcessInfo() );
        }

        if ( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
        {
            for ( ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; ++i )
            {
                if ( i->HasDofFor( DISPLACEMENT_X ) )
                {
                    i->GetSolutionStepValue( DISPLACEMENT_OLD_X )
                    = 0.0;
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_X )
                    = 0.0;
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_X )
                    = 0.0;
                    i->GetSolutionStepValue( ACCELERATION_NULL_X )
                    = 0.0;
                }

                if ( i->HasDofFor( DISPLACEMENT_Y ) )
                {
                    i->GetSolutionStepValue( DISPLACEMENT_OLD_Y )
                    = 0.0;
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_Y )
                    = 0.0;
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Y )
                    = 0.0;
                    i->GetSolutionStepValue( ACCELERATION_NULL_Y )
                    = 0.0;
                }

                if ( i->HasDofFor( DISPLACEMENT_Z ) )
                {
                    i->GetSolutionStepValue( DISPLACEMENT_OLD_Z )
                    = 0.0;
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_Z )
                    = 0.0;
                    i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Z )
                    = 0.0;
                    i->GetSolutionStepValue( ACCELERATION_NULL_Z )
                    = 0.0;
                }
            }
        }
        else
        {
            for ( ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; ++i )
            {
                if ( i->HasDofFor( DISPLACEMENT_X ) )
                {
                    i->GetSolutionStepValue( DISPLACEMENT_OLD_X ) = i->GetSolutionStepValue( DISPLACEMENT_X );

                    if ( CurrentProcessInfo[FIRST_TIME_STEP] )
                    {
                        i->GetSolutionStepValue( ACCELERATION_NULL_X ) = i->GetSolutionStepValue( ACCELERATION_X );
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_X ) = i->GetSolutionStepValue( DISPLACEMENT_DT_X );
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_X ) = i->GetSolutionStepValue( DISPLACEMENT_X );
                    }
                    else
                    {
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_X ) = i->GetSolutionStepValue( DISPLACEMENT_EINS_X );
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_X ) = i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_X );
                        i->GetSolutionStepValue( ACCELERATION_NULL_X ) = i->GetSolutionStepValue( ACCELERATION_EINS_X );
                    }
                }

                if ( i->HasDofFor( DISPLACEMENT_Y ) )
                {
                    i->GetSolutionStepValue( DISPLACEMENT_OLD_Y ) = i->GetSolutionStepValue( DISPLACEMENT_Y );

                    if ( CurrentProcessInfo[FIRST_TIME_STEP] )
                    {
                        i->GetSolutionStepValue( ACCELERATION_NULL_Y ) = i->GetSolutionStepValue( ACCELERATION_Y );
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Y ) = i->GetSolutionStepValue( DISPLACEMENT_DT_Y );
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_Y ) = i->GetSolutionStepValue( DISPLACEMENT_Y );
                    }
                    else
                    {
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_Y ) = i->GetSolutionStepValue( DISPLACEMENT_EINS_Y );
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Y ) = i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_Y );
                        i->GetSolutionStepValue( ACCELERATION_NULL_Y ) = i->GetSolutionStepValue( ACCELERATION_EINS_Y );
                    }
                }

                if ( i->HasDofFor( DISPLACEMENT_Z ) )
                {
                    i->GetSolutionStepValue( DISPLACEMENT_OLD_Z ) = i->GetSolutionStepValue( DISPLACEMENT_Z );

                    if ( CurrentProcessInfo[FIRST_TIME_STEP] )
                    {
                        i->GetSolutionStepValue( ACCELERATION_NULL_Z ) = i->GetSolutionStepValue( ACCELERATION_Z );
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Z ) = i->GetSolutionStepValue( DISPLACEMENT_DT_Z );
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_Z ) = i->GetSolutionStepValue( DISPLACEMENT_Z );
                    }
                    else
                    {
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_Z ) = i->GetSolutionStepValue( DISPLACEMENT_EINS_Z );
                        i->GetSolutionStepValue( DISPLACEMENT_NULL_DT_Z ) = i->GetSolutionStepValue( DISPLACEMENT_EINS_DT_Z );
                        i->GetSolutionStepValue( ACCELERATION_NULL_Z ) = i->GetSolutionStepValue( ACCELERATION_EINS_Z );
                    }
                }

                if ( i->HasDofFor( WATER_PRESSURE ) )
                {
                    if ( CurrentProcessInfo[FIRST_TIME_STEP] )
                    {
                        i->GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION ) = i->GetSolutionStepValue( WATER_PRESSURE_ACCELERATION );
                        i->GetSolutionStepValue( WATER_PRESSURE_NULL_DT ) = i->GetSolutionStepValue( WATER_PRESSURE_DT );
                        i->GetSolutionStepValue( WATER_PRESSURE_NULL ) = i->GetSolutionStepValue( WATER_PRESSURE );
                    }
                    else
                    {
                        i->GetSolutionStepValue( WATER_PRESSURE_NULL_DT )
                        = i->GetSolutionStepValue( WATER_PRESSURE_EINS_DT );
                        i->GetSolutionStepValue( WATER_PRESSURE_NULL )
                        = i->GetSolutionStepValue( WATER_PRESSURE_EINS );
                        i->GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION )
                        = i->GetSolutionStepValue( WATER_PRESSURE_EINS_ACCELERATION );
                    }
                }

                if ( i->HasDofFor( AIR_PRESSURE ) )
                {
                    if ( CurrentProcessInfo[FIRST_TIME_STEP] )
                    {
                        i->GetSolutionStepValue( AIR_PRESSURE_NULL_ACCELERATION )
                        = i->GetSolutionStepValue( AIR_PRESSURE_ACCELERATION );
                        i->GetSolutionStepValue( AIR_PRESSURE_NULL_DT )
                        = i->GetSolutionStepValue( AIR_PRESSURE_DT );
                        i->GetSolutionStepValue( AIR_PRESSURE_NULL )
                        = i->GetSolutionStepValue( AIR_PRESSURE );
                    }
                    else
                    {
                        i->GetSolutionStepValue( AIR_PRESSURE_NULL_DT )
                        = i->GetSolutionStepValue( AIR_PRESSURE_EINS_DT );
                        i->GetSolutionStepValue( AIR_PRESSURE_NULL )
                        = i->GetSolutionStepValue( AIR_PRESSURE_EINS );
                        i->GetSolutionStepValue( AIR_PRESSURE_NULL_ACCELERATION )
                        = i->GetSolutionStepValue( AIR_PRESSURE_EINS_ACCELERATION );
                    }
                }
            }
        }
    }

    //***************************************************************************
    //***************************************************************************

    /** this function is designed to be called in the builder and solver to introduce*/

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo ) override
    {
        KRATOS_TRY

        Matrix DampingMatrix;

        Matrix MassMatrix;

        ( rCurrentElement ) -> InitializeNonLinearIteration( CurrentProcessInfo );

        ( rCurrentElement )->
        CalculateLocalSystem( LHS_Contribution, RHS_Contribution, CurrentProcessInfo );

        ( rCurrentElement )->CalculateDampingMatrix( DampingMatrix, CurrentProcessInfo );

//          if(!(CurrentProcessInfo[QUASI_STATIC_ANALYSIS]))
//          {
//              (rCurrentElement)->CalculateMassMatrix(MassMatrix, CurrentProcessInfo);
//
//              AddDynamicsToRHS(rCurrentElement, RHS_Contribution, MassMatrix, CurrentProcessInfo);
//          }

        AssembleTimeSpaceLHS( rCurrentElement, LHS_Contribution,
                              DampingMatrix, MassMatrix, CurrentProcessInfo );

        ( rCurrentElement )->EquationIdVector( EquationId, CurrentProcessInfo );

        KRATOS_CATCH( "" )

    }

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo ) override
    {
        //Initializing the non linear iteration for the current element
        ( rCurrentElement ) -> InitializeNonLinearIteration( CurrentProcessInfo );
        //basic operations for the element considered
        ( rCurrentElement )->CalculateRightHandSide( RHS_Contribution, CurrentProcessInfo );
        ( rCurrentElement )->EquationIdVector( EquationId, CurrentProcessInfo );
    }


    /** functions totally analogous to the precedent but applied to
          the "condition" objects
    *       At the current status of implementation it does nothing
    */
    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo ) override
    {

        KRATOS_TRY

        Matrix DampingMatrix;

        Matrix MassMatrix;

        ( rCurrentCondition )-> CalculateLocalSystem( LHS_Contribution, RHS_Contribution, CurrentProcessInfo );

        ( rCurrentCondition )->CalculateDampingMatrix( DampingMatrix, CurrentProcessInfo );

        AssembleTimeSpaceLHS_Condition( rCurrentCondition, LHS_Contribution,
                                        DampingMatrix, MassMatrix, CurrentProcessInfo );

        ( rCurrentCondition )->EquationIdVector( EquationId, CurrentProcessInfo );

        KRATOS_CATCH( "" )
    }


    /**       At the current status of implementation it does nothing
     */

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo ) override
    {
        KRATOS_TRY
        ( rCurrentCondition ) ->
        CalculateRightHandSide( RHS_Contribution, CurrentProcessInfo );
        ( rCurrentCondition )->EquationIdVector( EquationId, CurrentProcessInfo );
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

    double mAlpha;
    double mAlpha_m;
    double mBeta;
    double mGamma;


    void AddDynamicsToRHS(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo )
    {
//          Vector acceleration;
//          //adding inertia contribution
//          if (M.size1() != 0)
//          {
//              rCurrentElement->GetSecondDerivativesVector(acceleration,0);
//              noalias(RHS_Contribution) -= prod(M, acceleration );
//          }
    }

    void AssembleTimeSpaceLHS(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& DampingMatrix,
        LocalSystemMatrixType& MassMatrix,
        ProcessInfo& CurrentProcessInfo )
    {
        // adding mass contribution to the dynamic stiffness
        for ( unsigned int prim = 0; prim < LHS_Contribution.size1(); prim++ )
            for ( unsigned int sec = 0; sec < LHS_Contribution.size2(); sec++ )
                LHS_Contribution( prim, sec ) = ( 1 - mAlpha ) * LHS_Contribution( prim, sec );

        if ( DampingMatrix.size1() == LHS_Contribution.size1()
                || DampingMatrix.size2() == LHS_Contribution.size2() )
        {
            for ( unsigned int prim = 0; prim < LHS_Contribution.size1(); prim++ )
                for ( unsigned int sec = 0; sec < LHS_Contribution.size2(); sec++ )
                    LHS_Contribution( prim, sec ) += DampingMatrix( prim, sec ) * ( 1 - mAlpha )
                                                     * mGamma / ( mBeta * CurrentProcessInfo[DELTA_TIME] );
        }


//          if(MassMatrix.size1() == LHS_Contribution.size1()
//              || MassMatrix.size2() == LHS_Contribution.size2())
//          {
//                 for(unsigned int prim=0; prim< LHS_Contribution.size1(); prim++)
//                       for(unsigned int sec=0; sec< LHS_Contribution.size2(); sec++)
//                         LHS_Contribution(prim,sec) += MassMatrix(prim,sec)*(1-mAlpha_m)
//                              *1/(mBeta*CurrentProcessInfo[DELTA_TIME]*CurrentProcessInfo[DELTA_TIME]);
//          }

    }


    void AssembleTimeSpaceLHS_Condition(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& DampingMatrix,
        LocalSystemMatrixType& MassMatrix,
        ProcessInfo& CurrentProcessInfo )
    {
        // adding mass contribution to the dynamic stiffness
        for ( unsigned int prim = 0; prim < LHS_Contribution.size1(); prim++ )
            for ( unsigned int sec = 0; sec < LHS_Contribution.size2(); sec++ )
                LHS_Contribution( prim, sec ) = ( 1 - mAlpha ) * LHS_Contribution( prim, sec );

        if ( DampingMatrix.size1() == LHS_Contribution.size1()
                || DampingMatrix.size2() == LHS_Contribution.size2() )
        {
            for ( unsigned int prim = 0; prim < LHS_Contribution.size1(); prim++ )
                for ( unsigned int sec = 0; sec < LHS_Contribution.size2(); sec++ )
                    LHS_Contribution( prim, sec ) += DampingMatrix( prim, sec ) * ( 1 - mAlpha )
                                                     * mGamma / ( mBeta * CurrentProcessInfo[DELTA_TIME] );
        }

//             if(MassMatrix.size1() == LHS_Contribution.size1()
//              || MassMatrix.size2() == LHS_Contribution.size2())
//          {
//                 for(unsigned int prim=0; prim< LHS_Contribution.size1(); prim++)
//                      for( unsigned int sec=0; sec< LHS_Contribution.size2(); sec++)
//                         LHS_Contribution(prim,sec) += MassMatrix(prim,sec)*(1-mAlpha_m)
//                              *1/(mBeta*CurrentProcessInfo[DELTA_TIME]*CurrentProcessInfo[DELTA_TIME]);
//          }
    }

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

private:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */
    DofsVectorType mElementalDofList;
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
}; /* Class Scheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_RESIDUALBASED_NEWMARK_SCHEME  defined */


