// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "utilities/parallel_utilities.h"

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

                      \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


 */
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class ResidualBasedRelaxationScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    /**@name Type Definitions */
    /*@{ */

    //typedef Kratos::shared_ptr< ResidualBasedRelaxationScheme<TSparseSpace,TDenseSpace> > Pointer;
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedRelaxationScheme );

    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename Element::DofsVectorType DofsVectorType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ResidualBasedRelaxationScheme(double NewAlphaBossak, double damping_factor)
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        //default values for the Newmark Scheme
        mAlphaBossak = NewAlphaBossak;
        mBetaNewmark = 0.25 * pow((1.00 - mAlphaBossak), 2);
        mGammaNewmark = 0.5 - mAlphaBossak;

        mdamping_factor = damping_factor;

        // Allocate auxiliary memory
        const int num_threads = ParallelUtilities::GetNumThreads();
        mMass.resize(num_threads);
        mDamp.resize(num_threads);
        mvel.resize(num_threads);

        //std::cout << "using the Relaxation Time Integration Scheme" << std::endl;
    }

    /** Destructor.
     */
    ~ResidualBasedRelaxationScheme() override
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    /**
        Performing the update of the solution.
     */

    //***************************************************************************

    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY

        //update of displacement (by DOF)
        for (typename DofsArrayType::iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
        {
            if (i_dof->IsFree())
            {
                i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
            }
        }

        //updating time derivatives (nodally for efficiency)
        array_1d<double, 3 > DeltaDisp;
        for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                i != r_model_part.NodesEnd(); ++i)
        {

            noalias(DeltaDisp) = (i)->FastGetSolutionStepValue(DISPLACEMENT) - (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3 > & CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 0);
            array_1d<double, 3 > & OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3 > & CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
            array_1d<double, 3 > & OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

            UpdateVelocity(CurrentVelocity, DeltaDisp, OldVelocity, OldAcceleration);

            UpdateAcceleration(CurrentAcceleration, DeltaDisp, OldVelocity, OldAcceleration);

        }

        KRATOS_CATCH( "" )

    }


    //***************************************************************************
    //predicts the solution at the current step as
    // x = xold + vold * Dt

    void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        std::cout << "prediction" << std::endl;
        array_1d<double, 3 > DeltaDisp;
        double DeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];


        for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                i != r_model_part.NodesEnd(); ++i)
        {

            array_1d<double, 3 > & OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 1);
            array_1d<double, 3 > & OldDisp = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            //predicting displacement = OldDisplacement + OldVelocity * DeltaTime;
            //ATTENTION::: the prediction is performed only on free nodes
            array_1d<double, 3 > & CurrentDisp = (i)->FastGetSolutionStepValue(DISPLACEMENT);

            if ((i->pGetDof(DISPLACEMENT_X))->IsFixed() == false)
                (CurrentDisp[0]) = OldDisp[0] + DeltaTime * OldVelocity[0];
            if (i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false)
                (CurrentDisp[1]) = OldDisp[1] + DeltaTime * OldVelocity[1];
            if (i->HasDofFor(DISPLACEMENT_Z))
                if (i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false)
                    (CurrentDisp[2]) = OldDisp[2] + DeltaTime * OldVelocity[2];

            //updating time derivatives ::: please note that displacements and its time derivatives
            //can not be consistently fixed separately
            noalias(DeltaDisp) = CurrentDisp - OldDisp;
            array_1d<double, 3 > & OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

            array_1d<double, 3 > & CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION);

            UpdateVelocity(CurrentVelocity, DeltaDisp, OldVelocity, OldAcceleration);

            UpdateAcceleration(CurrentAcceleration, DeltaDisp, OldVelocity, OldAcceleration);
        }

    }


    //***************************************************************************

    /** this function is designed to be called in the builder and solver
    to introduce
    the selected time integration scheme. It "asks" the matrix needed to
    the element and
    performs the operations needed to introduce the seected time
    integration scheme.

      this function calculates at the same time the contribution to the
    LHS and to the RHS
      of the system
     */
    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY
        int k = OpenMPUtils::ThisThread();
        //Initializing the non linear iteration for the current element
        //basic operations for the element considered
        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        rCurrentElement.CalculateMassMatrix(mMass[k], CurrentProcessInfo);
        rCurrentElement.CalculateDampingMatrix(mDamp[k], CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
        //adding the dynamic contributions (statics is already included)

        AddDynamicsToLHS(LHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

        AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

        KRATOS_CATCH( "" )

    }

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        int k = OpenMPUtils::ThisThread();
        //Initializing the non linear iteration for the current element

        //basic operations for the element considered
        rCurrentElement.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
        rCurrentElement.CalculateMassMatrix(mMass[k], CurrentProcessInfo);
        rCurrentElement.CalculateDampingMatrix(mDamp[k], CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        //adding the dynamic contributions (static is already included)
        AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

    }

    /** functions totally analogous to the precedent but applied to
    the "condition" objects
     */
    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY
        int k = OpenMPUtils::ThisThread();
        rCurrentCondition.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        rCurrentCondition.CalculateMassMatrix(mMass[k], CurrentProcessInfo);
        rCurrentCondition.CalculateDampingMatrix(mDamp[k], CurrentProcessInfo);
        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        AddDynamicsToLHS(LHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

        AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY
        int k = OpenMPUtils::ThisThread();
        //Initializing the non linear iteration for the current condition

        //basic operations for the element considered
        rCurrentCondition.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
        rCurrentCondition.CalculateMassMatrix(mMass[k], CurrentProcessInfo);
        rCurrentCondition.CalculateDampingMatrix(mDamp[k], CurrentProcessInfo);
        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        //adding the dynamic contributions (static is already included)

        AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        Scheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);

        double DeltaTime = CurrentProcessInfo[DELTA_TIME];

        if (DeltaTime == 0)
            KRATOS_THROW_ERROR( std::logic_error, "detected delta_time = 0 in the Bossak Scheme ... check if the time step is created correctly for the current model part", "" )

        //initializing constants
        ma0 = 1.0 / (mBetaNewmark * pow(DeltaTime, 2));
        ma1 = mGammaNewmark / (mBetaNewmark * DeltaTime);
        ma2 = 1.0 / (mBetaNewmark * DeltaTime);
        ma3 = 1.0 / (2.0 * mBetaNewmark) - 1.0;
        ma4 = mGammaNewmark / mBetaNewmark - 1.0;
        ma5 = DeltaTime * 0.5 * (mGammaNewmark / mBetaNewmark - 2.0);
        mam = (1.0 - mAlphaBossak) / (mBetaNewmark * pow(DeltaTime, 2));
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    int Check(const ModelPart& r_model_part) const override
    {
        KRATOS_TRY

        int err = Scheme<TSparseSpace, TDenseSpace>::Check(r_model_part);
        if (err != 0) return err;

        //check that variables are correctly allocated
        for (const auto& r_node : r_model_part.Nodes())
        {
            if (r_node.SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", r_node.Id() )
            if (r_node.SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", r_node.Id() )
            if (r_node.SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", r_node.Id() )
        }

        //check that dofs exist
        for (const auto& r_node : r_model_part.Nodes())
        {
            if (r_node.HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_THROW_ERROR( std::invalid_argument, "missing DISPLACEMENT_X dof on node ", r_node.Id() )
            if (r_node.HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR( std::invalid_argument, "missing DISPLACEMENT_Y dof on node ", r_node.Id() )
            if (r_node.HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR( std::invalid_argument, "missing DISPLACEMENT_Z dof on node ", r_node.Id() )
        }


        //check for admissible value of the AlphaBossak
        if (mAlphaBossak > 0.0 || mAlphaBossak < -0.3)
            KRATOS_THROW_ERROR( std::logic_error, "Value not admissible for AlphaBossak. Admissible values should be between 0.0 and -0.3. Current value is ", mAlphaBossak )

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


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */
    double mAlphaBossak;
    double mBetaNewmark;
    double mGammaNewmark;

    double ma0;
    double ma1;
    double ma2;
    double ma3;
    double ma4;
    double ma5;
    double mam;

    std::vector< Matrix >mMass;
    std::vector< Matrix >mDamp;
    std::vector< Vector >mvel;

    double mdamping_factor;


    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    //*********************************************************************************
    //Updating first time Derivative

    //*********************************************************************************

    inline void UpdateVelocity(array_1d<double, 3 > & CurrentVelocity,
                               const array_1d<double, 3 > & DeltaDisp,
                               const array_1d<double, 3 > & OldVelocity,
                               const array_1d<double, 3 > & OldAcceleration)
    {
        noalias(CurrentVelocity) = ma1 * DeltaDisp - ma4 * OldVelocity - ma5*OldAcceleration;
    }



    //**************************************************************************

    inline void UpdateAcceleration(array_1d<double, 3 > & CurrentAcceleration,
                                   const array_1d<double, 3 > & DeltaDisp,
                                   const array_1d<double, 3 > & OldVelocity,
                                   const array_1d<double, 3 > & OldAcceleration)
    {

        noalias(CurrentAcceleration) = ma0 * DeltaDisp - ma2 * OldVelocity - ma3*OldAcceleration;
    }




    //****************************************************************************

    /**
    Kdyn = a0*M + a1*D + K

     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& CurrentProcessInfo)
    {
        // adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += (mdamping_factor * ma1) * M;
        }
    }

    //****************************************************************************

    /**
    bdyn = b - M*acc - D*vel

     */
    void AddDynamicsToRHS(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& CurrentProcessInfo)
    {
        //adding inertia contribution
        if (M.size1() != 0)
        {
            int k = OpenMPUtils::ThisThread();
            const auto& r_const_elem_ref = rCurrentElement;
            r_const_elem_ref.GetFirstDerivativesVector(mvel[k], 0);
            noalias(RHS_Contribution) -= mdamping_factor * prod(M, mvel[k]);
        }

    }

    void AddDynamicsToRHS(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& CurrentProcessInfo)
    {
        //adding inertia contribution - DO NOT ADD
        if (M.size1() != 0)
        {
            int k = OpenMPUtils::ThisThread();
            const auto& r_const_cond_ref = rCurrentCondition;
            r_const_cond_ref.GetFirstDerivativesVector(mvel[k], 0);
            noalias(RHS_Contribution) -= mdamping_factor * prod(M, mvel[k]);
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

} /* namespace Kratos.*/

