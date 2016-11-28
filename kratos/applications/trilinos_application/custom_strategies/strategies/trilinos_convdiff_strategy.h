/* *********************************************************
 *
 *   Last Modified by:    $Author: rrossi $
 *   Date:                $Date: 2007-03-06 10:30:32 $
 *   Revision:            $Revision: 1.2 $
 *
 * ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY )
#define  KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/builder_and_solvers/trilinos_elimination_builder_and_solver.h"
// #include "convection_diffusion_application.h"
/*#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"*/
#include "includes/convection_diffusion_settings.h"
//#include "custom_utilities/convection_diffusion_settings.h"



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

/// This strategy is used to solve convection-diffusion problem.

/**   Detail class definition.

The convection-diffusion problem
        \f$ \rho C M \frac{\partial T}{\partial t} + \rho C S T  = - \kappa L T \f$ (1)

      is completed with the standard boundary conditions of prescribed temperature and prescribed normal
      heat flux in the thermal problem. For surfaces exposed to fire conditions, energy losses
      due to radiation and convection must be taken into account, and the thermal boundary condition is

      \f$ \kappa \frac{\partial T}{\partial n} + \overline{q_n} = 0 \f$ (2)
      where

      \f$ \overline{q_n} = q_n - \varepsilon \sigma (T^4 - T_0^4) - \alpha_c (T - T_0)\f$ (3)

Then, this strategy is employed to solve the following equation


Evaluates  \f$ L h s = \frac{\rho C}{\Delta t} (W, N) + (W, v. \nabla N) + \kappa (\nabla W, \nabla N) + 4 \epsilon \sigma T^3 \left\langle W, N \right\rangle + \alpha \left\langle
W, N \right\rangle \f$ and \f$ R h s = \rho (W, Q) + \frac{\rho C}{\Delta t} (W, T^n)- \left\langle W, q \right\rangle - \epsilon \sigma \left\langle W, T^4 - T_0^4 \right\rangle - \left\langle W, \alpha (T - T_0) \right\rangle - L h s \ast T \f$


 */
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class TrilinosConvectionDiffusionStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosConvectionDiffusionStrategy);
 
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;



    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /**
    @param model_part Reference to the ModelPart that contains the problem.
    @param pNewLinearSolver pointer to the solver for the temperature system.
    @paramReformDofAtEachIteration=true.
    @param time_order=2.
    @param prediction_order == 2.
     */


    TrilinosConvectionDiffusionStrategy(
        Epetra_MpiComm& Comm,
        ModelPart& model_part,
        typename TLinearSolver::Pointer pNewLinearSolver,
        bool ReformDofAtEachIteration = true,
        int time_order = 2,
        int prediction_order = 2
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, false)
    {
        KRATOS_TRY

        mtime_order = time_order;
        mprediction_order = prediction_order;

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
//            const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();




        //initializing fractional velocity solution step
        // 			typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
        // 			typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
        // 				( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

        bool CalculateReactions = false;
        bool CalculateNormDxFlag = true;

        // 			//choosing the solving strategy
        // //			mstep1 = typename BaseType::Pointer(
        // //				new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >
        // //				(model_part,pscheme,pNewLinearSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        // //			mstep1->SetEchoLevel(2);
        //
        //
        // 			typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;
        //
        // 			BuilderSolverTypePointer componentwise_build = BuilderSolverTypePointer(new	ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace,TDenseSpace,TLinearSolver,Variable<double> > (pNewLinearSolver,rUnknownVar) );
        // 			mstep1 = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > 				(model_part,pscheme,pNewLinearSolver,componentwise_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        // 			mstep1->SetEchoLevel(2);




        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
        typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

        //initializing fractional velocity solution step
        typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
        typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
                                               (new TrilinosResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());

        int guess_row_size = 30;

        //CONSTRUCTION OF VELOCITY
        BuilderSolverTypePointer builder = BuilderSolverTypePointer(new TrilinosResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (Comm, guess_row_size, pNewLinearSolver));
        this->mstep1 = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (model_part, pscheme, pNewLinearSolver, builder, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
        this->mstep1->SetEchoLevel(0);

        KRATOS_CATCH("")
    }

    /** Destructor.
     */
    virtual ~TrilinosConvectionDiffusionStrategy()
    {
    }

    /** Destructor.
     */

    //*********************************************************************************
    //**********************************************************************

    double Solve()
    {
        KRATOS_TRY

        //calculate the BDF coefficients
        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        double Dt = rCurrentProcessInfo[DELTA_TIME];

        if (mtime_order == 2)
        {
            if (BaseType::GetModelPart().GetBufferSize() < 3)
                KRATOS_THROW_ERROR(std::logic_error, "insufficient buffer size for BDF2", "")


                if (BaseType::GetModelPart().GetBufferSize() < 3)
                    KRATOS_THROW_ERROR(std::logic_error, "insufficient buffer size for BDF2", "")

                    double dt_old = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

            double rho = dt_old / Dt;
            double coeff = 1.0 / (Dt * rho * rho + Dt * rho);

            rCurrentProcessInfo[BDF_COEFFICIENTS].resize(3);
            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs[0] = coeff * (rho * rho + 2.0 * rho); //coefficient for step n+1
            BDFcoeffs[1] = -coeff * (rho * rho + 2.0 * rho + 1.0); //coefficient for step n
            BDFcoeffs[2] = coeff;
        }
        else
        {
            rCurrentProcessInfo[BDF_COEFFICIENTS].resize(2);
            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs[0] = 1.0 / Dt; //coefficient for step n+1
            BDFcoeffs[1] = -1.0 / Dt; //coefficient for step n
        }

        //second order prediction for the velocity
        if (mprediction_order == 2)
        {
            if (BaseType::GetModelPart().GetBufferSize() < 3)
                KRATOS_THROW_ERROR(std::logic_error, "insufficient buffer size for BDF2", "")

                ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
            const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

            for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                    i != BaseType::GetModelPart().NodesEnd(); ++i)
            {
                i->FastGetSolutionStepValue(rUnknownVar) = 2.00 * i->FastGetSolutionStepValue(rUnknownVar, 1) - i->FastGetSolutionStepValue(rUnknownVar, 2);
            }

        }

        //SOLVING THE PROBLEM
        rCurrentProcessInfo[FRACTIONAL_STEP] = 1;

        CalculateProjection();
        double Dp_norm = mstep1->Solve();
//            CalculateProjection();


        //Dp_norm = mstep1->Solve();
        //CalculateProjection();
        //Dp_norm = mstep1->Solve();
        //CalculateProjection();

        return Dp_norm;
        KRATOS_CATCH("")
    }





    //******************************************************************************************************
    //******************************************************************************************************
    //calculation of projection

    void CalculateProjection()
    {
        KRATOS_TRY;

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        //			const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();

        //first of all set to zero the nodal variables to be updated nodally
        for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                i != BaseType::GetModelPart().NodesEnd(); ++i)
        {
            (i)->GetSolutionStepValue(TEMP_CONV_PROJ) = 0.00;
            (i)->GetSolutionStepValue(NODAL_AREA) = 0.00;
        }

        //add the elemental contributions for the calculation of the velocity
        //and the determination of the nodal area
        rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
        for (ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin();
                i != BaseType::GetModelPart().ElementsEnd(); ++i)
        {
            (i)->InitializeSolutionStep(rCurrentProcessInfo);
        }

        //do syncronization
        BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(TEMP_CONV_PROJ);
        BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(NODAL_AREA);

        //solve nodally for the velocity
        for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                i != BaseType::GetModelPart().NodesEnd(); ++i)
        {
            double& conv_proj = (i)->GetSolutionStepValue(TEMP_CONV_PROJ);
            double temp = 1.00 / (i)->GetSolutionStepValue(NODAL_AREA);
            conv_proj *= temp;
        }

        KRATOS_CATCH("")
    }

    virtual void SetEchoLevel(int Level)
    {
        mstep1->SetEchoLevel(Level);
    }

    virtual void Clear()
    {
        mstep1->Clear();
    }

    /*@} */
    /**@name Operators
     */
    /*@{ */

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
    typename BaseType::Pointer mstep1;
    int mtime_order;
    int mprediction_order;





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

    /** Copy constructor.
     */
    TrilinosConvectionDiffusionStrategy(const TrilinosConvectionDiffusionStrategy& Other);


    /*@} */

}; /* Class TrilinosConvectionDiffusionStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY  defined */

