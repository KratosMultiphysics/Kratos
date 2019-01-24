// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

#if !defined(KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY_NONLINEAR )
#define  KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY_NONLINEAR

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "convection_diffusion_application.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "includes/convection_diffusion_settings.h"


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

/**   The convection-diffusion problem
        \f$ \rho C M \frac{\partial T}{\partial t} + \rho C S T  = - \kappa L T \f$ (1)

      is completed with the standard boundary conditions of prescribed temperature and prescribed normal
      heat flux in the thermal problem. For surfaces exposed to fire conditions, energy losses
      due to radiation and convection must be taken into account, and the thermal boundary condition is

      \f$ \kappa \frac{\partial T}{\partial n} + \overline{q_n} = 0 \f$ (2)
      where

      \f$ \overline{q_n} = q_n - \varepsilon \sigma (T^4 - T_0^4) - \alpha_c (T - T_0)\f$ (3)

Then, this strategy is employed to solve in iterative form the following equation

Evaluates  \f$ L h s = \frac{\rho C}{\Delta t} (W, N) + (W, v. \nabla N) + \kappa (\nabla W, \nabla N) + 4 \epsilon \sigma T^3 \left\langle W, N \right\rangle + \alpha \left\langle
W, N \right\rangle \f$ and \f$R h s = \rho (W, Q) + \frac{\rho C}{\Delta t} (W, T^n)- \left\langle W, q \right\rangle - \epsilon \sigma \left\langle W, T^4 -
T_0^4 \right\rangle - \left\langle W, \alpha (T - T_0) \right\rangle  - L h s \ast T \f$

until to achieve the convergence




 */
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class ResidualBasedConvectionDiffusionStrategyNonLinear
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedConvectionDiffusionStrategyNonLinear);

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

    /// Constructor

    /**
    @param model_part Reference to the ModelPart that contains the problem.
    @param pNewLinearSolver pointer to the solver for the temperature system.
    @paramReformDofAtEachIteration=true.
    @param time_order=2.
    @param int max_iter = 30.
    @param double toll = 1e-6.

     */


    ResidualBasedConvectionDiffusionStrategyNonLinear(
        ModelPart& model_part,
        typename TLinearSolver::Pointer pNewLinearSolver,
        bool ReformDofAtEachIteration = true,
        int time_order = 2,
        int max_iter = 30,
        double toll = 1e-6
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, false)
    {
        KRATOS_TRY

        mtime_order = time_order;
        mOldDt = 0.00;
        mtoll = toll;
        mmax_iter = max_iter;

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        //initializing fractional velocity solution step
        typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
        typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
                                               (new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());

        bool CalculateReactions = false;
        bool CalculateNormDxFlag = true;



        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;

        BuilderSolverTypePointer componentwise_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> > (pNewLinearSolver, rUnknownVar));
        mstep1 = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (model_part, pscheme, pNewLinearSolver, componentwise_build, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
        mstep1->SetEchoLevel(2);

        Check();
        KRATOS_CATCH("")
    }

    /** Destructor.
     */
    virtual ~ResidualBasedConvectionDiffusionStrategyNonLinear()
    {
    }

    /** Destructor.
     */

    //*********************************************************************************
    //**********************************************************************

    double Solve() override
    {
      KRATOS_TRY

        //calculate the BDF coefficients
        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
      double Dt = rCurrentProcessInfo[DELTA_TIME];
      int stationary= rCurrentProcessInfo[STATIONARY];

      unsigned int iter = 0;
      double ratio;
      bool is_converged = false;
      double dT_norm = 0.0;
      double T_norm = 0.0;

      if(stationary==1){
	iter = 0;
        //ratio;
        is_converged = false;
        dT_norm = 0.0;
        T_norm = 0.0;
	while (iter++ < mmax_iter && is_converged == false)
	  {
            rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
            dT_norm = mstep1->Solve();
            T_norm = CalculateTemperatureNorm();
            CalculateProjection();

            ratio = 1.00;
            if (T_norm != 0.00)
	      ratio = dT_norm / T_norm;
            else
	      {
                std::cout << "T_norm = " << T_norm << " dT_norm = " << dT_norm << std::endl;
	      }

            if (dT_norm < 1e-11)
	      ratio = 0; //converged

            if (ratio < mtoll)
	      is_converged = true;

            std::cout << " iter = " << iter << " ratio = " << ratio << std::endl;
	  }
      }
      else{

        if (mOldDt == 0.00) //needed for the first step
	  mOldDt = Dt;
        if (mtime_order == 2)
	  {
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

	iter = 0;
        //ratio;
        is_converged = false;
        dT_norm = 0.0;
        T_norm = 0.0;

        while (iter++ < mmax_iter && is_converged == false)
	  {
            rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
            dT_norm = mstep1->Solve();
            T_norm = CalculateTemperatureNorm();
            CalculateProjection();

            ratio = 1.00;
            if (T_norm != 0.00)
	      ratio = dT_norm / T_norm;
            else
            {
	      std::cout << "T_norm = " << T_norm << " dT_norm = " << dT_norm << std::endl;
            }

            if (dT_norm < 1e-11)
	      ratio = 0; //converged

            if (ratio < mtoll)
	      is_converged = true;

            std::cout << " iter = " << iter << " ratio = " << ratio << std::endl;
	  }
      }

      return dT_norm;
      KRATOS_CATCH("")
	}
    //******************************************************************************************************
    //******************************************************************************************************
    ///calculation of temperature norm

    double CalculateTemperatureNorm()
    {
        KRATOS_TRY;

        double norm = 0.00;

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();


        for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                i != BaseType::GetModelPart().NodesEnd(); ++i)
        {
            norm += pow(i->FastGetSolutionStepValue(rUnknownVar), 2);
        }

        return sqrt(norm);

        KRATOS_CATCH("")
    }

    //******************************************************************************************************
    //******************************************************************************************************
    ///calculation of projection

    void CalculateProjection()
    {
        KRATOS_TRY;

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

        const Variable<double>& rProjectionVariable = my_settings->GetProjectionVariable();

        //first of all set to zero the nodal variables to be updated nodally
        for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                i != BaseType::GetModelPart().NodesEnd(); ++i)
        {
            if ((i->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
            {


                (i)->FastGetSolutionStepValue(rProjectionVariable) = 0.00;
                (i)->FastGetSolutionStepValue(NODAL_AREA) = 0.00;
            }
            else
            {
                (i)->FastGetSolutionStepValue(NODAL_AREA) = 1.00;
            }
        }

        //add the elemental contributions for the calculation of the velocity
        //and the determination of the nodal area
        rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
        for (ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin();
                i != BaseType::GetModelPart().ElementsEnd(); ++i)
        {
            (i)->InitializeSolutionStep(rCurrentProcessInfo);
        }

        //solve nodally for the velocity
        for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                i != BaseType::GetModelPart().NodesEnd(); ++i)
        {
            double& conv_proj = (i)->FastGetSolutionStepValue(rProjectionVariable);
            double temp = 1.00 / (i)->FastGetSolutionStepValue(NODAL_AREA);
            conv_proj *= temp;
        }

        KRATOS_CATCH("")
    }

    void SetEchoLevel(int Level) override
    {
        mstep1->SetEchoLevel(Level);
    }

    void Clear() override
    {
        mstep1->Clear();
    }

    int Check() override
    {
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
        const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
        const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();
        const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
        const Variable<double>& rProjectionVariable = my_settings->GetProjectionVariable();


        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(rUnknownVar) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----UnknownVar---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(rDensityVar) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----rDensityVar---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(rDiffusionVar) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----rDiffusionVar---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(rSourceVar) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----rSourceVar---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(rMeshVelocityVar) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----rMeshVelocityVar---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(rProjectionVariable) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----rProjectionVariable---- variable!!!!!! ERROR", "");


        //if(BaseType::GetModelPart().GetMesh().GetProperties(0)[EMISSIVITY] == false)
        //    KRATOS_THROW_ERROR(std::logic_error, "Add  ----EMISSIVITY---- variable!!!!!! ERROR", "");
        //if(BaseType::GetModelPart().GetMesh().GetProperties(0)[CONVECTION_COEFFICIENT] == false)
        //    KRATOS_THROW_ERROR(std::logic_error, "Add  ----CONVECTION_COEFFICIENT---- variable!!!!!! ERROR", "");
        //if(BaseType::GetModelPart().GetMesh().GetProperties(0)[AMBIENT_TEMPERATURE] == false)
        //    KRATOS_THROW_ERROR(std::logic_error, "Add  ----AMBIENT_TEMPERATURE---- variable!!!!!! ERROR", "");


        return 0;

        KRATOS_CATCH("")

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
    double mOldDt;
    int mtime_order;
    unsigned int mmax_iter;
    double mtoll;





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
    ResidualBasedConvectionDiffusionStrategyNonLinear(const ResidualBasedConvectionDiffusionStrategyNonLinear& Other);


    /*@} */

}; /* Class ResidualBasedConvectionDiffusionStrategyNonLinear */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY_NONLINEAR  defined */

