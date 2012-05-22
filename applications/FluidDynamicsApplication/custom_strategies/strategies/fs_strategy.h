#ifndef KRATOS_FS_STRATEGY_H
#define KRATOS_FS_STRATEGY_H

#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_elements/fractional_step.h"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

#include "custom_utilities/solver_settings.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{


///@}
///@name Type Definitions
///@{

///@}


///@name  Enum's
///@{


///@}
///@name  Functions
///@{



///@}
///@name Kratos Classes
///@{

template<class TSparseSpace,
class TDenseSpace,
class TLinearSolver
>
class FSStrategy : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of FSStrategy
    typedef boost::shared_ptr< FSStrategy<TSparseSpace, TDenseSpace, TLinearSolver> > Pointer;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

    typedef SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver> SolverSettingsType;

    ///@}
    ///@name Life Cycle
    ///@{

    FSStrategy(ModelPart& rModelPart,
               SolverSettingsType& rSolverConfig,
               bool PredictorCorrector):
        BaseType(rModelPart,false)
    {
        KRATOS_TRY;

        mDomainSize = rSolverConfig.GetDomainSize();
        mTimeOrder = rSolverConfig.GetTimeOrder();

        mPredictorCorrector = PredictorCorrector;

        mUseSlipConditions = rSolverConfig.UseSlipConditions();

        BaseType::SetEchoLevel(rSolverConfig.GetEchoLevel());

        // Initialize strategies for each step
        bool HaveVelStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Velocity,mpMomentumStrategy);

        if (HaveVelStrategy)
        {
            rSolverConfig.FindTolerance(SolverSettingsType::Velocity,mVelocityTolerance);
            rSolverConfig.FindMaxIter(SolverSettingsType::Velocity,mMaxVelocityIter);
        }
        else
        {
            KRATOS_ERROR(std::runtime_error,"FS_Strategy error: No Velocity strategy defined in FractionalStepSettings","");
        }

        bool HavePressStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Pressure,mpPressureStrategy);

        if (HavePressStrategy)
        {
            rSolverConfig.FindTolerance(SolverSettingsType::Pressure,mPressureTolerance);
            rSolverConfig.FindMaxIter(SolverSettingsType::Pressure,mMaxPressueIter);
        }
        else
        {
            KRATOS_ERROR(std::runtime_error,"FS_Strategy error: No Pressure strategy defined in FractionalStepSettings","");
        }

        Process::Pointer pTurbulenceProcess;
        bool HaveTurbulence = rSolverConfig.GetTurbulenceModel(pTurbulenceProcess);

        if (HaveTurbulence)
            mExtraIterationSteps.push_back(pTurbulenceProcess);

        // Set up nodes to use slip conditions if needed.
        if (mUseSlipConditions)
        {
#pragma omp parallel
            {
                ModelPart::ConditionIterator CondBegin;
                ModelPart::ConditionIterator CondEnd;
                OpenMPUtils::PartitionedIterators(rModelPart.Conditions(),CondBegin,CondEnd);

                for (ModelPart::ConditionIterator itCond = CondBegin; itCond != CondEnd; ++itCond)
                {
                    const double FlagValue = itCond->GetValue(IS_STRUCTURE);
                    if (FlagValue != 0.0)
                    {

                        Condition::GeometryType& rGeom = itCond->GetGeometry();
                        for (unsigned int i = 0; i < rGeom.PointsNumber(); ++i)
                        {
                            rGeom[i].SetLock();
                            rGeom[i].SetValue(IS_STRUCTURE,FlagValue);
                            rGeom[i].UnSetLock();
                        }
                    }
                }
            }
            rModelPart.GetCommunicator().AssembleNonHistoricalData(IS_STRUCTURE);
        }

        KRATOS_CATCH("");
    }

    FSStrategy(ModelPart& rModelPart,
               /*SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& rSolverConfig,*/
               typename TLinearSolver::Pointer pVelocityLinearSolver,
               typename TLinearSolver::Pointer pPressureLinearSolver,
               bool MoveMeshFlag, ///@todo: Read from solver configuration? Should match the one passed to vel/pre strategies?
               bool ReformDofAtEachIteration = true,
               double VelTol = 0.01,
               double PresTol = 0.01,
               int MaxVelocityIterations = 3,
               int MaxPressureIterations = 1,// Only for predictor-corrector
               unsigned int TimeOrder = 2, ///@todo check if really needed
               unsigned int DomainSize = 2,
               bool PredictorCorrector= true):
        BaseType(rModelPart,MoveMeshFlag), // Move Mesh flag, pass as input?
        mVelocityTolerance(VelTol),
        mPressureTolerance(PresTol),
        mMaxVelocityIter(MaxVelocityIterations),
        mMaxPressueIter(MaxPressureIterations),
        mDomainSize(DomainSize),
        mTimeOrder(TimeOrder),
        mPredictorCorrector(PredictorCorrector),
        mUseSlipConditions(true), ///@todo initialize somehow
        mExtraIterationSteps()
    {
        KRATOS_TRY;

        BaseType::SetEchoLevel(1);

        // Check that input parameters are reasonable and sufficient.
        this->Check();

//        // Get solution strategies from solver setup
//        this->mpMomentumStrategy = rSolverConfig.pGetStrategy(std::string("vel_strategy"));
//        this->mpPressureStrategy = rSolverConfig.pGetStrategy(std::string("pressure_strategy"));

        bool CalculateReactions = false;
        bool CalculateNormDxFlag = true;

        //computation of the fractional vel velocity (first step)
        //3 dimensional case
        typedef typename Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3 > > > VarComponent;
        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
        typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

        //initializing fractional velocity solution step
        typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
        typename SchemeType::Pointer pScheme;
        if (mUseSlipConditions)
        {
            typename SchemeType::Pointer Temp = typename SchemeType::Pointer(new ResidualBasedIncrementalUpdateStaticSchemeSlip< TSparseSpace, TDenseSpace > (mDomainSize,mDomainSize));
            pScheme.swap(Temp);
        }
        else
        {
            typename SchemeType::Pointer Temp = typename SchemeType::Pointer(new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());
            pScheme.swap(Temp);
        }

        //CONSTRUCTION OF VELOCITY
        BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (pVelocityLinearSolver));
//        BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverSlip<TSparseSpace, TDenseSpace, TLinearSolver, VarComponent > (pNewVelocityLinearSolver, this->mDomainSize, VELOCITY_X, VELOCITY_Y, VELOCITY_Z));
        this->mpMomentumStrategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (rModelPart, pScheme, pVelocityLinearSolver, vel_build, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
        this->mpMomentumStrategy->SetEchoLevel( BaseType::GetEchoLevel() );

        BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(
                    //new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(pPressureLinearSolver));
                new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pPressureLinearSolver, PRESSURE));

        this->mpPressureStrategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (rModelPart, pScheme, pPressureLinearSolver, pressure_build, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
        this->mpPressureStrategy->SetEchoLevel( BaseType::GetEchoLevel() );

        if (mUseSlipConditions)
        {
#pragma omp parallel
            {
                ModelPart::ConditionIterator CondBegin;
                ModelPart::ConditionIterator CondEnd;
                OpenMPUtils::PartitionedIterators(rModelPart.Conditions(),CondBegin,CondEnd);

                for (ModelPart::ConditionIterator itCond = CondBegin; itCond != CondEnd; ++itCond)
                {
                    const double FlagValue = itCond->GetValue(IS_STRUCTURE);
                    if (FlagValue != 0.0)
                    {

                        Condition::GeometryType& rGeom = itCond->GetGeometry();
                        for (unsigned int i = 0; i < rGeom.PointsNumber(); ++i)
                        {
                            rGeom[i].SetLock();
                            rGeom[i].SetValue(IS_STRUCTURE,FlagValue);
                            rGeom[i].UnSetLock();
                        }
                    }
                }
            }
            rModelPart.GetCommunicator().AssembleNonHistoricalData(IS_STRUCTURE);
        }


        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~FSStrategy(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual int Check()
    {
        KRATOS_TRY;

        // Check elements and conditions in the model part
        int ierr = BaseType::Check();
        if (ierr != 0) return ierr;

        ModelPart& rModelPart = BaseType::GetModelPart();

        if ( mTimeOrder == 2 && rModelPart.GetBufferSize() < 3 )
            KRATOS_ERROR(std::logic_error,"Buffer size too small for fractional step strategy (BDF2), needed 3, got ",rModelPart.GetBufferSize());
        if ( mTimeOrder == 1 && rModelPart.GetBufferSize() < 2 )
            KRATOS_ERROR(std::logic_error,"Buffer size too small for fractional step strategy (Backward Euler), needed 2, got ",rModelPart.GetBufferSize());

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        for ( ModelPart::ElementIterator itEl = rModelPart.ElementsBegin(); itEl != rModelPart.ElementsEnd(); ++itEl )
        {
            ierr = itEl->Check(rCurrentProcessInfo);
            if (ierr != 0) break;
        }

        for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond)
        {
            ierr = itCond->Check(rCurrentProcessInfo);
            if (ierr != 0) break;
        }

        return ierr;

        KRATOS_CATCH("");
    }

    virtual double Solve()
    {
        // Initialize BDF2 coefficients
        ModelPart& rModelPart = BaseType::GetModelPart();
        this->SetTimeCoefficients(rModelPart.GetProcessInfo());

        double NormDp = 0.0;

        if (mPredictorCorrector)
        {
            // Iterative solution for pressure
            for(unsigned int it = 0; it < mMaxPressueIter; ++it)
            {
                if ( BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
                    std::cout << "Pressure iteration " << it << std::endl;

                NormDp = this->SolveStep();

                bool Converged = this->CheckPressureConvergence(NormDp);

                if ( Converged ) break;
            }
        }
        else
        {
            // Solve for fractional step velocity, then update pressure once
            NormDp = this->SolveStep();
        }

        return NormDp;
    }


    virtual void CalculateReactions()
    {
        ModelPart& rModelPart = BaseType::GetModelPart();
        ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        // Set fractional step index to the momentum equation step
        int OriginalStep = rCurrentProcessInfo[FRACTIONAL_STEP];
        rCurrentProcessInfo.SetValue(FRACTIONAL_STEP,1);

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            const array_1d<double,3> Zero(3,0.0);

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
            {
                itNode->FastGetSolutionStepValue(REACTION) = Zero;
            }
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator ElemBegin;
            ModelPart::ElementIterator ElemEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(),ElemBegin,ElemEnd);

            LocalSystemVectorType RHS_Contribution;
            LocalSystemMatrixType LHS_Contribution;

            for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
            {

                //itElem->InitializeNonLinearIteration(rCurrentProcessInfo);

                // Build local system
                itElem->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);

                Element::GeometryType& rGeom = itElem->GetGeometry();
                unsigned int NumNodes = rGeom.PointsNumber();
                unsigned int index = 0;

                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    rGeom[i].SetLock();
                    array_1d<double,3>& rReaction = rGeom[i].FastGetSolutionStepValue(REACTION);
                    for (unsigned int d = 0; d < mDomainSize; ++d)
                        rReaction[d] -= RHS_Contribution[index++];
                    rGeom[i].UnSetLock();
                }
            }
        }

        rModelPart.GetCommunicator().AssembleCurrentData(REACTION);

        // Reset original fractional step index
        rCurrentProcessInfo.SetValue(FRACTIONAL_STEP,OriginalStep);
    }

    virtual void AddIterationStep(Process::Pointer pNewStep)
    {
        mExtraIterationSteps.push_back(pNewStep);
    }


    ///@}
    ///@name Access
    ///@{

    virtual void SetEchoLevel(int Level)
    {
        BaseType::SetEchoLevel(Level);
        mpMomentumStrategy->SetEchoLevel(Level);
        mpPressureStrategy->SetEchoLevel(Level);
    }


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "FSStrategy" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "FSStrategy";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:

    ///@name Protected Life Cycle
    ///@{

    /// Protected constructor to use in derived classes, does not initialize solution strategy-related members.
    FSStrategy(ModelPart& rModelPart,
               bool MoveMeshFlag, ///@todo: Read from solver configuration? Should match the one passed to vel/pre strategies?
               double VelTol,
               double PresTol,
               int MaxVelocityIterations,
               int MaxPressureIterations,// Only for predictor-corrector
               unsigned int TimeOrder, ///@todo check if really needed
               unsigned int DomainSize,
               bool PredictorCorrector,
               bool UseSlipConditions):
        BaseType(rModelPart,MoveMeshFlag), // Move Mesh flag, pass as input?
        mVelocityTolerance(VelTol),
        mPressureTolerance(PresTol),
        mMaxVelocityIter(MaxVelocityIterations),
        mMaxPressueIter(MaxPressureIterations),
        mDomainSize(DomainSize),
        mTimeOrder(TimeOrder),
        mPredictorCorrector(PredictorCorrector),
        mUseSlipConditions(UseSlipConditions),
        mExtraIterationSteps()
    {}

    ///@}
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    /// Calculate the coefficients for time iteration.
    /**
     * @param rCurrentProcessInfo ProcessInfo instance from the fluid ModelPart. Must contain DELTA_TIME and BDF_COEFFICIENTS variables.
     */
    void SetTimeCoefficients(ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if (mTimeOrder == 2)
        {
            //calculate the BDF coefficients
            double Dt = rCurrentProcessInfo[DELTA_TIME];
            double OldDt = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

            double Rho = OldDt / Dt;
            double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs.resize(3, false);

            BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho); //coefficient for step n+1 (3/2Dt if Dt is constant)
            BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
            BDFcoeffs[2] = TimeCoeff; //coefficient for step n-1 (1/2Dt if Dt is constant)
        }
        else if (mTimeOrder == 1)
        {
            double Dt = rCurrentProcessInfo[DELTA_TIME];
            double TimeCoeff = 1.0 / Dt;

            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs.resize(2, false);

            BDFcoeffs[0] = TimeCoeff; //coefficient for step n+1 (1/Dt)
            BDFcoeffs[1] = -TimeCoeff; //coefficient for step n (-1/Dt)
        }

        KRATOS_CATCH("");
    }

    double SolveStep()
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

        // 1. Fractional step momentum iteration
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);

        for(unsigned int it = 0; it < mMaxVelocityIter; ++it)
        {
            if ( BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
                std::cout << "Momentum iteration " << it << std::endl;

            // build momentum system and solve for fractional step velocity increment
            rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);
            double NormDv = mpMomentumStrategy->Solve();

//            // Compute projections (for stabilization)
//            rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
//            this->ComputeSplitOssProjections();

            // Additional steps
            for (std::vector<Process::Pointer>::iterator iExtraSteps = mExtraIterationSteps.begin();
                 iExtraSteps != mExtraIterationSteps.end(); ++iExtraSteps)
                (*iExtraSteps)->Execute();

            // Check convergence
            bool Converged = this->CheckFractionalStepConvergence(NormDv);

            if (Converged) break;
        }

        // Compute projections (for stabilization)
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
        this->ComputeSplitOssProjections();

        // 2. Pressure solution (store pressure variation in PRESSURE_OLD_IT)
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,5);

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
            {
                const double OldPress = itNode->FastGetSolutionStepValue(PRESSURE);
                itNode->FastGetSolutionStepValue(PRESSURE_OLD_IT) = -OldPress;
            }
        }

        double NormDp = mpPressureStrategy->Solve();

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
                itNode->FastGetSolutionStepValue(PRESSURE_OLD_IT) += itNode->FastGetSolutionStepValue(PRESSURE);
        }

        // 3. Compute end-of-step velocity
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,6);

        this->CalculateEndOfStepVelocity();


        return NormDp;
    }

    bool CheckFractionalStepConvergence(const double NormDv)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

        double NormV = 0.00;

#pragma omp parallel reduction(+:NormV)
        {
            ModelPart::NodeIterator NodeBegin;
            ModelPart::NodeIterator NodeEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodeBegin,NodeEnd);

            for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
            {
                const array_1d<double,3> &Vel = itNode->FastGetSolutionStepValue(VELOCITY);

                for (unsigned int d = 0; d < 3; ++d)
                    NormV += Vel[d] * Vel[d];
            }
        }

        BaseType::GetModelPart().GetCommunicator().SumAll(NormV);

        NormV = sqrt(NormV);

        if (NormV == 0.0) NormV = 1.00;

        double Ratio = NormDv / NormV;

        if ( BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
            std::cout << "Relative error: " << Ratio << std::endl;

        if (Ratio < mVelocityTolerance)
        {
            return true;
        }
        else
            return false;
    }

    bool CheckPressureConvergence(const double NormDp)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

        double NormP = 0.00;

#pragma omp parallel reduction(+:NormP)
        {
            ModelPart::NodeIterator NodeBegin;
            ModelPart::NodeIterator NodeEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodeBegin,NodeEnd);

            for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
            {
                const double Pr = itNode->FastGetSolutionStepValue(PRESSURE);
                NormP += Pr * Pr;
            }
        }

        BaseType::GetModelPart().GetCommunicator().SumAll(NormP);

        NormP = sqrt(NormP);

        if (NormP == 0.0) NormP = 1.00;

        double Ratio = NormDp / NormP;

        if ( BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
            std::cout << "Relative error: " << Ratio << std::endl;

        if (Ratio < mPressureTolerance)
        {
            return true;
        }
        else
            return false;
    }


    void ComputeSplitOssProjections()
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

        const array_1d<double,3> Zero(3,0.0);

        array_1d<double,3> Out(3,0.0);

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for ( ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode )
            {
                itNode->FastGetSolutionStepValue(CONV_PROJ) = Zero;
                itNode->FastGetSolutionStepValue(PRESS_PROJ) = Zero;
                itNode->FastGetSolutionStepValue(DIVPROJ) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
            }
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator ElemBegin;
            ModelPart::ElementIterator ElemEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(),ElemBegin,ElemEnd);

            for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
            {
                itElem->Calculate(CONV_PROJ,Out,rModelPart.GetProcessInfo());
            }
        }

        rModelPart.GetCommunicator().AssembleCurrentData(CONV_PROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(PRESS_PROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);


#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for ( ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode )
            {
                const double NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);
                itNode->FastGetSolutionStepValue(CONV_PROJ) /= NodalArea;
                itNode->FastGetSolutionStepValue(PRESS_PROJ) /= NodalArea;
                itNode->FastGetSolutionStepValue(DIVPROJ) /= NodalArea;
            }
        }
    }

    void CalculateEndOfStepVelocity()
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

        const array_1d<double,3> Zero(3,0.0);
        array_1d<double,3> Out(3,0.0);

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for ( ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode )
            {
                itNode->FastGetSolutionStepValue(FRACT_VEL) = Zero;
            }
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator ElemBegin;
            ModelPart::ElementIterator ElemEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(),ElemBegin,ElemEnd);

            for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
            {
                itElem->Calculate(VELOCITY,Out,rModelPart.GetProcessInfo());
            }
        }

        rModelPart.GetCommunicator().AssembleCurrentData(FRACT_VEL);

        // Force the end of step velocity to verify slip conditions in the model
        if (mUseSlipConditions)
            this->EnforceSlipCondition(IS_STRUCTURE);

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for ( ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode )
            {
                const double NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);
                if ( ! itNode->IsFixed(VELOCITY_X) )
                    itNode->FastGetSolutionStepValue(VELOCITY_X) += itNode->FastGetSolutionStepValue(FRACT_VEL_X) / NodalArea;
                if ( ! itNode->IsFixed(VELOCITY_Y) )
                    itNode->FastGetSolutionStepValue(VELOCITY_Y) += itNode->FastGetSolutionStepValue(FRACT_VEL_Y) / NodalArea;
                if ( mDomainSize > 2 && ( ! itNode->IsFixed(VELOCITY_Z) ) )
                    itNode->FastGetSolutionStepValue(VELOCITY_Z) += itNode->FastGetSolutionStepValue(FRACT_VEL_Z) / NodalArea;
            }
        }
    }

    /**
     * @brief Substract wall-normal component of velocity update to ensure that the final velocity satisfies slip conditions.
     * @param rSlipWallFlag If Node.GetValue(rSlipWallFlag) != 0, the node is in the wall.
     */
    void EnforceSlipCondition(Variable<double>& rSlipWallFlag)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

#pragma omp parallel
        {
            ModelPart::NodeIterator NodeBegin; // = rModelPart.NodesBegin();
            ModelPart::NodeIterator NodeEnd; // = rModelPart.NodesEnd();
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodeBegin,NodeEnd);

            for ( ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode )
            {
                if ( itNode->GetValue(rSlipWallFlag) != 0.0 )
                {
                    const array_1d<double,3>& rNormal = itNode->FastGetSolutionStepValue(NORMAL);
                    array_1d<double,3>& rDeltaVelocity = itNode->FastGetSolutionStepValue(FRACT_VEL);

                    double Proj = rNormal[0] * rDeltaVelocity[0];
                    double Norm = rNormal[0] * rNormal[0];

                    for (unsigned int d = 1; d < mDomainSize; ++d)
                    {
                        Proj += rNormal[d] * rDeltaVelocity[d];
                        Norm += rNormal[d] * rNormal[d];
                    }

                    Proj /= Norm;
                    rDeltaVelocity -= Proj * rNormal;
                }
            }
        }
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    double mVelocityTolerance;

    double mPressureTolerance;

    unsigned int mMaxVelocityIter;

    unsigned int mMaxPressueIter;

    unsigned int mDomainSize;

    unsigned int mTimeOrder;

    bool mPredictorCorrector;

    bool mUseSlipConditions;

    // Fractional step index.
    /*  1 : Momentum step (calculate fractional step velocity)
      * 2-3 : Unused (reserved for componentwise calculation of frac step velocity)
      * 4 : Pressure step
      * 5 : Computation of projections
      * 6 : End of step velocity
      */
//    unsigned int mStepId;

    /// Scheme for the solution of the momentum equation
    StrategyPointerType mpMomentumStrategy;

    /// Scheme for the solution of the mass equation
    StrategyPointerType mpPressureStrategy;

    std::vector< Process::Pointer > mExtraIterationSteps;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FSStrategy& operator=(FSStrategy const& rOther){}

    /// Copy constructor.
    FSStrategy(FSStrategy const& rOther){}


    ///@}

}; /// Class FStepStrategy

///@}
///@name Type Definitions
///@{


///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_FS_STRATEGY_H
