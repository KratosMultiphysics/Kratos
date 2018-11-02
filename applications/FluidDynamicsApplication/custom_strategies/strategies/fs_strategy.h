//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#ifndef KRATOS_FS_STRATEGY_H
#define KRATOS_FS_STRATEGY_H

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"
#include "utilities/openmp_utils.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"

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
    KRATOS_CLASS_POINTER_DEFINITION(FSStrategy);

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
        BaseType(rModelPart,false),
        mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {
        InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    FSStrategy(ModelPart& rModelPart,
               SolverSettingsType& rSolverConfig,
               bool PredictorCorrector,
               const Kratos::Variable<int>& PeriodicVar):
        BaseType(rModelPart,false),
        mrPeriodicIdVar(PeriodicVar)
    {
        InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    FSStrategy(ModelPart& rModelPart,
               /*SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& rSolverConfig,*/
               typename TLinearSolver::Pointer pVelocityLinearSolver,
               typename TLinearSolver::Pointer pPressureLinearSolver,
               bool MoveMeshFlag, ///@todo: Read from solver configuration? Should match the one passed to vel/pre strategies?
               bool ReformDofSet = true,
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
        mMaxPressureIter(MaxPressureIterations),
        mDomainSize(DomainSize),
        mTimeOrder(TimeOrder),
        mPredictorCorrector(PredictorCorrector),
        mUseSlipConditions(true), ///@todo initialize somehow
        mReformDofSet(ReformDofSet),
        mExtraIterationSteps(),
        mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {
        KRATOS_TRY;

        BaseType::SetEchoLevel(1);

        // Check that input parameters are reasonable and sufficient.
        this->Check();

        bool CalculateReactions = false;
        bool CalculateNormDxFlag = true;

        bool ReformDofAtEachIteration = false; // DofSet modifiaction is managed by the fractional step strategy, auxiliary strategies should not modify the DofSet directly.

        // Additional Typedefs
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
        this->mpMomentumStrategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (rModelPart, pScheme, pVelocityLinearSolver, vel_build, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
        this->mpMomentumStrategy->SetEchoLevel( BaseType::GetEchoLevel() );

        BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(
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
    ~FSStrategy() override{}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {}

    int Check() override
    {
        KRATOS_TRY;

        // Check elements and conditions in the model part
        int ierr = BaseType::Check();
        if (ierr != 0) return ierr;

        if(DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::runtime_error,"DELTA_TIME Key is 0. Check that the application was correctly registered.","");
        if(BDF_COEFFICIENTS.Key() == 0)
            KRATOS_THROW_ERROR(std::runtime_error,"BDF_COEFFICIENTS Key is 0. Check that the application was correctly registered.","");

        ModelPart& rModelPart = BaseType::GetModelPart();

        if ( mTimeOrder == 2 && rModelPart.GetBufferSize() < 3 )
            KRATOS_THROW_ERROR(std::invalid_argument,"Buffer size too small for fractional step strategy (BDF2), needed 3, got ",rModelPart.GetBufferSize());
        if ( mTimeOrder == 1 && rModelPart.GetBufferSize() < 2 )
            KRATOS_THROW_ERROR(std::invalid_argument,"Buffer size too small for fractional step strategy (Backward Euler), needed 2, got ",rModelPart.GetBufferSize());

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

    double Solve() override
    {
        // Initialize BDF2 coefficients
        ModelPart& rModelPart = BaseType::GetModelPart();
        this->SetTimeCoefficients(rModelPart.GetProcessInfo());

        double NormDp = 0.0;

        if (mPredictorCorrector)
        {
            bool Converged = false;

            // Iterative solution for pressure
            for(unsigned int it = 0; it < mMaxPressureIter; ++it)
            {
                if ( BaseType::GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)
                    std::cout << "Pressure iteration " << it << std::endl;

                NormDp = this->SolveStep();

                Converged = this->CheckPressureConvergence(NormDp);

                if ( Converged )
                {
                    if ( BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
                        std::cout << "Predictor-corrector converged in " << it+1 << " iterations." << std::endl;
                    break;
                }
            }
            if (!Converged && BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
                std::cout << "Predictor-correctior iterations did not converge." << std::endl;

        }
        else
        {
            // Solve for fractional step velocity, then update pressure once
            NormDp = this->SolveStep();
        }

        if (mReformDofSet)
            this->Clear();



        return NormDp;
    }

    bool SolveSolutionStep() override
    {
        double norm_dp = this->Solve();
        /* If not doing predictor corrector iterations, norm_dp will
         * typically be "large" since we are not iterating on pressure.
         * It makes no sense to report that the iteration didn't converge
         * based on this.
         */
        return mPredictorCorrector ? this->CheckPressureConvergence(norm_dp) : true;
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

    virtual void ClearExtraIterationSteps()
    {
        mExtraIterationSteps.clear();
    }

    void Clear() override
    {
        mpMomentumStrategy->Clear();
        mpPressureStrategy->Clear();
    }


    ///@}
    ///@name Access
    ///@{

    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
        int StrategyLevel = Level > 0 ? Level - 1 : 0;
        mpMomentumStrategy->SetEchoLevel(StrategyLevel);
        mpPressureStrategy->SetEchoLevel(StrategyLevel);
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

        bool Converged = false;
        int Rank = rModelPart.GetCommunicator().MyPID();

        for(unsigned int it = 0; it < mMaxVelocityIter; ++it)
        {
            if ( BaseType::GetEchoLevel() > 1 && Rank == 0)
                std::cout << "Momentum iteration " << it << std::endl;

            // build momentum system and solve for fractional step velocity increment
            rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);
            double NormDv = mpMomentumStrategy->Solve();

//            // Compute projections (for stabilization)
//            rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
//            this->ComputeSplitOssProjections(rModelPart);

//            // Additional steps // Moved to end of step
//            for (std::vector<Process::Pointer>::iterator iExtraSteps = mExtraIterationSteps.begin();
//                 iExtraSteps != mExtraIterationSteps.end(); ++iExtraSteps)
//                (*iExtraSteps)->Execute();

            // Check convergence
            Converged = this->CheckFractionalStepConvergence(NormDv);

            if (Converged)
            {
                if ( BaseType::GetEchoLevel() > 0 && Rank == 0)
                    std::cout << "Fractional velocity converged in " << it+1 << " iterations." << std::endl;
                break;
            }
        }

        if (!Converged && BaseType::GetEchoLevel() > 0 && Rank == 0)
            std::cout << "Fractional velocity iterations did not converge." << std::endl;

        // Compute projections (for stabilization)
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
        this->ComputeSplitOssProjections(rModelPart);

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

        if (BaseType::GetEchoLevel() > 0 && Rank == 0)
            std::cout << "Calculating Pressure." << std::endl;
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
        if (BaseType::GetEchoLevel() > 0 && Rank == 0)
            std::cout << "Updating Velocity." << std::endl;
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,6);

        this->CalculateEndOfStepVelocity();

        /*
        mpPressureStrategy->Clear();
        double NormDu = mpPressureStrategy->Solve();
        mpPressureStrategy->Clear();
        */

        // Additional steps
        for (std::vector<Process::Pointer>::iterator iExtraSteps = mExtraIterationSteps.begin();
             iExtraSteps != mExtraIterationSteps.end(); ++iExtraSteps)
            (*iExtraSteps)->Execute();


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
            std::cout << "Fractional velocity relative error: " << Ratio << std::endl;

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
            std::cout << "Pressure relative error: " << Ratio << std::endl;

        if (Ratio < mPressureTolerance)
        {
            return true;
        }
        else
            return false;
    }


    void ComputeSplitOssProjections(ModelPart& rModelPart)
    {
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

        // If there are periodic conditions, add contributions from both sides to the periodic nodes
        //this->PeriodicConditionProjectionCorrection(rModelPart);

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
        //this->PeriodicConditionVelocityCorrection(rModelPart);

        // Force the end of step velocity to verify slip conditions in the model
        if (mUseSlipConditions)
            this->EnforceSlipCondition(IS_STRUCTURE);

        if (mDomainSize > 2)
        {
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
                    if ( ! itNode->IsFixed(VELOCITY_Z) )
                        itNode->FastGetSolutionStepValue(VELOCITY_Z) += itNode->FastGetSolutionStepValue(FRACT_VEL_Z) / NodalArea;
                }
            }
        }
        else
        {
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
                }
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

        const int num_nodes_in_model_part = rModelPart.NumberOfNodes();

        #pragma omp parallel for
        for (int i = 0; i < num_nodes_in_model_part; i++)
        {
            ModelPart::NodeIterator itNode = rModelPart.NodesBegin() + i;
            const Node<3>& r_const_node = *itNode;

            if ( r_const_node.GetValue(rSlipWallFlag) != 0.0 )
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

    /** On periodic boundaries, the nodal area and the values to project need to take into account contributions from elements on
     * both sides of the boundary. This is done using the conditions and the non-historical nodal data containers as follows:\n
     * 1- The partition that owns the PeriodicCondition adds the values on both nodes to their non-historical containers.\n
     * 2- The non-historical containers are added across processes, transmiting the right value from the condition owner to all partitions.\n
     * 3- The value on all periodic nodes is replaced by the one received in step 2.
     */
     void PeriodicConditionProjectionCorrection(ModelPart& rModelPart)
     {
         if (mrPeriodicIdVar.Key() != 0)
         {
             int GlobalNodesNum = rModelPart.GetCommunicator().LocalMesh().Nodes().size();
             rModelPart.GetCommunicator().SumAll(GlobalNodesNum);

             for (typename ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); itCond++ )
             {
                 ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
                 if (rGeom.PointsNumber() == 2)
                 {
                     Node<3>& rNode0 = rGeom[0];
                     int Node0Pair = rNode0.FastGetSolutionStepValue(mrPeriodicIdVar);

                     Node<3>& rNode1 = rGeom[1];
                     int Node1Pair = rNode1.FastGetSolutionStepValue(mrPeriodicIdVar);

                     // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                     if ( ( static_cast<int>(rNode0.Id()) == Node1Pair ) && (static_cast<int>(rNode1.Id()) == Node0Pair ) )
                     {
                         double NodalArea = rNode0.FastGetSolutionStepValue(NODAL_AREA) + rNode1.FastGetSolutionStepValue(NODAL_AREA);
                         array_1d<double,3> ConvProj = rNode0.FastGetSolutionStepValue(CONV_PROJ) + rNode1.FastGetSolutionStepValue(CONV_PROJ);
                         array_1d<double,3> PressProj = rNode0.FastGetSolutionStepValue(PRESS_PROJ) + rNode1.FastGetSolutionStepValue(PRESS_PROJ);
                         double DivProj = rNode0.FastGetSolutionStepValue(DIVPROJ) + rNode1.FastGetSolutionStepValue(DIVPROJ);

                         rNode0.GetValue(NODAL_AREA) = NodalArea;
                         rNode0.GetValue(CONV_PROJ) = ConvProj;
                         rNode0.GetValue(PRESS_PROJ) = PressProj;
                         rNode0.GetValue(DIVPROJ) = DivProj;
                         rNode1.GetValue(NODAL_AREA) = NodalArea;
                         rNode1.GetValue(CONV_PROJ) = ConvProj;
                         rNode1.GetValue(PRESS_PROJ) = PressProj;
                         rNode1.GetValue(DIVPROJ) = DivProj;
                     }
                 }
                 else if (rGeom.PointsNumber() == 4 && rGeom[0].FastGetSolutionStepValue(mrPeriodicIdVar) > GlobalNodesNum)
                 {
                     double NodalArea = rGeom[0].FastGetSolutionStepValue(NODAL_AREA);
                     array_1d<double,3> ConvProj = rGeom[0].FastGetSolutionStepValue(CONV_PROJ);
                     array_1d<double,3> PressProj = rGeom[0].FastGetSolutionStepValue(PRESS_PROJ);
                     double DivProj = rGeom[0].FastGetSolutionStepValue(DIVPROJ);

                     for (unsigned int i = 1; i < 4; i++)
                     {
                         NodalArea += rGeom[i].FastGetSolutionStepValue(NODAL_AREA);
                         ConvProj += rGeom[i].FastGetSolutionStepValue(CONV_PROJ);
                         PressProj += rGeom[i].FastGetSolutionStepValue(PRESS_PROJ);
                         DivProj += rGeom[i].FastGetSolutionStepValue(DIVPROJ);
                     }

                     for (unsigned int i = 0; i < 4; i++)
                     {
                         rGeom[i].GetValue(NODAL_AREA) = NodalArea;
                         rGeom[i].GetValue(CONV_PROJ) = ConvProj;
                         rGeom[i].GetValue(PRESS_PROJ) = PressProj;
                         rGeom[i].GetValue(DIVPROJ) = DivProj;
                     }
                 }
             }

             rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
             rModelPart.GetCommunicator().AssembleNonHistoricalData(CONV_PROJ);
             rModelPart.GetCommunicator().AssembleNonHistoricalData(PRESS_PROJ);
             rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);

             for (typename ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
             {
                 if (itNode->GetValue(NODAL_AREA) != 0.0)
                 {
                     itNode->FastGetSolutionStepValue(NODAL_AREA) = itNode->GetValue(NODAL_AREA);
                     itNode->FastGetSolutionStepValue(CONV_PROJ) = itNode->GetValue(CONV_PROJ);
                     itNode->FastGetSolutionStepValue(PRESS_PROJ) = itNode->GetValue(PRESS_PROJ);
                     itNode->FastGetSolutionStepValue(DIVPROJ) = itNode->GetValue(DIVPROJ);

                     // reset for next iteration
                     itNode->GetValue(NODAL_AREA) = 0.0;
                     itNode->GetValue(CONV_PROJ) = array_1d<double,3>(3,0.0);
                     itNode->GetValue(PRESS_PROJ) = array_1d<double,3>(3,0.0);
                     itNode->GetValue(DIVPROJ) = 0.0;
                 }
             }
         }
     }

     void PeriodicConditionVelocityCorrection(ModelPart& rModelPart)
     {
         if (mrPeriodicIdVar.Key() != 0)
         {
             int GlobalNodesNum = rModelPart.GetCommunicator().LocalMesh().Nodes().size();
             rModelPart.GetCommunicator().SumAll(GlobalNodesNum);

             for (typename ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); itCond++ )
             {
                 ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
                 if (rGeom.PointsNumber() == 2)
                 {
                     Node<3>& rNode0 = rGeom[0];
                     int Node0Pair = rNode0.FastGetSolutionStepValue(mrPeriodicIdVar);

                     Node<3>& rNode1 = rGeom[1];
                     int Node1Pair = rNode1.FastGetSolutionStepValue(mrPeriodicIdVar);

                     // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                     if ( ( static_cast<int>(rNode0.Id()) == Node1Pair ) && (static_cast<int>(rNode1.Id()) == Node0Pair ) )
                     {
                         array_1d<double,3> DeltaVel = rNode0.FastGetSolutionStepValue(FRACT_VEL) + rNode1.FastGetSolutionStepValue(FRACT_VEL);

                         rNode0.GetValue(FRACT_VEL) = DeltaVel;
                         rNode1.GetValue(FRACT_VEL) = DeltaVel;
                     }
                 }
                 else if (rGeom.PointsNumber() == 4 && rGeom[0].FastGetSolutionStepValue(mrPeriodicIdVar) > GlobalNodesNum)
                 {
                     array_1d<double,3> DeltaVel = rGeom[0].FastGetSolutionStepValue(FRACT_VEL);
                     for (unsigned int i = 1; i < 4; i++)
                     {
                         DeltaVel += rGeom[i].FastGetSolutionStepValue(FRACT_VEL);
                     }

                     for (unsigned int i = 0; i < 4; i++)
                     {
                         rGeom[i].GetValue(FRACT_VEL) = DeltaVel;
                     }
                 }
             }

             rModelPart.GetCommunicator().AssembleNonHistoricalData(FRACT_VEL);

             for (typename ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
             {
                 array_1d<double,3>& rDeltaVel = itNode->GetValue(FRACT_VEL);
                 if ( rDeltaVel[0]*rDeltaVel[0] + rDeltaVel[1]*rDeltaVel[1] + rDeltaVel[2]*rDeltaVel[2] != 0.0)
                 {
                     itNode->FastGetSolutionStepValue(FRACT_VEL) = itNode->GetValue(FRACT_VEL);
                     rDeltaVel = array_1d<double,3>(3,0.0);
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

    unsigned int mMaxPressureIter;

    unsigned int mDomainSize;

    unsigned int mTimeOrder;

    bool mPredictorCorrector;

    bool mUseSlipConditions;

    bool mReformDofSet;

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

    const Kratos::Variable<int>& mrPeriodicIdVar;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    void InitializeStrategy(SolverSettingsType& rSolverConfig,
            bool PredictorCorrector)
    {
        KRATOS_TRY;

        mTimeOrder = rSolverConfig.GetTimeOrder();

        // Check that input parameters are reasonable and sufficient.
        this->Check();

        ModelPart& rModelPart = this->GetModelPart();

        mDomainSize = rSolverConfig.GetDomainSize();

        mPredictorCorrector = PredictorCorrector;

        mUseSlipConditions = rSolverConfig.UseSlipConditions();

        mReformDofSet = rSolverConfig.GetReformDofSet();

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
            KRATOS_THROW_ERROR(std::runtime_error,"FS_Strategy error: No Velocity strategy defined in FractionalStepSettings","");
        }

        bool HavePressStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Pressure,mpPressureStrategy);

        if (HavePressStrategy)
        {
            rSolverConfig.FindTolerance(SolverSettingsType::Pressure,mPressureTolerance);
            rSolverConfig.FindMaxIter(SolverSettingsType::Pressure,mMaxPressureIter);
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"FS_Strategy error: No Pressure strategy defined in FractionalStepSettings","");
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
                    const Condition& rCond = *itCond;
                    const double& FlagValue = rCond.GetValue(IS_STRUCTURE);
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

        // Check input parameters
        this->Check();

        KRATOS_CATCH("");
    }


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
