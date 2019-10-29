
#include "fractional_step_settings_for_chimera.h"

namespace Kratos
{

template< class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
FractionalStepSettingsForChimera<TSparseSpace, TDenseSpace, TLinearSolver>::FractionalStepSettingsForChimera(ModelPart &rModelPart,
                                                                   const std::size_t ThisDomainSize,
                                                                   const std::size_t ThisTimeOrder,
                                                                   const bool UseSlip,
                                                                   const bool MoveMeshFlag,
                                                                   const bool ReformDofSet)
    : BaseType(rModelPart, ThisDomainSize, ThisTimeOrder, UseSlip, MoveMeshFlag, ReformDofSet)
{
}

template< class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
FractionalStepSettingsForChimera<TSparseSpace, TDenseSpace, TLinearSolver>::~FractionalStepSettingsForChimera() {}

template< class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
void FractionalStepSettingsForChimera<TSparseSpace, TDenseSpace, TLinearSolver>::SetStrategy(StrategyLabel const &rStrategyLabel,
                                                   typename TLinearSolver::Pointer pLinearSolver,
                                                   const double Tolerance,
                                                   const unsigned int MaxIter)
{
    KRATOS_TRY;

    // pointer types for solution strategy construcion
    typedef typename Scheme<TSparseSpace, TDenseSpace>::Pointer SchemePointerType;
    typedef ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera<TSparseSpace, TDenseSpace, TLinearSolver> ResidualBasedBlockBuilderAndSolverWithConstraintsForChimeraType;

    // Default, fixed flags
    bool CalculateReactions = false;
    bool CalculateNormDxFlag = true;

    ModelPart &rModelPart = BaseType::GetModelPart();

    bool UseSlip = BaseType::UseSlipConditions();
    // Modification of the DofSet is managed by the fractional step strategy, not the auxiliary velocity and pressure strategies.
    bool ReformDofSet = false; //BaseType::GetReformDofSet();
    std::size_t EchoLevel = BaseType::GetEchoLevel();
    std::size_t StrategyEchoLevel = (EchoLevel > 0) ? (EchoLevel - 1) : 0;

    if (rStrategyLabel == BaseType::Velocity)
    {
        ModelPart &r_fs_velocity_model_part = rModelPart.CreateSubModelPart("fs_velocity_model_part");
        FastTransferBetweenModelPartsProcess(r_fs_velocity_model_part, rModelPart).Execute();
        // Velocity Builder and Solver
        //ResidualBasedBlockBuilderAndSolverPointerType pBuildAndSolver = new ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera<TSparseSpace, TDenseSpace, TLinearSolver >(pLinearSolver);
        ResidualBasedBlockBuilderAndSolverPointerType pBuildAndSolver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolverWithConstraintsForChimeraType>(pLinearSolver);

        SchemePointerType pScheme;
        //initializing fractional velocity solution step
        if (UseSlip)
        {
            std::size_t DomainSize = BaseType::GetDomainSize();
            SchemePointerType Temp = SchemePointerType(new ResidualBasedIncrementalUpdateStaticSchemeSlip<TSparseSpace, TDenseSpace>(DomainSize, DomainSize));
            pScheme.swap(Temp);
        }
        else
        {
            SchemePointerType Temp = SchemePointerType(new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>());
            pScheme.swap(Temp);
        }

        // Strategy
        BaseType::mStrategies[rStrategyLabel] = StrategyPointerType(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(r_fs_velocity_model_part, pScheme, pLinearSolver, pBuildAndSolver, CalculateReactions, ReformDofSet, CalculateNormDxFlag));
    }
    else if (rStrategyLabel == BaseType::Pressure)
    {
        ModelPart &r_fs_pressure_model_part = rModelPart.CreateSubModelPart("fs_pressure_model_part");
        FastTransferBetweenModelPartsProcess(r_fs_pressure_model_part, rModelPart).Execute();
        // Pressure Builder and Solver
        ResidualBasedBlockBuilderAndSolverPointerType pBuildAndSolver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolverWithConstraintsForChimeraType>(pLinearSolver);
        SchemePointerType pScheme = SchemePointerType(new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>());

        // Strategy
        BaseType::mStrategies[rStrategyLabel] = StrategyPointerType(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(r_fs_pressure_model_part, pScheme, pLinearSolver, pBuildAndSolver, CalculateReactions, ReformDofSet, CalculateNormDxFlag));
    }
    else
    {
        KRATOS_THROW_ERROR(std::runtime_error, "Error in FractionalStepSettingsForChimera: Unknown strategy label.", "");
    }

    BaseType::mTolerances[rStrategyLabel] = Tolerance;

    BaseType::mMaxIter[rStrategyLabel] = MaxIter;

    BaseType::mStrategies[rStrategyLabel]->SetEchoLevel(StrategyEchoLevel);

    KRATOS_CATCH("");
}

template< class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
bool FractionalStepSettingsForChimera<TSparseSpace, TDenseSpace, TLinearSolver>::FindStrategy(StrategyLabel const &rStrategyLabel,
                                                    StrategyPointerType &pThisStrategy)
{
    pThisStrategy = BaseType::mStrategies[rStrategyLabel];
    if (pThisStrategy != nullptr)
        return true;
    return false;
}

template< class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
std::string FractionalStepSettingsForChimera<TSparseSpace, TDenseSpace, TLinearSolver>::Info() const
{
    std::stringstream buffer;
    buffer << "FractionalStepSettingsForChimera";
    return buffer.str();
}

/// Print information about this object.
template< class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
void FractionalStepSettingsForChimera<TSparseSpace, TDenseSpace, TLinearSolver>::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "FractionalStepSettingsForChimera";
}

/// Print object's data.
template< class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
void FractionalStepSettingsForChimera<TSparseSpace, TDenseSpace, TLinearSolver>::PrintData(std::ostream &rOStream) const {}

} // namespace Kratos.
