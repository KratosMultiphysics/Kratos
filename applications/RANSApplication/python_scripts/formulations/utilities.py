import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS

from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.RANSApplication.model_part_factory import CreateDuplicateModelPart

if (IsDistributedRun()
        and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory as linear_solver_factory
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericResidualBasedSimpleSteadyScalarScheme as steady_scalar_scheme
    # from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIResidualBasedSimpleSteadyVelocityScheme as steady_velocity_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericResidualBasedBossakVelocityDynamicScalarScheme as bossak_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericScalarConvergenceCriteria as scalar_convergence_criteria
    from KratosMultiphysics.TrilinosApplication import TrilinosResidualCriteria as residual_criteria
    from KratosMultiphysics.TrilinosApplication import TrilinosNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosPeriodicBlockBuilderAndSolver as periodic_block_builder_and_solver
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosBlockBuilderAndSolver as block_builder_and_solver
    from KratosMultiphysics.TrilinosApplication import TrilinosResidualBasedIncrementalUpdateStaticScheme as incemental_update_static_scheme
elif (not IsDistributedRun()):
    from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
    # from KratosMultiphysics.RANSApplication import AlgebraicFluxCorrectedScalarSteadyScheme as steady_scalar_scheme
    # from KratosMultiphysics.RANSApplication import ResidualBasedSimpleSteadyVelocityScheme as steady_velocity_scheme
    from KratosMultiphysics.RANSApplication import GenericResidualBasedBossakVelocityDynamicScalarScheme as bossak_scheme
    from KratosMultiphysics.RANSApplication import GenericScalarConvergenceCriteria as scalar_convergence_criteria
    from Kratos import ResidualCriteria as residual_criteria
    from Kratos import ResidualBasedNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import PeriodicBlockBuilderAndSolver as periodic_block_builder_and_solver
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import BlockBuilderAndSolver as block_builder_and_solver
    from KratosMultiphysics import ResidualBasedIncrementalUpdateStaticScheme as incemental_update_static_scheme
else:
    raise Exception("Distributed run requires TrilinosApplication")


def CreateResidualBasedBlockBuilderAndSolver(linear_solver, is_periodic,
                                             communicator):
    if (is_periodic):
        return periodic_block_builder_and_solver(linear_solver, communicator)
    else:
        return block_builder_and_solver(linear_solver, communicator)


def CreateResidualBasedNewtonRaphsonStrategy(model_part, scheme, linear_solver,
                                             convergence_criteria,
                                             builder_and_solver,
                                             max_iterations, compute_reactions,
                                             reform_dofs_at_each_step,
                                             move_mesh_flag):
    return newton_raphson_strategy(model_part, scheme, linear_solver,
                                   convergence_criteria, builder_and_solver,
                                   max_iterations, compute_reactions,
                                   reform_dofs_at_each_step, move_mesh_flag)


def CreateLinearSolver(solver_settings):
    return linear_solver_factory.ConstructSolver(solver_settings)


def CreateResidualCriteria(relative_tolerance, absolute_tolerance):
    return residual_criteria(relative_tolerance, absolute_tolerance)


# def CreateScalarSteadyScheme(relaxation_factor):
#     return steady_scalar_scheme(relaxation_factor)


def CreateIncremantalUpdateScheme():
    return incemental_update_static_scheme()


def CreateScalarScheme(scheme_settings, relaxation_rate, scalar_variable,
                       scalar_variable_rate, relaxed_scalar_variable_rate):
    scheme_type = scheme_settings["scheme_type"].GetString()
    if (scheme_type == "bossak"):
        return bossak_scheme(
            scheme_settings["scheme_settings"]["alpha_bossak"].GetDouble(),
            relaxation_rate, scalar_variable, scalar_variable_rate,
            relaxed_scalar_variable_rate)
    elif (scheme_type == "steady"):
        return steady_scheme(relaxation_rate)
    else:
        raise Exception("Unsupported scheme type. [ scheme_type = \"" +
                        scheme_type + "\" ]")


def CalculateNormalsOnConditions(model_part):
    domain_size = model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
    if (not RansVariableUtilities.IsAnalysisStepCompleted(
            model_part, "CONDITION_NORMAL_CALCULATION")):
        Kratos.NormalCalculationUtils().CalculateOnSimplex(
            model_part.Conditions, domain_size)
        RansVariableUtilities.AddAnalysisStep(model_part,
                                              "CONDITION_NORMAL_CALCULATION")

        # RansVariableUtilities.AssignConditionVariableValuesToNodes(
        #     model_part, Kratos.NORMAL, Kratos.SLIP)

        Kratos.Logger.PrintInfo("NormalCalculationUtils",
                                "Condition normals calculated.")


def CreateFormulationModelPart(formulation, element_name, condition_name):
    formulation.domain_size = formulation.GetBaseModelPart().ProcessInfo[
        Kratos.DOMAIN_SIZE]

    element_suffix = str(
        formulation.domain_size) + "D" + str(formulation.domain_size + 1) + "N"
    condition_suffix = str(formulation.domain_size) + "D" + str(
        formulation.domain_size) + "N"

    element_name = element_name + element_suffix
    condition_name = condition_name + condition_suffix
    return CreateDuplicateModelPart(formulation.GetBaseModelPart(),
                                    formulation.GetName(), element_name,
                                    condition_name, "")
