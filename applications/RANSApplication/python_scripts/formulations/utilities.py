import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.RANSApplication.model_part_factory import CreateDuplicateModelPart

if (IsDistributedRun()
        and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory as linear_solver_factory
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericResidualBasedSimpleSteadyScalarScheme as steady_scalar_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIAlgebraicFluxCorrectedScalarSteadyScheme as afc_steady_scalar_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericResidualBasedBossakVelocityDynamicScalarScheme as bossak_scheme
    from KratosMultiphysics.TrilinosApplication import TrilinosResidualCriteria as residual_criteria
    from KratosMultiphysics.TrilinosApplication import TrilinosNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosPeriodicBlockBuilderAndSolver as periodic_block_builder_and_solver
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosBlockBuilderAndSolver as block_builder_and_solver
    from KratosMultiphysics.TrilinosApplication import TrilinosResidualBasedIncrementalUpdateStaticScheme as incemental_update_static_scheme
elif (not IsDistributedRun()):
    from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
    from KratosMultiphysics.RANSApplication import GenericResidualBasedSimpleSteadyScalarScheme as steady_scalar_scheme
    from KratosMultiphysics.RANSApplication import AlgebraicFluxCorrectedScalarSteadyScheme as afc_steady_scalar_scheme
    from KratosMultiphysics.RANSApplication import GenericResidualBasedBossakVelocityDynamicScalarScheme as bossak_scheme
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

def CreateVelocityPressureCriteria():
    pass

def CreateIncremantalUpdateScheme():
    return incemental_update_static_scheme()


def CreateSteadyScalarScheme(relaxation_factor):
    return steady_scalar_scheme(relaxation_factor)

def CreateSteadyAlgeraicFluxCorrectedTransportScheme(relaxation_factor):
    return afc_steady_scalar_scheme(relaxation_factor)

def CreateBossakScalarScheme(bossak_value, relaxation_factor, scalar_variable,
                             scalar_rate_variable,
                             relaxed_scalar_rate_variable):
    return bossak_scheme(bossak_value, relaxation_factor, scalar_variable,
                         scalar_rate_variable, relaxed_scalar_rate_variable)


def CalculateNormalsOnConditions(model_part):
    domain_size = model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
    if (not RansVariableUtilities.IsAnalysisStepCompleted(
            model_part, "CONDITION_NORMAL_CALCULATION")):
        Kratos.NormalCalculationUtils().CalculateOnSimplex(
            model_part.Conditions, domain_size)
        RansVariableUtilities.AddAnalysisStep(model_part,
                                              "CONDITION_NORMAL_CALCULATION")

        RansVariableUtilities.AssignConditionVariableValuesToNodes(
            model_part, Kratos.NORMAL, Kratos.SLIP)

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


def GetConvergenceInfo(variable,
                       relative_error,
                       relative_tolerance,
                       absolute_error=-1.0,
                       absolute_tolerance=-1.0):
    info = "[ Obtained ratio: {0:6e}; Expected ratio: {1:6e}".format(
        relative_error, relative_tolerance)
    if (absolute_error >= 0.0 and absolute_tolerance >= 0.0):
        info += "; Absolute norm: {0:6e}; Expected norm: {1:6e}".format(
            absolute_error, absolute_tolerance)
    info += " ] - " + str(variable.Name())
    return info

def IsBufferInitialized(formulation):
    return (formulation.GetBaseModelPart().ProcessInfo[Kratos.STEP] + 1 >= formulation.GetMinimumBufferSize())

def InitializePeriodicConditions(base_model_part, model_part, variables_list):
    properties = model_part.CreateNewProperties(
        model_part.NumberOfProperties() + 1)
    pcu = KratosCFD.PeriodicConditionUtilities(
        model_part, model_part.ProcessInfo[Kratos.DOMAIN_SIZE])
    for variable in variables_list:
        pcu.AddPeriodicVariable(properties, variable)

    index = model_part.NumberOfConditions()
    for condition in base_model_part.Conditions:
        if condition.Is(Kratos.PERIODIC):
            index += 1
            node_id_list = [node.Id for node in condition.GetNodes()]
            periodic_condition = model_part.CreateNewCondition(
                "PeriodicCondition", index, node_id_list, properties)
            periodic_condition.Set(Kratos.PERIODIC)