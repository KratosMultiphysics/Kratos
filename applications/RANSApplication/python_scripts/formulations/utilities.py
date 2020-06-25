import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics import VariableUtils
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

from KratosMultiphysics.RANSApplication import RansVariableUtilities

if (IsDistributedRun()
        and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory as linear_solver_factory
    # from KratosMultiphysics.TrilinosApplication import TrilinosResidualCriteria as residual_criteria
    # the core residual criteria is not used since, the ratio is calcualted based on ratio between solution step start and end error,
    # for pseudo time stepping steady problems, this does not work
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericScalarConvergenceCriteria as residual_criteria
    from KratosMultiphysics.TrilinosApplication import TrilinosNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosPeriodicBlockBuilderAndSolver as periodic_block_builder_and_solver
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosBlockBuilderAndSolver as block_builder_and_solver
    from KratosMultiphysics.TrilinosApplication import TrilinosResidualBasedIncrementalUpdateStaticScheme as incremental_update_static_scheme
elif (not IsDistributedRun()):
    from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
    from KratosMultiphysics.RANSApplication import GenericScalarConvergenceCriteria as residual_criteria
    # from Kratos import ResidualCriteria as residual_criteria
    from Kratos import ResidualBasedNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import PeriodicBlockBuilderAndSolver as periodic_block_builder_and_solver
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import BlockBuilderAndSolver as block_builder_and_solver
    from KratosMultiphysics import ResidualBasedIncrementalUpdateStaticScheme as incremental_update_static_scheme

else:
    raise Exception("Distributed run requires TrilinosApplication")

def CreateDuplicateModelPart(origin_modelpart, destination_modelpart_name,
                             element_name, condition_name,
                             original_condition_name):
    # domain_size = origin_modelpart.ProcessInfo[Kratos.DOMAIN_SIZE]
    model = origin_modelpart.GetModel()
    connectivity_preserve_modeler = Kratos.ConnectivityPreserveModeler()

    if not model.HasModelPart(destination_modelpart_name):
        model_part = model.CreateModelPart(destination_modelpart_name)

        # TODO: Remove this line once the warnings from connectivity preserve modeller is gone, otherwise,
        #       output will be cluttered with lots of missing variable warnings
        RansVariableUtilities.CopyNodalSolutionStepVariablesList(
            origin_modelpart, model_part)

        # TODO: [PeriodicCondition]
        #       Currently, all the conditions will be replaced with the given new condition. This is an issue
        #       in the case of periodic cases in mpi, there we have to put PeriodicConditions in the mdpa file,
        #       where MetisParitioner will use that condition list to properly partition it. Therefore, "original_condition_name"
        #       is not used in this method at the moment.
        #       Following is one of the proposals to make PeriodicConditions to work with connectivity_preserve_modeller.
        # connectivity_preserve_modeler.GenerateModelPart(
        #     origin_modelpart, model_part, element_name, condition_name,
        #     original_condition_name + str(domain_size) + "D" + str(domain_size)
        #     + "N")

        connectivity_preserve_modeler.GenerateModelPart(
            origin_modelpart, model_part, element_name, condition_name)

    Kratos.Logger.PrintInfo("RANSModelPartFactory",
                            "Created " + destination_modelpart_name)
    return model.GetModelPart(destination_modelpart_name)


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


def CreateIncrementalUpdateScheme():
    return incremental_update_static_scheme()


def CalculateNormalsOnConditions(model_part):
    domain_size = model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
    if (not RansVariableUtilities.IsAnalysisStepCompleted(
            model_part, "CONDITION_NORMAL_CALCULATION")):

        # this calculates normals on whole model part, and assigns
        # NORMAL variable in NodalSolutionStepDataValue container.
        # NORMAL on conditions is required by some formulations such as inlet condition for
        # incompressible potential flow velocity formulation, and all boundaries for incompressible
        # potential flow pressure formulation.
        Kratos.NormalCalculationUtils().CalculateOnSimplex(
            model_part.Conditions, domain_size)
        RansVariableUtilities.AddAnalysisStep(model_part,
                                              "CONDITION_NORMAL_CALCULATION")

        # This reverts incorrectly calculated nodal NORMALS from previous method
        # since it spreads condition NORMAL to all nodes of model part, but from this
        # method, it again spreads condition NORMALs to nodes where condition is applied
        # with SLIP flag.
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
    return CreateDuplicateModelPart(
        formulation.GetBaseModelPart(),
        formulation.GetName() + "_" + element_name + "_" + condition_name,
        element_name, condition_name, "")
