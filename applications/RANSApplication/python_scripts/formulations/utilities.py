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
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPISteadyScalarScheme as steady_scalar_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIAlgebraicFluxCorrectedSteadyScalarScheme as afc_steady_scalar_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIBossakScalarScheme as bossak_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import TrilinosRansWallDistanceCalculationProcess as wall_distance_calculation_process
elif (not IsDistributedRun()):
    from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
    from KratosMultiphysics.RANSApplication import GenericScalarConvergenceCriteria as residual_criteria
    # from Kratos import ResidualCriteria as residual_criteria
    from Kratos import ResidualBasedNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import PeriodicBlockBuilderAndSolver as periodic_block_builder_and_solver
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import BlockBuilderAndSolver as block_builder_and_solver
    from KratosMultiphysics import ResidualBasedIncrementalUpdateStaticScheme as incremental_update_static_scheme
    from KratosMultiphysics.RANSApplication import SteadyScalarScheme as steady_scalar_scheme
    from KratosMultiphysics.RANSApplication import AlgebraicFluxCorrectedSteadyScalarScheme as afc_steady_scalar_scheme
    from KratosMultiphysics.RANSApplication import BossakScalarScheme as bossak_scheme
    from KratosMultiphysics.RANSApplication import RansWallDistanceCalculationProcess as wall_distance_calculation_process
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


def IsBufferInitialized(formulation):
    return (formulation.GetBaseModelPart().ProcessInfo[Kratos.STEP] + 1 >=
            formulation.GetMinimumBufferSize())


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


def CreateSteadyScalarScheme(relaxation_factor):
    return steady_scalar_scheme(relaxation_factor)


def CreateSteadyAlgebraicFluxCorrectedScheme(relaxation_factor,
                                                     boundary_flags,
                                                     is_periodic):
    if (is_periodic):
        return afc_steady_scalar_scheme(relaxation_factor, boundary_flags,
                                        KratosCFD.PATCH_INDEX)
    else:
        return afc_steady_scalar_scheme(relaxation_factor, boundary_flags)


def CreateBossakScalarScheme(bossak_value, relaxation_factor, scalar_variable,
                             scalar_rate_variable,
                             relaxed_scalar_rate_variable):
    return bossak_scheme(bossak_value, relaxation_factor, scalar_variable,
                         scalar_rate_variable, relaxed_scalar_rate_variable)


def GetBoundaryFlags(boundary_flags_parameters):
    if (boundary_flags_parameters.size == 0):
        raise Exception("No boundary flags were found")

    flags = Kratos.KratosGlobals.GetFlag(
        boundary_flags_parameters[0].GetString())
    for i in range(1, boundary_flags_parameters.size()):
        flags |= Kratos.KratosGlobals.GetFlag(
            boundary_flags_parameters[i].GetString())

    return (flags)


def GetDefaultConditionName(model_part):
    process_info = model_part.ProcessInfo
    if (process_info.Has(Kratos.DOMAIN_SIZE)):
        domain_size = process_info[Kratos.DOMAIN_SIZE]
        if (domain_size == 2):
            return "LineCondition"
        elif (domain_size == 3):
            return "SurfaceCondition"
        else:
            raise Exception("Unsupported domain size. [ DOMAIN_SIZE = " +
                            str(domain_size) + " ].")
    else:
        raise Exception("DOMAIN_SIZE is not found in process info of " +
                        model_part.Name() + ".")


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

def InitializeYPlusVariablesInConditions(model_part):
    if (not RansVariableUtilities.IsAnalysisStepCompleted(
            model_part, "CONDITION_TURBULENCE_VARIABLE_INITIALIZATION")):
        VariableUtils().SetNonHistoricalVariableToZero(KratosRANS.RANS_Y_PLUS, model_part.Conditions)
        VariableUtils().SetNonHistoricalVariableToZero(KratosRANS.FRICTION_VELOCITY, model_part.Conditions)
        RansVariableUtilities.AddAnalysisStep(model_part,
                                              "CONDITION_TURBULENCE_VARIABLE_INITIALIZATION")
        Kratos.Logger.PrintInfo("Initialization",
                                "Initialized condition variables.")

def CreateWallDistanceCalculationProcess(model, settings):
    return wall_distance_calculation_process(model, settings)