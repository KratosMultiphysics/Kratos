import os, sys
from importlib import import_module
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics import VariableUtils
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.RANSApplication import RansCalculationUtilities

if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.TrilinosApplication import TrilinosBlockBuilderAndSolverPeriodic
    from KratosMultiphysics.TrilinosApplication import TrilinosBlockBuilderAndSolver
elif (not IsDistributedRun()):
    from KratosMultiphysics import ResidualBasedBlockBuilderAndSolver
    from KratosMultiphysics.FluidDynamicsApplication import ResidualBasedBlockBuilderAndSolverPeriodic
else:
    raise Exception("Distributed run requires TrilinosApplication")


def GetKratosObjectPrototype(type_name):
    type_dict = {
        "LinearSolverFactory": [
            "KratosMultiphysics.python_linear_solver_factory.ConstructSolver",
            "KratosMultiphysics.TrilinosApplication.trilinos_linear_solver_factory.ConstructSolver"
        ],
        "ResidualBasedLinearStrategy": [
            "KratosMultiphysics.ResidualBasedLinearStrategy",
            "KratosMultiphysics.TrilinosApplication.TrilinosLinearStrategy"
        ],
        "ResidualBasedNewtonRaphsonStrategy": [
            "KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy",
            "KratosMultiphysics.TrilinosApplication.TrilinosNewtonRaphsonStrategy"
        ],
        "MixedGenericCriteria": [
            "KratosMultiphysics.MixedGenericCriteria",
            "KratosMultiphysics.TrilinosApplication.TrilinosMixedGenericCriteria"
        ],
        "ResidualBasedIncrementalUpdateStaticScheme": [
            "KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme",
            "KratosMultiphysics.TrilinosApplication.TrilinosResidualBasedIncrementalUpdateStaticScheme"
        ],
        "SteadyScalarScheme": [
            "KratosMultiphysics.RANSApplication.SteadyScalarScheme",
            "KratosMultiphysics.RANSApplication.TrilinosExtension.MPISteadyScalarScheme"
        ],
        "AlgebraicFluxCorrectedSteadyScalarScheme": [
            "KratosMultiphysics.RANSApplication.AlgebraicFluxCorrectedSteadyScalarScheme",
            "KratosMultiphysics.RANSApplication.TrilinosExtension.MPIAlgebraicFluxCorrectedSteadyScalarScheme"
        ],
        "BossakRelaxationScalarScheme": [
            "KratosMultiphysics.RANSApplication.BossakRelaxationScalarScheme",
            "KratosMultiphysics.RANSApplication.TrilinosExtension.MPIBossakRelaxationScalarScheme"
        ],
        "ResidualBasedSimpleSteadyScheme": [
            "KratosMultiphysics.FluidDynamicsApplication.ResidualBasedSimpleSteadyScheme",
            "KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension.TrilinosResidualBasedSimpleSteadyScheme"
        ],
        "ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent":[
            "KratosMultiphysics.FluidDynamicsApplication.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent",
            "KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension.TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent"
        ],
        "FractionalStepSettingsPeriodic":[
            "KratosMultiphysics.FluidDynamicsApplication.FractionalStepSettingsPeriodic",
            "KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension.TrilinosFractionalStepSettingsPeriodic"
        ],
        "FractionalStepSettings":[
            "KratosMultiphysics.FluidDynamicsApplication.FractionalStepSettings",
            "KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension.TrilinosFractionalStepSettings"
        ],
        "FractionalStepStrategy":[
            "KratosMultiphysics.FluidDynamicsApplication.FractionalStepStrategy",
            "KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension.TrilinosFractionalStepStrategy"
        ],
        "StrategyLabel":[
            "KratosMultiphysics.FluidDynamicsApplication.StrategyLabel",
            "KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension.TrilinosStrategyLabel"
        ],
        "SimpleSteadyAdjointScheme":[
            "KratosMultiphysics.FluidDynamicsApplication.SimpleSteadyAdjointScheme",
            "KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension.TrilinosSimpleSteadyAdjointScheme"
        ],
        "VelocityBossakAdjointScheme":[
            "KratosMultiphysics.FluidDynamicsApplication.VelocityBossakAdjointScheme",
            "KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension.TrilinosVelocityBossakAdjointScheme"
        ]
    }

    if (type_name not in type_dict.keys()):
        raise Exception(type_name + " not found in type_dict. Followings are allowed type_names:\n\t" + "\n\t".join(sorted(type_dict.keys())))

    module_info = type_dict[type_name][IsDistributedRun()]
    index = module_info.rfind(".")
    module_name = module_info[:index]
    attribute_name = module_info[index + 1:]

    module = import_module(module_name)
    return getattr(module, attribute_name)


def CreateDuplicateModelPart(
    origin_modelpart,
    destination_modelpart_name,
    element_name,
    condition_name):
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

        if (condition_name != ""):
            connectivity_preserve_modeler.GenerateModelPart(
                origin_modelpart, model_part, element_name, condition_name)
        else:
            connectivity_preserve_modeler.GenerateModelPart(
                origin_modelpart, model_part, element_name)

    Kratos.Logger.PrintInfo("RANSModelPartFactory",
                            "Created " + destination_modelpart_name)
    return model.GetModelPart(destination_modelpart_name)


def CreateRansFormulationModelPart(
    original_model_part,
    model_part_name_suffix,
    domain_size,
    element_name,
    condition_name = ""):

    element_suffix = str(domain_size) + "D" + str(domain_size + 1) + "N"
    element_name = element_name + element_suffix

    new_model_part_name = model_part_name_suffix + "_" + element_name

    if (condition_name != ""):
        condition_suffix = str(domain_size) + "D" + str(
                               domain_size) + "N"
        condition_name = condition_name + condition_suffix
        new_model_part_name += "_" + condition_name

    return CreateDuplicateModelPart(original_model_part,
                                    new_model_part_name, element_name,
                                    condition_name)


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


def InitializeYPlusVariablesInConditions(model_part):
    if (not RansVariableUtilities.IsAnalysisStepCompleted(
            model_part, "CONDITION_TURBULENCE_VARIABLE_INITIALIZATION")):
        VariableUtils().SetNonHistoricalVariableToZero(KratosRANS.RANS_Y_PLUS, model_part.Conditions)
        VariableUtils().SetNonHistoricalVariableToZero(KratosRANS.FRICTION_VELOCITY, model_part.Conditions)
        RansVariableUtilities.AddAnalysisStep(model_part,
                                              "CONDITION_TURBULENCE_VARIABLE_INITIALIZATION")
        Kratos.Logger.PrintInfo("Initialization",
                                "Initialized condition variables.")


def InitializePeriodicConditions(
    base_model_part,
    model_part,
    variables_list,
    periodic_condition_name = "PeriodicCondition",
    create_new_conditions = True):

    properties = model_part.CreateNewProperties(
        model_part.NumberOfProperties() + 1)

    pcu = KratosCFD.PeriodicConditionUtilities(
        model_part, model_part.ProcessInfo[Kratos.DOMAIN_SIZE])
    for variable in variables_list:
        pcu.AddPeriodicVariable(properties, variable)

    def modify_condition_properties(condition, _):
        condition.SetProperties(properties)

    def create_new_condition(condition, index):
        node_id_list = [node.Id for node in condition.GetNodes()]
        periodic_condition = model_part.CreateNewCondition(
            periodic_condition_name, index, node_id_list, properties)
        periodic_condition.Set(Kratos.PERIODIC)

    if (create_new_conditions):
        periodic_condition_update = create_new_condition
    else:
        periodic_condition_update = modify_condition_properties

    index = model_part.GetCommunicator().GlobalNumberOfConditions()
    for condition in base_model_part.Conditions:
        if condition.Is(Kratos.PERIODIC):
            index += 1
            periodic_condition_update(condition, index)

def InitializeWallLawProperties(model):
    for model_part_name in model.GetModelPartNames():
        model_part = model[model_part_name]
        process_info = model_part.ProcessInfo
        for properties in model_part.Properties:
            # logarithmic wall law
            if (properties.Has(KratosRANS.WALL_SMOOTHNESS_BETA) and process_info.Has(KratosRANS.VON_KARMAN)):
                von_karman = model_part.ProcessInfo[KratosRANS.VON_KARMAN]
                beta = properties[KratosRANS.WALL_SMOOTHNESS_BETA]
                y_plus_limit = RansCalculationUtilities.CalculateLogarithmicYPlusLimit(von_karman, beta)
                properties.SetValue(KratosRANS.RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT, y_plus_limit)


def GetBoundaryFlags(boundary_flags_parameters):
    if (boundary_flags_parameters.size == 0):
        raise Exception("No boundary flags were found")

    flags = Kratos.KratosGlobals.GetFlag(
        boundary_flags_parameters[0].GetString())
    for i in range(1, boundary_flags_parameters.size()):
        flags |= Kratos.KratosGlobals.GetFlag(
            boundary_flags_parameters[i].GetString())

    return (flags)


def GetConvergenceInfo(
    variable,
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


def CreateBlockBuilderAndSolver(
    linear_solver,
    is_periodic,
    communicator):
    if (IsDistributedRun()):
        if (is_periodic):
            return TrilinosBlockBuilderAndSolverPeriodic(
                communicator, 30, linear_solver,
                KratosCFD.PATCH_INDEX)
        else:
            return TrilinosBlockBuilderAndSolver(
                communicator, 30, linear_solver)
    else:
        if (is_periodic):
            return ResidualBasedBlockBuilderAndSolverPeriodic(
                linear_solver, KratosCFD.PATCH_INDEX)
        else:
            return ResidualBasedBlockBuilderAndSolver(linear_solver)

def GetTimeDerivativeVariable(variable):
    if (variable == Kratos.VELOCITY_X):
        return Kratos.ACCELERATION_X
    elif (variable == Kratos.VELOCITY_Y):
        return Kratos.ACCELERATION_Y
    elif (variable == Kratos.VELOCITY_Z):
        return Kratos.ACCELERATION_Z
    elif (variable == KratosRANS.TURBULENT_KINETIC_ENERGY):
        return KratosRANS.TURBULENT_KINETIC_ENERGY_RATE
    elif (variable == KratosRANS.TURBULENT_KINETIC_ENERGY_RATE):
        return KratosRANS.RANS_AUXILIARY_VARIABLE_1
    elif (variable == KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE):
        return KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2
    elif (variable == KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2):
        return KratosRANS.RANS_AUXILIARY_VARIABLE_2
    elif (variable == KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE):
        return KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2
    elif (variable == KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2):
        return KratosRANS.RANS_AUXILIARY_VARIABLE_2
    elif (variable == KratosRANS.VELOCITY_POTENTIAL):
        return KratosRANS.VELOCITY_POTENTIAL_RATE
    elif (variable == KratosRANS.VELOCITY_POTENTIAL_RATE):
        return KratosRANS.RANS_AUXILIARY_VARIABLE_1
    else:
        return None

def GetTimeDerivativeVariablesRecursively(var):
    time_derivative_var = GetTimeDerivativeVariable(var)
    if (time_derivative_var is None):
        return [var]
    else:
        v = [var]
        v.extend(GetTimeDerivativeVariablesRecursively(time_derivative_var))
        return v

def AddFileLoggerOutput(log_file_name):
    file_logger = Kratos.FileLoggerOutput(log_file_name)
    default_severity = Kratos.Logger.GetDefaultOutput().GetSeverity()
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    Kratos.Logger.AddOutput(file_logger)

    return default_severity, file_logger

def RemoveFileLoggerOutput(default_severity, file_logger):
    Kratos.Logger.Flush()
    Kratos.Logger.RemoveOutput(file_logger)
    Kratos.Logger.GetDefaultOutput().SetSeverity(default_severity)

def SolveProblem(analysis_class_type, kratos_parameters, execution_prefix):
    # set the loggers
    default_severity, file_logger = AddFileLoggerOutput(execution_prefix + ".log")

    # run the primal analysis
    model = Kratos.Model()
    primal_simulation = analysis_class_type(model, kratos_parameters)
    primal_simulation.Run()

    with open(execution_prefix + ".json", "w") as file_output:
        file_output.write(kratos_parameters.PrettyPrintJsonString())

    # flush the primal output
    RemoveFileLoggerOutput(default_severity, file_logger)
    Kratos.Logger.PrintInfo("SolvePrimalProblem", "Solved primal evaluation at {}.".format(execution_prefix + ".json"))

    return model, primal_simulation

class ExecutionScope:
    def __init__(self, execution_path):
        self.currentPath = Path.cwd()
        self.scope = Path(execution_path)

    def __enter__(self):
        self.scope.mkdir(parents=True, exist_ok=True)
        sys.path.append(str(self.scope.absolute()))
        os.chdir(str(self.scope))

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)
        sys.path.remove(str(self.scope.absolute()))


def AddWallPropertiesUpdateProcess(rans_formulation, settings):
    wall_properties_update_execution_points = ["initialize"]
    if settings.Has("wall_properties_update_execution_points"):
        wall_properties_update_execution_points = settings["wall_properties_update_execution_points"].GetStringArray()

    update_wall_normals = True
    if settings.Has("update_wall_normals"):
        update_wall_normals = settings["update_wall_normals"].GetBool()

    update_wall_condition_heights = True
    if settings.Has("update_wall_condition_heights"):
        update_wall_condition_heights = settings["update_wall_condition_heights"].GetBool()

    echo_level = 0
    if settings.Has("echo_level"):
        echo_level = settings["echo_level"].GetInt()

    wall_model_part_name = "ALL_WALL_MODEL_PART"
    if settings.Has("wall_model_part_name"):
        wall_model_part_name = settings["wall_model_part_name"].GetString()

    base_model_part = rans_formulation.GetBaseModelPart()
    wall_model_part_full_name = f"{base_model_part.FullName()}.{wall_model_part_name}"
    if len(wall_properties_update_execution_points) > 0:
        wall_properties_update_process = KratosRANS.RansWallPropertiesUpdateProcess(
            rans_formulation.GetBaseModelPart().GetModel(),
            wall_model_part_full_name,
            update_wall_normals,
            update_wall_condition_heights,
            wall_properties_update_execution_points,
            echo_level)

        rans_formulation.AddProcess(wall_properties_update_process)
