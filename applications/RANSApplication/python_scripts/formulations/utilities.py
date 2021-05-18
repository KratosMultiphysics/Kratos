from importlib import import_module

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
    periodic_condition_name = "PeriodicCondition"):

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
                periodic_condition_name, index, node_id_list, properties)
            periodic_condition.Set(Kratos.PERIODIC)

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
