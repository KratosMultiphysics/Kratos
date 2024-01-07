from math import sqrt

import KratosMultiphysics as Kratos
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializeWallLawProperties

try:
    import KratosMultiphysics.HDF5Application as KratosHDF5
except ImportError:
    KratosHDF5 = None

try:
    import KratosMultiphysics.RANSApplication as KratosRANS
except ImportError:
    KratosRANS = None

def _Array3DInnerProduct(vec_1, vec_2):
    return vec_1[0]*vec_2[0] + vec_1[1]*vec_2[1] + vec_1[2]*vec_2[2]

def _GetListOfVariableNames(list_of_variables):
    list_of_variable_names = []
    for var in list_of_variables:
        if isinstance(var, str):
            list_of_variable_names.append(var)
        else:
            list_of_variable_names.append(var.Name())
    return list_of_variable_names

def _GetListOfVariables(list_of_variable_names):
    list_of_variables = []
    for var in list_of_variable_names:
        if isinstance(var, str):
            list_of_variables.append(Kratos.KratosGlobals.GetVariable(var))
        else:
            list_of_variables.append(var.Name())
    return list_of_variables

def _SetParameterValue(parameter, value):
    if isinstance(value, str):
        parameter.SetString(value)
    elif isinstance(value, int):
        parameter.SetInt(value)
    elif isinstance(value, float):
        parameter.SetDouble(value)
    elif isinstance(value, bool):
        parameter.SetBool(value)
    else:
        raise Exception("Unsupported value type")

def _SetSubParameter(parameters, sub_parameter_name, value):
    if (not parameters.Has(sub_parameter_name)):
        parameters.AddEmptyValue(sub_parameter_name)
    _SetParameterValue(parameters[sub_parameter_name], value)

def CreateModelPartFromMDPA(model, mdpa_file_name, list_of_solution_step_variables, model_part_name = "rans_post_processing", buffer_size = 1):
    model_part = model.CreateModelPart(model_part_name)
    model_part.SetBufferSize(buffer_size)

    for var in _GetListOfVariables(list_of_solution_step_variables):
        model_part.AddNodalSolutionStepVariable(var)

    Kratos.ModelPartIO(mdpa_file_name, Kratos.IO.READ).ReadModelPart(model_part)
    return model_part

def CalculateNodalDomainSize(model_part, output_variable, is_historical = True):
    output_variable_type = Kratos.KratosGlobals.GetVariableType(output_variable.Name())
    if (output_variable_type != "Double"):
        raise Exception("Output variable should be of type double. [ Output variable name = {:s}, output variable type = {:s}".format(output_variable.Name(), output_variable_type))

    if is_historical:
        Kratos.VariableUtils().SetHistoricalVariableToZero(output_variable, model_part.Nodes)
        def update_method(node, nodal_mass):
            node.SetSolutionStepValue(output_variable, node.GetSolutionStepValue(output_variable) + nodal_mass)
    else:
        Kratos.VariableUtils().SetNonHistoricalVariableToZero(output_variable, model_part.Nodes)
        def update_method(node, nodal_mass):
            node.SetValue(output_variable, node.GetValue(output_variable) + nodal_mass)

    for element in model_part.Elements:
        geometry = element.GetGeometry()
        nodal_mass = geometry.Area() / geometry.PointsNumber()
        for node in geometry:
            update_method(node, nodal_mass)

    if is_historical:
        model_part.GetCommunicator().AssembleCurrentData(output_variable)
    else:
        model_part.GetCommunicator().AssembleNonHistoricalData(output_variable)

def CalculateConditionNeighboursAndNormals(model_part):
    tmoc = Kratos.TetrahedralMeshOrientationCheck
    throw_errors = False
    flags = (tmoc.COMPUTE_NODAL_NORMALS) | (tmoc.COMPUTE_CONDITION_NORMALS) | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
    Kratos.TetrahedralMeshOrientationCheck(model_part, throw_errors, flags).Execute()

def GetHDF5File(hdf5_file_name, file_access_mode):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_file_settings = Kratos.Parameters("""{
        "file_name"       : "",
        "file_access_mode": "read_only"
    }""")
    hdf5_file_settings["file_name"].SetString(hdf5_file_name)
    hdf5_file_settings["file_access_mode"].SetString(file_access_mode)
    return KratosHDF5.HDF5FileSerial(hdf5_file_settings)

def GetHDF5FileForReading(hdf5_file_name):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_file_settings = Kratos.Parameters("""{
        "file_name"       : "",
        "file_access_mode": "read_only"
    }""")
    hdf5_file_settings["file_name"].SetString(hdf5_file_name)
    return KratosHDF5.HDF5FileSerial(hdf5_file_settings)

def GetHDF5FileForWriting(hdf5_file_name):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_file_settings = Kratos.Parameters("""{
        "file_name"       : "",
        "file_access_mode": "truncate"
    }""")
    hdf5_file_settings["file_name"].SetString(hdf5_file_name)
    return KratosHDF5.HDF5FileSerial(hdf5_file_settings)

def InputNodalResultsFromHDF5(model_part, hdf5_file, list_of_variables, is_historical = True):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_input_parameters = Kratos.Parameters("""
    {
        "prefix": "/ResultsData",
        "list_of_variables":[]
    }""")
    hdf5_input_parameters["list_of_variables"].SetStringArray(_GetListOfVariableNames(list_of_variables))

    if is_historical:
        nodal_io = KratosHDF5.HDF5NodalSolutionStepDataIO(
            hdf5_input_parameters, hdf5_file)
        nodal_io.ReadNodalResults(model_part, 0)
    else:
        nodal_io = KratosHDF5.HDF5NodalDataValueIO(
            hdf5_input_parameters,
            hdf5_file)
        nodal_io.ReadNodalResults(model_part.Nodes, model_part.GetCommunicator())

def InputNodalFlagsFromHDF5(model_part, hdf5_file, list_of_variables):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_input_parameters = Kratos.Parameters("""
    {
        "prefix": "/ResultsData",
        "list_of_variables":[]
    }""")
    hdf5_input_parameters["list_of_variables"].SetStringArray(_GetListOfVariableNames(list_of_variables))

    nodal_io = KratosHDF5.HDF5NodalFlagValueIO(
        hdf5_input_parameters,
        hdf5_file)
    nodal_io.ReadNodalFlags(model_part.Nodes, model_part.GetCommunicator())

def InputConditionResultsFromHDF5(model_part, hdf5_file, list_of_variables):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_input_parameters = Kratos.Parameters("""
    {
        "prefix": "/ResultsData",
        "list_of_variables":[]
    }""")
    hdf5_input_parameters["list_of_variables"].SetStringArray(_GetListOfVariableNames(list_of_variables))

    nodal_io = KratosHDF5.HDF5ConditionDataValueIO(
        hdf5_input_parameters,
        hdf5_file)
    nodal_io.ReadConditionResults(model_part.Conditions, model_part.GetCommunicator())

def InputConditionFlagsFromHDF5(model_part, hdf5_file, list_of_variables):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_input_parameters = Kratos.Parameters("""
    {
        "prefix": "/ResultsData",
        "list_of_variables":[]
    }""")
    hdf5_input_parameters["list_of_variables"].SetStringArray(_GetListOfVariableNames(list_of_variables))

    nodal_io = KratosHDF5.HDF5ConditionFlagValueIO(
        hdf5_input_parameters,
        hdf5_file)
    nodal_io.ReadConditionFlags(model_part.Conditions, model_part.GetCommunicator())

def InputElementResultsFromHDF5(model_part, hdf5_file, list_of_variables):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_input_parameters = Kratos.Parameters("""
    {
        "prefix": "/ResultsData",
        "list_of_variables":[]
    }""")
    hdf5_input_parameters["list_of_variables"].SetStringArray(_GetListOfVariableNames(list_of_variables))

    nodal_io = KratosHDF5.HDF5ElementDataValueIO(
        hdf5_input_parameters,
        hdf5_file)
    nodal_io.ReadElementResults(model_part.Elements, model_part.GetCommunicator())

def InputElementFlagsFromHDF5(model_part, hdf5_file, list_of_variables):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_input_parameters = Kratos.Parameters("""
    {
        "prefix": "/ResultsData",
        "list_of_variables":[]
    }""")
    hdf5_input_parameters["list_of_variables"].SetStringArray(_GetListOfVariableNames(list_of_variables))

    nodal_io = KratosHDF5.HDF5ElementFlagValueIO(
        hdf5_input_parameters,
        hdf5_file)
    nodal_io.ReadElementFlags(model_part.Elements, model_part.GetCommunicator())

def OutputModelPartToHDF5(model_part, hdf5_file):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    KratosHDF5.HDF5ModelPartIO(hdf5_file, "/ModelData").WriteModelPart(model_part)

def OutputNodalResultsToHDF5(model_part, hdf5_file, list_of_variables, is_historical = True, step = 0):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_output_parameters = Kratos.Parameters("""
    {
        "prefix": "/ResultsData",
        "list_of_variables":[]
    }""")
    hdf5_output_parameters["list_of_variables"].SetStringArray(_GetListOfVariableNames(list_of_variables))

    if is_historical:
        nodal_io = KratosHDF5.HDF5NodalSolutionStepDataIO(
            hdf5_output_parameters,
            hdf5_file)
        nodal_io.WriteNodalResults(model_part, step)
    else:
        nodal_io = KratosHDF5.HDF5NodalDataValueIO(
            hdf5_output_parameters,
            hdf5_file)
        nodal_io.WriteNodalResults(model_part.Nodes)

def OutputConditionResultsToHDF5(model_part, hdf5_file, list_of_variables):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_output_parameters = Kratos.Parameters("""
    {
        "prefix": "/ResultsData",
        "list_of_variables":[]
    }""")
    hdf5_output_parameters["list_of_variables"].SetStringArray(_GetListOfVariableNames(list_of_variables))

    nodal_io = KratosHDF5.HDF5ConditionDataValueIO(
        hdf5_output_parameters,
        hdf5_file)
    nodal_io.WriteConditionResults(model_part.Conditions)

def OutputElementResultsToHDF5(model_part, hdf5_file, list_of_variables):
    if (KratosHDF5 is None):
        raise Exception("Please compile and install the HDF5 application first.")

    hdf5_output_parameters = Kratos.Parameters("""
    {
        "prefix": "/ResultsData",
        "list_of_variables":[]
    }""")
    hdf5_output_parameters["list_of_variables"].SetStringArray(_GetListOfVariableNames(list_of_variables))

    nodal_io = KratosHDF5.HDF5ElementDataValueIO(
        hdf5_output_parameters,
        hdf5_file)
    nodal_io.WriteElementResults(model_part.Elements)

def CreateLineOutput(model_part, line_output_parameters):
    if (KratosRANS is None):
        raise Exception("Please compile and install the RANS application first.")

    model_part.ProcessInfo[Kratos.STEP] = 1

    _SetSubParameter(line_output_parameters, "model_part_name", model_part.FullName())
    _SetSubParameter(line_output_parameters, "output_step_interval", 0)
    _SetSubParameter(line_output_parameters, "output_step_control_variable_name", "STEP")

    line_output_process = KratosRANS.RansLineOutputProcess(model_part.GetModel(), line_output_parameters)
    line_output_process.Check()
    line_output_process.ExecuteInitialize()
    line_output_process.ExecuteInitializeSolutionStep()
    line_output_process.ExecuteFinalizeSolutionStep()

def CalculateReactions(model_part, wall_model_part_name, materials_filename):
    # In here, we assume model_part process info constants are properly defined.

    # find parent elements
    tmoc = Kratos.TetrahedralMeshOrientationCheck
    throw_errors = False
    flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS) | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
    Kratos.TetrahedralMeshOrientationCheck(model_part, throw_errors, flags).Execute()

    # calculate normals
    KratosRANS.RansVariableUtilities.CalculateNodalNormal(model_part.GetModel()[wall_model_part_name])

    # populate material properties
    material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": ""}} """)
    material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
    Kratos.ReadMaterialsUtility(material_settings, model_part.GetModel())

    # add wall law properties
    InitializeWallLawProperties(model_part.GetModel())

    # initialize constitutive laws
    KratosRANS.RansVariableUtilities.SetElementConstitutiveLaws(model_part.Elements)

    # calculate reactions
    params = Kratos.Parameters("""
    {
        "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "execution_points": ["execute"],
        "echo_level"      : 1
    }""")
    params["model_part_name"].SetString(wall_model_part_name)
    process = KratosRANS.RansComputeReactionsProcess(model_part.GetModel(), params)
    process.Check()
    process.Execute()


def CalculateFirstElementYPlusValues(model_part, output_variable, density, kinmeatic_viscosity):
    if (KratosRANS is None):
        raise Exception("Please compile and install the RANS application first.")

    KratosRANS.RansVariableUtilities.CalculateNodalNormal(model_part)

    output_variable = _GetListOfVariables([output_variable])[0]
    for condition in model_part.Conditions:
        normal = condition.GetValue(Kratos.NORMAL)
        normal_mag = sqrt(_Array3DInnerProduct(normal, normal))
        normal = normal / normal_mag
        y = KratosRANS.RansCalculationUtilities.CalculateWallHeight(condition, normal)

        reaction_per_area = Kratos.Array3(0.0)
        for node in condition.GetGeometry():
            nodal_normal = node.GetSolutionStepValue(Kratos.NORMAL)
            nodal_normal_mag = sqrt(_Array3DInnerProduct(nodal_normal, nodal_normal))
            reaction_per_area += node.GetSolutionStepValue(Kratos.REACTION) / nodal_normal_mag

        reaction_per_area /= len(condition.GetGeometry())

        perpendicular_reaction = normal * _Array3DInnerProduct(reaction_per_area, normal)
        tangential_reaction = (reaction_per_area - perpendicular_reaction)
        shear_stress = sqrt(_Array3DInnerProduct(tangential_reaction, tangential_reaction))
        u_tau = sqrt(shear_stress / density)
        condition.SetValue(output_variable, u_tau * y / kinmeatic_viscosity)
