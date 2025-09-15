import KratosMultiphysics as Kratos

try:
    import KratosMultiphysics.HDF5Application as KratosHDF5
except ImportError:
    KratosHDF5 = None

def _GetListOfVariableNames(list_of_variables):
    list_of_variable_names = []
    for var in list_of_variables:
        if isinstance(var, str):
            list_of_variable_names.append(var)
        else:
            list_of_variable_names.append(var.Name())
    return list_of_variable_names

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