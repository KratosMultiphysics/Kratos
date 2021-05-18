import KratosMultiphysics as Kratos
from random import uniform


def GetRandomValue(min_value=-10.0, max_value=10.0):
    return uniform(min_value, max_value)


def GetRandomVector(vector_size=3, min_value=-10.0, max_value=10.0):
    v = Kratos.Vector(vector_size)

    for i in range(vector_size):
        v[i] = GetRandomValue(min_value, max_value)

    return v


def GetRandomMatrix(size_1=5, size_2=5, min_value=-10.0, max_value=10.0):
    m = Kratos.Matrix(size_1, size_2)

    for i in range(size_1):
        for j in range(size_1):
            m[i, j] = GetRandomValue(min_value, max_value)

    return m

def HistoricalRetrievalMethod(item, variable):
    return item.GetSolutionStepValue(variable)

def NonHistoricalRetrievalMethod(item, variable):
    return item.GetValue(variable)

def InitializeContainerArrays(container):
    scalar_list = []
    vector_3d_list = []
    vector_list = []
    matrix_list = []
    for _ in container:
        scalar_list.append([])
        vector_3d_list.append([])
        vector_list.append([])
        matrix_list.append([])

    return scalar_list, vector_3d_list, vector_list, matrix_list

def CheckValues(test_unit, value_a, value_b, tolerance):
    if (isinstance(value_a, tuple)):
        for i, v_a in enumerate(value_a):
            CheckValues(test_unit, v_a, value_b[i], tolerance)
    else:
        if (isinstance(value_a, Kratos.Matrix)):
            test_unit.assertMatrixAlmostEqual(value_a, value_b, tolerance)
        elif (isinstance(value_a, Kratos.Vector)):
            test_unit.assertVectorAlmostEqual(value_a, value_b, tolerance)
        elif (isinstance(value_a, list)):
            for i, v_a in enumerate(value_a):
                test_unit.assertAlmostEqual(v_a, value_b[i], tolerance)
        else:
            test_unit.assertAlmostEqual(value_a, value_b, tolerance)

def CreateModelPart(model_part):
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
    model_part.CreateNewNode(4, 2.0, 1.0, 0.0)
    model_part.CreateNewNode(5, 1.5, 1.0, 0.0)
    model_part.CreateNewNode(6, 0.5, 1.0, 0.0)
    model_part.CreateNewNode(7, 0.0, 1.0, 0.0)

    prop = model_part.GetProperties()[0]

    model_part.CreateNewElement("Element2D3N", 1, [1, 6, 7], prop)
    model_part.CreateNewElement("Element2D3N", 2, [1, 2, 6], prop)
    model_part.CreateNewElement("Element2D3N", 3, [6, 2, 5], prop)
    model_part.CreateNewElement("Element2D3N", 4, [2, 3, 5], prop)
    model_part.CreateNewElement("Element2D3N", 5, [5, 3, 4], prop)

    model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], prop)
    model_part.CreateNewCondition("LineCondition2D2N", 2, [2, 3], prop)
    model_part.CreateNewCondition("LineCondition2D2N", 3, [3, 4], prop)
    model_part.CreateNewCondition("LineCondition2D2N", 4, [4, 5], prop)
    model_part.CreateNewCondition("LineCondition2D2N", 5, [5, 6], prop)
    model_part.CreateNewCondition("LineCondition2D2N", 6, [6, 7], prop)
    model_part.CreateNewCondition("LineCondition2D2N", 7, [7, 1], prop)

def InitializeModelPartVariables(model_part, is_random = True):
    if is_random:
        scalar_method = lambda item : GetRandomValue()
        array_3d_method = lambda item : GetRandomVector()
        vector_method = lambda item : GetRandomVector(5)
        matrix_method = lambda item : GetRandomMatrix()
    else:
        t = model_part.ProcessInfo[Kratos.STEP]
        scalar_method = lambda item : GetInitialVariableValue(Kratos.PRESSURE, "none", (t+1) * (item.Id+1))
        array_3d_method = lambda item : GetInitialVariableValue(Kratos.VELOCITY, "none", (t+1) * (item.Id +1) * 2.0)
        vector_method = lambda item : GetInitialVariableValue(Kratos.LOAD_MESHES, "none", (t+1) * (item.Id + 2.0))
        matrix_method = lambda item : GetInitialVariableValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, "none", (t+1) * (item.Id+3)* 2.0)

    communicator = model_part.GetCommunicator()
    container = communicator.LocalMesh()

    for node in container.Nodes:
        node.SetValue(Kratos.PRESSURE, scalar_method(node))
        node.SetValue(Kratos.VELOCITY, array_3d_method(node))
        node.SetValue(Kratos.LOAD_MESHES, vector_method(node))
        node.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, matrix_method(node))

        node.SetSolutionStepValue(Kratos.PRESSURE, 0, scalar_method(node))
        node.SetSolutionStepValue(Kratos.VELOCITY, 0, array_3d_method(node))
        node.SetSolutionStepValue(Kratos.LOAD_MESHES, 0, vector_method(node))
        node.SetSolutionStepValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, 0, matrix_method(node))

    communicator.GetDataCommunicator().Barrier()
    model_part.GetCommunicator().SynchronizeVariable(Kratos.PRESSURE)
    model_part.GetCommunicator().SynchronizeVariable(Kratos.VELOCITY)
    model_part.GetCommunicator().SynchronizeVariable(Kratos.LOAD_MESHES)
    model_part.GetCommunicator().SynchronizeVariable(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
    model_part.GetCommunicator().SynchronizeNonHistoricalVariable(Kratos.PRESSURE)
    model_part.GetCommunicator().SynchronizeNonHistoricalVariable(Kratos.VELOCITY)
    model_part.GetCommunicator().SynchronizeNonHistoricalVariable(Kratos.LOAD_MESHES)
    model_part.GetCommunicator().SynchronizeNonHistoricalVariable(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

    for condition in container.Conditions:
        condition.SetValue(Kratos.PRESSURE, scalar_method(condition))
        condition.SetValue(Kratos.VELOCITY, array_3d_method(condition))
        condition.SetValue(Kratos.LOAD_MESHES, vector_method(condition))
        condition.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, matrix_method(condition))

    for element in container.Elements:
        element.SetValue(Kratos.PRESSURE, scalar_method(element))
        element.SetValue(Kratos.VELOCITY, array_3d_method(element))
        element.SetValue(Kratos.LOAD_MESHES, vector_method(element))
        element.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, matrix_method(element))

    communicator.GetDataCommunicator().Barrier()

def GetInitialVariableValue(variable, norm_type, value = 0.0):
    if (norm_type == "none"):
        if (variable == Kratos.PRESSURE):
            return value
        elif (variable == Kratos.VELOCITY):
            return Kratos.Vector(3, value)
        elif (variable == Kratos.LOAD_MESHES):
            return Kratos.Vector(5, value)
        elif (variable == Kratos.GREEN_LAGRANGE_STRAIN_TENSOR):
            return Kratos.Matrix(5, 5, value)
    else:
        return 0.0

def InitializeProcesses(test):
    communicator = test.model_part.GetCommunicator().GetDataCommunicator()
    for process in test.process_list:
        communicator.Barrier()
        process.Check()

    for process in test.process_list:
        communicator.Barrier()
        process.ExecuteInitialize()

def ExecuteProcessFinalizeSolutionStep(test):
    communicator = test.model_part.GetCommunicator().GetDataCommunicator()
    for process in test.process_list:
        communicator.Barrier()
        process.ExecuteFinalizeSolutionStep()
