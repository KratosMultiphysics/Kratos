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
        for i in range(len(value_a)):
            test_unit.__CheckValues(value_a[i], value_b[i], tolerance)
    else:
        if (isinstance(value_a, Kratos.Matrix)):
            test_unit.assertMatrixAlmostEqual(value_a, value_b, tolerance)
        elif (isinstance(value_a, Kratos.Vector)):
            test_unit.assertVectorAlmostEqual(value_a, value_b, tolerance)
        elif (isinstance(value_a, list)):
            for i in range(len(value_a)):
                test_unit.assertAlmostEqual(value_a[i], value_b[i], tolerance)
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

    model_part.CreateNewCondition("Condition2D2N", 1, [1, 2], prop)
    model_part.CreateNewCondition("Condition2D2N", 2, [2, 3], prop)
    model_part.CreateNewCondition("Condition2D2N", 3, [3, 4], prop)
    model_part.CreateNewCondition("Condition2D2N", 4, [4, 5], prop)
    model_part.CreateNewCondition("Condition2D2N", 5, [5, 6], prop)
    model_part.CreateNewCondition("Condition2D2N", 6, [6, 7], prop)
    model_part.CreateNewCondition("Condition2D2N", 7, [7, 1], prop)

def InitializeModelPartVariables(model_part):
    for node in model_part.Nodes:
        node.SetValue(Kratos.PRESSURE, GetRandomValue())
        node.SetValue(Kratos.VELOCITY, GetRandomVector())
        node.SetValue(Kratos.LOAD_MESHES, GetRandomVector(5))
        node.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, GetRandomMatrix())

        node.SetSolutionStepValue(Kratos.PRESSURE, 0, GetRandomValue())
        node.SetSolutionStepValue(Kratos.VELOCITY, 0, GetRandomVector())
        node.SetSolutionStepValue(Kratos.LOAD_MESHES, 0, GetRandomVector(5))
        node.SetSolutionStepValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, 0, GetRandomMatrix())

    for condition in model_part.Conditions:
        condition.SetValue(Kratos.PRESSURE, GetRandomValue())
        condition.SetValue(Kratos.VELOCITY, GetRandomVector())
        condition.SetValue(Kratos.LOAD_MESHES, GetRandomVector(5))
        condition.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, GetRandomMatrix())

    for element in model_part.Elements:
        element.SetValue(Kratos.PRESSURE, GetRandomValue())
        element.SetValue(Kratos.VELOCITY, GetRandomVector())
        element.SetValue(Kratos.LOAD_MESHES, GetRandomVector(5))
        element.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, GetRandomMatrix())

def GetInitialVariableValue(variable, norm_type):
    if (norm_type == "none"):
        if (variable == Kratos.PRESSURE):
            return 0.0
        elif (variable == Kratos.VELOCITY):
            return Kratos.Vector(3, 0.0)
        elif (variable == Kratos.LOAD_MESHES):
            return Kratos.Vector(5, 0.0)
        elif (variable == Kratos.GREEN_LAGRANGE_STRAIN_TENSOR):
            return Kratos.Matrix(5, 5, 0.0)
    else:
        return 0.0

def InitializeProcesses(test):
    for process in test.process_list:
        process.Check()
    for process in test.process_list:
        process.ExecuteInitialize()

def ExecuteProcessFinalizeSolutionStep(test):
    for process in test.process_list:
        process.ExecuteFinalizeSolutionStep()