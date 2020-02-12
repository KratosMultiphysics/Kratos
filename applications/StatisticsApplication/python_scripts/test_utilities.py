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