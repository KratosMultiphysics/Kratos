import sympy


def DefineMatrix(name, n, m):
    return sympy.Matrix(n, m, lambda i, j: sympy.Symbol("{}({},{})".format(name, i, j), real=True))


def DefineVector(name, n):
    return sympy.Matrix(n, 1, lambda i, _: sympy.Symbol("{}({})".format(name, i), real=True))


def ZeroMatrix(rows, cols):
    return sympy.Matrix(rows, cols, lambda *_: 0.0)


def ZeroVector(rows):
    return ZeroMatrix(rows, 1)
