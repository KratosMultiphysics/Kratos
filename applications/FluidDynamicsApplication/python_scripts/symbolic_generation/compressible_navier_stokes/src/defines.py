import sympy

class CompressibleNavierStokesDefines:
    matrix_format = "{name}({i},{j})"
    vector_format = "{name}({i})"

    @classmethod
    def Matrix(cls, name, n, m, **kwargs):
        return sympy.Matrix(n, m, lambda i, j: sympy.Symbol(
            cls.matrix_format.format(name=name, i=i, j=j), **kwargs))

    @classmethod
    def Vector(cls, name, n, **kwargs):
        return sympy.Matrix(n, 1, lambda i, _: sympy.Symbol(
            cls.vector_format.format(name=name, i=i), **kwargs))

    @classmethod
    def ZeroMatrix(cls, rows, cols):
        return sympy.Matrix(rows, cols, lambda *_: 0.0)

    @classmethod
    def ZeroVector(cls, rows):
        return cls.ZeroMatrix(rows, 1)
