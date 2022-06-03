import sympy

class CompressibleNavierStokesDefines:
    matrix_format = None
    vector_format = None

    @classmethod
    def SetFormat(cls, language):
        formats = {
            'python': ("{name}[{i},{j}]", "{name}[{i}]"),
            'c'     : ("{name}({i},{j})", "{name}({i})")
        }

        if language in formats:
            cls.matrix_format, cls.vector_format = formats[language]
            return

        msg = "Only formats implemented are for languages:\n"
        for lang in formats:
            msg += " - " + lang + "\n"
        raise ValueError(msg)

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
