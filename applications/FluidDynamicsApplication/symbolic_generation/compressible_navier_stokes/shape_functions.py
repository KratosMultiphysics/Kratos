import sympy
from sympy.functions.combinatorial.numbers import catalan


def DefineShapeFunctionsMatrix(dim, n_nodes, n_gauss):
    quadrature_data = {
        (2, 3): TriangleShapeFunctions,
        (2, 4): QuadrilateralShapeFunctions,
        (3, 4): TetrahedronShapeFunctions
    }
    try:
        return quadrature_data[(dim, n_nodes)](n_gauss)
    except KeyError as e:
        msg = "Geomerty {}D {}N is not implemented".format(dim, n_nodes)
        raise NotImplementedError(msg) from e


def TriangleShapeFunctions(n_gauss):
    mat_N = sympy.zeros(n_gauss, 3)
    if n_gauss == 1:
        return sympy.Matrix(n_gauss, 3, lambda *_: 1/3)
    if n_gauss == 3:
        mat_N[0, 0] = 2.0 / 3.0
        mat_N[0, 1] = 1.0 / 6.0
        mat_N[0, 2] = 1.0 / 6.0
        mat_N[1, 0] = 1.0 / 6.0
        mat_N[1, 1] = 2.0 / 3.0
        mat_N[1, 2] = 1.0 / 6.0
        mat_N[2, 0] = 1.0 / 6.0
        mat_N[2, 1] = 1.0 / 6.0
        mat_N[2, 2] = 2.0 / 3.0
        return mat_N
    msg = "Triangle with %d gauss points not implemented".format(n_gauss)
    raise NotImplementedError(msg)


def QuadrilateralShapeFunctions(n_gauss):
    N = [lambda xi, eta: 0.25 * (1 - xi) * (1 - eta),
         lambda xi, eta: 0.25 * (1 + xi) * (1 - eta),
         lambda xi, eta: 0.25 * (1 - xi) * (1 + eta),
         lambda xi, eta: 0.25 * (1 + xi) * (1 + eta)]

    if n_gauss == 1:
        return sympy.Matrix(n_gauss, 4, lambda _, j: N[j](0, 0))
    if n_gauss == 4:
        xi = [-3**-0.5,  3**0.5, -3**-0.5, 3**0.5]
        eta = [-3**-0.5, -3**0.5,  3**-0.5, 3**0.5]
        return sympy.Matrix(n_gauss, 4, lambda i, j: N[j](xi[i], eta[i]))

    msg = "Quadrilateral with %d gauss points not implemented".format(n_gauss)
    raise NotImplementedError(msg)


def TetrahedronShapeFunctions(n_gauss):
    if n_gauss == 1:
        return sympy.Matrix(n_gauss, 4, lambda *_: 0.25)
    if n_gauss == 4:
        mat_N = sympy.zeros(n_gauss,  4)
        mat_N[0, 0] = 0.58541020
        mat_N[0, 1] = 0.13819660
        mat_N[0, 2] = 0.13819660
        mat_N[0, 3] = 0.13819660
        mat_N[1, 0] = 0.13819660
        mat_N[1, 1] = 0.58541020
        mat_N[1, 2] = 0.13819660
        mat_N[1, 3] = 0.13819660
        mat_N[2, 0] = 0.13819660
        mat_N[2, 1] = 0.13819660
        mat_N[2, 2] = 0.58541020
        mat_N[2, 3] = 0.13819660
        mat_N[3, 0] = 0.13819660
        mat_N[3, 1] = 0.13819660
        mat_N[3, 2] = 0.13819660
        mat_N[3, 3] = 0.58541020
        return mat_N
    msg = "Tetrahedron with %d gauss points not implemented".format(n_gauss)
    raise NotImplementedError(msg)
