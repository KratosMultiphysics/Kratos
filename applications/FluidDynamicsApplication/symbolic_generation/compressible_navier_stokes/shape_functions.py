import sympy


def TriangleShapeFunctions(n_gauss):
    mat_N = sympy.Matrix(n_gauss, 3, lambda *_: 0)
    if n_gauss == 1:
        mat_N[0, 0] = 1.0 / 3.0
        mat_N[0, 1] = 1.0 / 3.0
        mat_N[0, 2] = 1.0 / 3.0
        return mat_N
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
    msg = "Invalid quadrature: tetrahedron with %d integration points".format(n_gauss)
    raise NotImplementedError(msg)


def QuadrilateralShapeFunctions(n_gauss):
    if n_gauss == 1:
        return sympy.Matrix(n_gauss, 4, lambda *_: 0.25)
    if n_gauss == 4:
        mat_N = sympy.Matrix(n_gauss, 4, lambda *_: 0)
        mat_N[0, 0] = 0 # TODO
        mat_N[0, 1] = 0 # TODO
        mat_N[0, 2] = 0 # TODO
        mat_N[0, 3] = 0 # TODO
        mat_N[1, 0] = 0 # TODO
        mat_N[1, 1] = 0 # TODO
        mat_N[1, 2] = 0 # TODO
        mat_N[1, 3] = 0 # TODO
        mat_N[2, 0] = 0 # TODO
        mat_N[2, 1] = 0 # TODO
        mat_N[2, 2] = 0 # TODO
        mat_N[2, 3] = 0 # TODO
        mat_N[3, 0] = 0 # TODO
        mat_N[3, 1] = 0 # TODO
        mat_N[3, 2] = 0 # TODO
        mat_N[3, 3] = 0 # TODO
        return mat_N
    msg = "Invalid quadrature: tetrahedron with %d integration points".format(n_gauss)
    raise NotImplementedError(msg)    


def TetrahedronShapeFunctions(n_gauss):
    mat_N = sympy.Matrix(n_gauss,  4, lambda *_: 0)
    if n_gauss == 1:
        mat_N[0, 0] = 1.0 / 4.0
        mat_N[0, 1] = 1.0 / 4.0
        mat_N[0, 2] = 1.0 / 4.0
        mat_N[0, 3] = 1.0 / 4.0
        return mat_N
    if n_gauss == 4:
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
    msg = "Invalid quadrature: tetrahedron with %d integration points".format(n_gauss)
    raise NotImplementedError(msg)


def DefineShapeFunctionsMatrix(dim, n_nodes, n_gauss):
    quadrature_data = {
        (2, 3): TriangleShapeFunctions,
        (2, 4): QuadrilateralShapeFunctions,
        (3, 4): TetrahedronShapeFunctions
    }
    return quadrature_data[(dim, n_nodes)](n_gauss)
