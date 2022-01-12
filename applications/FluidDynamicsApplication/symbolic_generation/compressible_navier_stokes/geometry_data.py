import sympy
from KratosMultiphysics.sympy_fe_utilities import DefineMatrix

def GeometryDataFactory(geometry_name, ngauss = None):
    geodata_dict = {
        "triangle":      TriangleData,
        "quadrilateral": QuadrilateralData,
        "hexahedron":    TetrahedronData
    }

    if geometry_name not in geodata_dict:
        msg = "Unknown geometry type. Available geometries are:\n"
        for gname in geodata_dict.keys():
            msg += " - {}\n".format(gname)
        raise ValueError(msg)

    geodata = geodata_dict[geometry_name]

    if ngauss is None:
        return geodata()
    return geodata(ngauss)


class GeometryData:
    nnodes: int
    ndims: int
    blocksize: int
    ndofs: int
    def __init__(self, ngauss):
        self.ngauss = ngauss
        self._N = None
        self._DN = None

    def N(self):
        if self._N is None:
            self._N = self._ComputeShapeFunctions()
        return self._N

    def _ComputeShapeFunctions(self):
        raise NotImplementedError("Calling base class GeometryData's ShapeFunctions.")

    def DN(self):
        if self._DN is None:
            self._DN = self._ComputeShapeFunctionGradients()
        return self._DN

    def _ComputeShapeFunctionsGradients(self):
        raise NotImplementedError("Calling base class GeometryData's ShapeFunctionsGradients.")


class TriangleData(GeometryData):
    nnodes = 3
    ndims = 2
    blocksize = ndims+2
    ndofs = blocksize*nnodes

    def __init__(self, ngauss = 3):
        super().__init__(ngauss)

    def _ComputeShapeFunctions(self):
        mat_N = sympy.zeros(self.ngauss, 3)
        if self.ngauss == 1:
            return sympy.Matrix(self.ngauss, 3, lambda *_: 1/3)
        if self.ngauss == 3:
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
        msg = "Triangle with {} gauss points not implemented".format(self.ngauss)
        raise NotImplementedError(msg)

    def _ComputeShapeFunctionsGradients(self):
        return DefineMatrix('DN', self.nnodes, self.ndims)



class QuadrilateralData(GeometryData):
    nnodes = 4
    ndims = 2
    blocksize = ndims+2
    ndofs = blocksize*nnodes

    def __init__(self, ngauss = 4):
        super().__init__(ngauss)

    def _ComputeShapeFunctions(self):
        N = [lambda xi, eta: 0.25 * (1 - xi) * (1 - eta),
            lambda xi, eta: 0.25 * (1 + xi) * (1 - eta),
            lambda xi, eta: 0.25 * (1 - xi) * (1 + eta),
            lambda xi, eta: 0.25 * (1 + xi) * (1 + eta)]

        if self.ngauss == 1:
            return sympy.Matrix(self.ngauss, 4, lambda _, j: N[j](0, 0))
        if self.ngauss == 4:
            xi = [-3**-0.5,  3**0.5, -3**-0.5, 3**0.5]
            eta = [-3**-0.5, -3**0.5,  3**-0.5, 3**0.5]
            return sympy.Matrix(self.ngauss, 4, lambda i, j: N[j](xi[i], eta[i]))

        msg = "Quadrilateral with {} gauss points not implemented".format(self.ngauss)
        raise NotImplementedError(msg)

    def _ComputeShapeFunctionsGradients(self):
        return DefineMatrix('DN', self.nnodes, self.ndims)



class TetrahedronData(GeometryData):
    nnodes = 4
    ndims = 3
    blocksize = ndims+2
    ndofs = blocksize*nnodes

    def __init__(self, ngauss = 4):
        super().__init__(ngauss)

    def _ComputeShapeFunctions(self):
        if self.ngauss == 1:
            return sympy.Matrix(self.ngauss, self.nnodes, lambda *_: 0.25)
        if self.ngauss == 4:
            mat_N = sympy.zeros(self.ngauss,  self.nnodes)
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
        msg = "Tetrahedron with {} gauss points not implemented".format(self.ngauss)
        raise NotImplementedError(msg)

    def _ComputeShapeFunctionsGradients(self):
        return DefineMatrix('DN', self.nnodes, self.ndims)

