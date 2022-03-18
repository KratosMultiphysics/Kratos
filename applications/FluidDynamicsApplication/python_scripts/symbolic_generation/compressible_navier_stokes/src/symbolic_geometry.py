import sympy
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src.defines \
    import CompressibleNavierStokesDefines as defs


def GeometryDataFactory(geometry_name, ngauss=None):
    geodata_dict = {
        "triangle":      TriangleData,
        "quadrilateral": QuadrilateralData,
        "tetrahedron":   TetrahedronData
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
    nnodes = None                   # Number of nodes
    ndims = None                    # Dimension the geometry is embeded in
    blocksize = None                # Number of Dof per node. Generally = ndims+2
    ndofs = None                    # Number of Dof per element. Generally = blocksize*nnodes

    symbolic_integration = None
    # True if the integration loop is to be performed in the generator (only valid for simplex geometries)
    # False if the integration loop is in the template

    name = None

    def __init__(self, ngauss):
        self.ngauss = ngauss
        self._N = None
        self._DN = None

    def __str__(self):
        return self.name

    def SymbolicIntegrationPoints(self):
        "Returns the number of gauss points to be evaluated at symbolic time"
        return range(self.ngauss) if self.symbolic_integration else [0]

    def N(self):
        "Shape functions evaluated at all gauss points"
        if self._N is None:
            self._N = self._ComputeShapeFunctions()
        return self._N

    def N_gauss(self, i_gauss):
        "Returns the shape functions at a gauss point as a vertical vector"
        if self.symbolic_integration:
            return sympy.Matrix(self.N()[i_gauss, :]).T
        else:
            return self.N()

    def _ComputeShapeFunctions(self):
        raise NotImplementedError("Calling base class GeometryData's _ComputeShapeFunctions.")

    def DN(self):
        "Shape functions gradients evaluated at all gauss points"
        if self._DN is None:
            self._DN = self._ComputeShapeFunctionsGradients()
        return self._DN

    def _ComputeShapeFunctionsGradients(self):
        raise NotImplementedError("Calling base class GeometryData's _ComputeShapeFunctionsGradients.")


class TriangleData(GeometryData):
    nnodes = 3
    ndims = 2
    blocksize = ndims+2
    ndofs = blocksize*nnodes
    symbolic_integration = True
    name = "triangle (2D3N)"

    def __init__(self, ngauss=3):
        super().__init__(ngauss)

    def _ComputeShapeFunctions(self):
        mat_N = defs.ZeroMatrix(self.ngauss, self.nnodes)
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
        return defs.Matrix('DN_DX', self.nnodes, self.ndims)


class QuadrilateralData(GeometryData):
    nnodes = 4
    ndims = 2
    blocksize = ndims+2
    ndofs = blocksize*nnodes
    symbolic_integration = False
    name = "quadrilateral (2D4N)"

    def __init__(self, ngauss=4):
        super().__init__(ngauss)

    def _ComputeShapeFunctions(self):
        return defs.Vector('N', self.nnodes, positive=True)

    def _ComputeShapeFunctionsGradients(self):
        return defs.Matrix('DN_DX', self.nnodes, self.ndims)


class TetrahedronData(GeometryData):
    nnodes = 4
    ndims = 3
    blocksize = ndims+2
    ndofs = blocksize*nnodes
    symbolic_integration = True
    name = "tetrahedron (3D4N)"

    def __init__(self, ngauss=4):
        super().__init__(ngauss)

    def _ComputeShapeFunctions(self):
        if self.ngauss == 1:
            return sympy.Matrix(self.ngauss, self.nnodes, lambda *_: 0.25)
        if self.ngauss == 4:
            mat_N = defs.ZeroMatrix(self.ngauss,  self.nnodes)
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
        return defs.Matrix('DN_DX', self.nnodes, self.ndims)
