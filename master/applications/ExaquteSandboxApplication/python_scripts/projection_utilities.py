from math import floor, ceil
from collections import namedtuple

def isclose(a, b, rel_tol=1e-9, abs_tol=0.):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def _DotProduct(A,B):
    return sum(i[0]*i[1] for i in zip(A, B))

Extent = namedtuple('Extent', ['lower', 'upper'])
IndexSpan = namedtuple('IndexSpan', ['begin', 'end'])

class RegularGrid1D:

    @property
    def lower_bound(self):
        return self.extent.lower

    @property
    def upper_bound(self):
        return self.extent.upper

    def __init__(self, start_pos, length, size):
        self.extent = Extent(start_pos, start_pos+length)
        self.size = size
        self.step_size = (self.upper_bound-self.lower_bound) / (self.size-1)

    def __getitem__(self, index):
        return self.lower_bound + self.step_size*index

    def __len__(self):
        return self.size

    def floor_index(self, coord):
        if isclose(coord, self.lower_bound, abs_tol=1e-4):
            # Guard against negative index due to floating point representation.
            return 0
        local_coord = coord - self.lower_bound
        return int(floor(local_coord / self.step_size))

    def ceil_index(self, coord):
        if isclose(coord, self.upper_bound, abs_tol=1e-4):
            return len(self) - 1
        local_coord = coord - self.lower_bound
        return int(ceil(local_coord / self.step_size))

    def has(self, coord):
        return (self.floor_index(coord) >= 0
                and self.ceil_index(coord) < len(self))


class DomainPanel3D:

    def __init__(self, grid, data):
        self.grid = grid
        self.data = data
        self.dx = self.grid.x.step_size
        self.x0 = self.grid.x.lower_bound
        self.dy = self.grid.y.step_size
        self.y0 = self.grid.y.lower_bound
        self.dz = self.grid.z.step_size
        self.z0 = self.grid.z.lower_bound

    def interpolate(self, node):
        # xi and gamma and eta are local mesh cell coordinates in the interval [-1, 1].
        xi = ((node.X - self.x0) % self.dx) / self.dx
        gamma = ((node.Y - self.y0) % self.dy) / self.dy
        eta = ((node.Z - self.z0) % self.dz) / self.dz
        xi = 2.0 * (xi-0.5)
        gamma = 2.0 * (gamma-0.5)
        eta = 2.0 * (eta-0.5)
        # Interpolate using trilinear shape functions.
        weights = (
            0.125 * (1.0-xi) * (1.0-gamma) * (1.0-eta), # - - -
            0.125 * (1.0+xi) * (1.0-gamma) * (1.0-eta), # + - -
            0.125 * (1.0+xi) * (1.0+gamma) * (1.0-eta), # + + -
            0.125 * (1.0+xi) * (1.0+gamma) * (1.0+eta), # + + +
            0.125 * (1.0+xi) * (1.0-gamma) * (1.0+eta), # + - +
            0.125 * (1.0-xi) * (1.0-gamma) * (1.0+eta), # - - +
            0.125 * (1.0-xi) * (1.0+gamma) * (1.0-eta), # - + -
            0.125 * (1.0-xi) * (1.0+gamma) * (1.0+eta)  # - + +
        )
        i = self.grid.x.floor_index(node.X)
        j = self.grid.y.floor_index(node.Y)
        k = self.grid.z.floor_index(node.Z)
        return (
            weights[0] * self.data[i,j,k]
            + weights[1] * self.data[i+1,j,k]
            + weights[2] * self.data[i+1,j+1,k]
            + weights[3] * self.data[i+1,j+1,k+1]
            + weights[4] * self.data[i+1,j,k+1]
            + weights[5] * self.data[i,j,k+1]
            + weights[6] * self.data[i,j+1,k]
            + weights[7] * self.data[i,j+1,k+1]
        )


class DomainPanel2D:

    def __init__(self, grid, data):
        self.grid = grid
        self.data = data
        self.dx = self.grid.x.step_size
        self.x0 = self.grid.x.lower_bound
        self.dz = self.grid.z.step_size
        self.z0 = self.grid.z.lower_bound

    def interpolate(self, node):
        # xi and eta are local mesh cell coordinates in the interval [-1, 1].
        xi = ((node.X - self.x0) % self.dx) / self.dx
        eta = ((node.Y - self.z0) % self.dz) / self.dz
        xi = 2.0 * (xi-0.5)
        eta = 2.0 * (eta-0.5)
        # Interpolate using bilinear shape functions.
        weights = (
            0.25 * (1.0-xi) * (1.0-eta), # - -
            0.25 * (1.0+xi) * (1.0-eta), # + -
            0.25 * (1.0+xi) * (1.0+eta), # + +
            0.25 * (1.0-xi) * (1.0+eta)  # - +
        )
        j = self.grid.x.floor_index(node.X)
        k = self.grid.z.floor_index(node.Y)
        return (
            weights[0] * self.data[j, k]
            + weights[1] * self.data[j+1, k]
            + weights[2] * self.data[j+1, k+1]
            + weights[3] * self.data[j, k+1]
        )
