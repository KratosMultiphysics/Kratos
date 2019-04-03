from KratosMultiphysics.FluidDynamicsApplication import *
import cube_mesher

BaseClass = cube_mesher.box_data

class altair_box_data(BaseClass):
    def __init__(self, xmin, ymin, zmin, xmax, ymax, zmax, nx, ny, nz):
        BaseClass.__init__(self, xmin, ymin, zmin, xmax, ymax, zmax, nx, ny, nz)

    def cube_vertices(self, ix, iy, iz):
        """ Identify 8 contiguous nodes forming a cube.
            Note: Even and odd levels are rotated 90 degrees along Z
            to ensure that they are conformant.
        """

        # node0 = (n+1)*(n+1)*k+(n+1)*j+i+1
        # node1 = (n+1)*(n+1)*(k+1)+(n+1)*j+i+1
        # node2 = (n+1)*(n+1)*(k+1)+(n+1)*j+i+2
        # node3 = (n+1)*(n+1)*k+(n+1)*j+i+2
        # node4 = (n+1)*(n+1)*k+(n+1)*(j+1)+i+1
        # node5 = (n+1)*(n+1)*(k+1)+(n+1)*(j+1)+i+1
        # node6 = (n+1)*(n+1)*(k+1)+(n+1)*(j+1)+i+2
        # node7 = (n+1)*(n+1)*k+(n+1)*(j+1)+i+2

        n0 = self.get_id(ix, iy, iz)
        n1 = self.get_id(ix, iy, iz + 1)
        n2 = self.get_id(ix + 1, iy, iz + 1)
        n3 = self.get_id(ix + 1, iy, iz)

        n4 = self.get_id(ix, iy + 1, iz)
        n5 = self.get_id(ix, iy + 1, iz + 1)
        n6 = self.get_id(ix + 1, iy + 1, iz + 1)
        n7 = self.get_id(ix + 1, iy + 1, iz)

        return n0, n1, n2, n3, n4, n5, n6, n7

    def generate_elements(mdpa, box, elemtype="FractionalStep3D", prop_id=0):
        index = 1
        nx = box.nx()
        ny = box.ny()
        nz = box.nz()

        print("Generating tetrahedral {0} elements.".format(elemtype))
        mdpa.write("Begin Elements {0}\n".format(elemtype))
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    n0, n1, n2, n3, n4, n5, n6, n7 = box.cube_vertices(ix, iy, iz)

                    # fill the cube with 6 tetrahedra
                    mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index, prop_id, n3, n2, n6, n0))
                    mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index + 1, prop_id, n0, n3, n7, n6))
                    mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index + 2, prop_id, n1, n5, n6, n0))
                    mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index + 3, prop_id, n0, n4, n5, n6))
                    mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index + 4, prop_id, n0, n1, n2, n6))
                    mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index + 5, prop_id, n4, n7, n6, n0))
                    index += 6

        mdpa.write("End Elements\n\n")