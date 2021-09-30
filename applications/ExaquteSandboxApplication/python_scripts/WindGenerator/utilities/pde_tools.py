

 
#=====================================================================================================
#=====================================================================================================
#
#                                           PDE FRAMEWORK
#
#=====================================================================================================
#=====================================================================================================

import dolfin as dl

#######################################################################################################
#	Mesh
#######################################################################################################

def generate_Mesh(shape):
    ndim = len(shape)
    Nd = tuple(shape)
    if ndim == 1:
        mesh = dl.UnitIntervalMesh(shape[0])
    elif ndim == 2:
        mesh = dl.UnitSquareMesh.create(shape[0], shape[1], dl.CellType.Type.quadrilateral)
        # mesh = dl.UnitSquareMesh(shape[0], shape[1])
    elif ndim == 3:
        mesh = dl.UnitCubeMesh.create(shape[0], shape[1], shape[2], dl.CellType.Type.hexahedron)
    else:
        raise Exception("The case of Dimension={0:d} is not inplemented!".format(ndim))
    return mesh


#######################################################################################################
#	Periodic boundary conditions
#######################################################################################################

class PeriodicBoundary(dl.SubDomain):

	def __init__(self, ndim):
		dl.SubDomain.__init__(self)
		self.ndim = ndim

	def inside(self, x, on_boundary):
		return bool( 	any([ dl.near(x[j], 0) for j in range(self.ndim) ]) and 
					not any([ dl.near(x[j], 1) for j in range(self.ndim) ]) and 
						on_boundary
					)

	def map(self, x, y):
		for j in range(self.ndim):
			if dl.near(x[j], 1):
				y[j] = x[j] - 1
			else:
				y[j] = x[j]

#######################################################################################################