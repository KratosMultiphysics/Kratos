import sys
sys.path.append('/home/gcasas/kratos')
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from cube_mesher import *
from altair_cube_mesher import *

if __name__ == "__main__":

    # mesh dimensions
    nx = int(sys.argv[1])
    ny = int(sys.argv[2])
    nz = int(sys.argv[3])

    # name tag for the mdpa
    tag = 'Kratos'

    if len(sys.argv) > 4:
        tag = sys.argv[4]

    # Channel dimensions
    xmin = 0.0
    xmax = 4.0
    ymin = 0.0
    ymax = 4.0
    zmin = 0.0
    zmax = 4.0

    if tag == 'Altair':
        box = altair_box_data(xmin, ymin, zmin, xmax, ymax, zmax, nx, ny, nz)
    elif tag == 'Kratos':
        box = box_data(xmin, ymin, zmin, xmax, ymax, zmax, nx, ny, nz)
    else:
        raise ValueError('the provided tag: ', tag, 'must be either \'Altair\' or \'Kratos\'')
    box.x_periodic = False
    # box.y_periodic = True
    def x_scale_func(box, position):
        amplitude = 0.5
        omega = 2.
        return node_x_cos(box, position, amplitude, omega)

    def y_scale_func(box, position):
        amplitude = 0.5
        omega = 2.
        return node_y_cos(box, position, amplitude, omega)

    def z_scale_func(box, position):
        amplitude = 0.5
        omega = 2.
        return node_z_cos(box, position, amplitude, omega)

    # element type
    if len(sys.argv) > 5:
        elemtype = sys.argv[5]
    else:
        elemtype = "VMS3D"
    condtype = "WallCondition3D"

    filename = "comparative_analysis_ndiv_{0}".format(nx) + tag + "Fluid.mdpa"

    with open(filename, "w") as mdpa:
        write_header(mdpa)
        generate_nodes(mdpa,box)
        # generate_nodes(mdpa,box,z_scale=z_scale_func)
        # generate_nodes(mdpa, box, x_scale=x_scale_func, y_scale=y_scale_func, z_scale=z_scale_func)
        generate_elements(mdpa, box, elemtype)
        generate_conditions(mdpa, box, condtype)
        generate_mesh_groups(mdpa, box)
