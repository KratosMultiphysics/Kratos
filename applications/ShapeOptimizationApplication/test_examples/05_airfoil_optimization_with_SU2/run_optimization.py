# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from interface_su2 import InterfaceSU2

# Create and use SU2 translator
interface_parameters = Parameters("""
{
    "su2_related":
    {
        "config_file"        : "inv_NACA0012.cfg",
        "number_of_cores"    : 14,
        "gradient_method"    : "DISCRETE_ADJOINT",
        "mesh_file"          : "mesh_NACA0012_inv.su2",
        "mesh_motion_file"   : "mesh_motion.dat",
        "design_surface_tag" : "airfoil"
    },
    "kratos_related":
    {
        "mdpa_file"      : "mesh_NACA0012_inv.mdpa",
        "write_elements" : false
    }
}""")

interface_su2 = InterfaceSU2(interface_parameters)

print("> Start writing mdpa file...")
interface_su2.WriteSU2MeshAsMDPA()

# read mdpa from su2
print("> Start reading mdpa file...")
mdpa_from_su2 = ModelPart("mesh_NACA0012_inv")
mdpa_from_su2.ProcessInfo.SetValue(DOMAIN_SIZE, 2)

model_part_io = ModelPartIO(("mesh_NACA0012_inv"))
model_part_io.ReadModelPart(mdpa_from_su2)

# Work with su2 interface
print("\n> Initializing SU2 project...")
interface_su2.InitializeNewSU2Project()

design_surface = mdpa_from_su2.GetSubModelPart("airfoil")

print("> Start writing SU2 file for mesh motion...")
interface_su2.WriteNodesAsSU2MeshMotionFile(design_surface.GetNodes())

for iteration in range(1,4):

    print("#############################################")
    print("> Starting iteration ",iteration)
    print("#############################################")

    print("> Starting SU2 response calculations...")
    drag,lift = interface_su2.ComputeValues(["DRAG","LIFT"],True)

    print("> Start writing SU2 file for mesh motion...")
    for node in design_surface.Nodes:
        node.Y0 = node.Y0+0.1*iteration*node.X0
    interface_su2.WriteNodesAsSU2MeshMotionFile(design_surface.GetNodes(),"DESIGNS/DSN_00"+str(iteration))

    # print("> Starting SU2 gradient calculations...")
    # ddrag = interface_su2.ComputeGradients(["DRAG"],False,True)
    # ddlift = interface_su2.ComputeGradients(["LIFT"],False,True)