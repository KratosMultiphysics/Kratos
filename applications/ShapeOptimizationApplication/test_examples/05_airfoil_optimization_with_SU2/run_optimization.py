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
interface_su2.WriteSU2MeshAsMDPA()

# read mdpa from su2
mdpa_from_su2 = ModelPart("mesh_NACA0012_inv")
mdpa_from_su2.ProcessInfo.SetValue(DOMAIN_SIZE, 2)

model_part_io = ModelPartIO(("mesh_NACA0012_inv"))
model_part_io.ReadModelPart(mdpa_from_su2)

# Work with mdpa from su2
interface_su2.WriteNodesAsSU2MeshMotionFile(mdpa_from_su2.GetSubModelPart("airfoil").GetNodes())

