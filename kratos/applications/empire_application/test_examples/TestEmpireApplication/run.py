import sys

if len(sys.argv) < 2:
    print "ERROR: An input filename is missing"
    print ""

from KratosMultiphysics import *
from KratosMultiphysics.EmpireApplication import *
from KratosMultiphysics.StructuralApplication import *

model_part = ModelPart("ModelPart")

model_part_io = ModelPartIO(str(sys.argv[1]))
model_part_io.ReadModelPart(model_part)

print "######## MODEL ########"
print model_part
print "#######################"

# Creating GidIO
gid_mode = GiDPostMode.GiD_PostBinary    # or GiDPostMode.GiD_PostAscii
use_multi_file = MultiFileFlag.MultipleFiles    # or MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed    # or WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly   # or WriteConditionsFlag.WriteConditions
gid_io = GidIO(sys.argv[1],gid_mode,use_multi_file,deformed_mesh_flag, write_conditions)

test_process = TestProcess(model_part)
test_process.TestMemberFunction()
