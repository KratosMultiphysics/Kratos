#import the configuration data as read from the GiD
import fsi_structure_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = fsi_structure_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = fsi_structure_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = fsi_structure_var.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_StructuralApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)
from KratosStructuralApplication import *

## from now on the order is not anymore crucial
##################################################################
##################################################################

#defining a model part
model_part = ModelPart("StructurePart");  

#adding of Variables to Model Part should be here when the "very fix container will be ready"
if(fsi_structure_var.SolverType == "StaticSolver"): 
    import structural_solver_static
    structural_solver_static.AddVariables(model_part)
elif(fsi_structure_var.SolverType == "DynamicSolver"):
    import structural_solver_dynamic
    structural_solver_dynamic.AddVariables(model_part)


#reading a model
name = fsi_structure_var.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io = ModelPartIO(name)
model_part_io.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()

print model_part
print model_part.Properties

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
if(fsi_structure_var.SolverType == "StaticSolver"):
    print ">>>>>>>>>>>>>>>>< 1yquqee "
    structural_solver_static.AddDofs(model_part)
    solver = structural_solver_static.StaticStructuralSolver(model_part,domain_size)
elif(fsi_structure_var.SolverType == "DynamicSolver"):
    print " *       ASRDF +********************************  FF GDFSG DG G"
    structural_solver_dynamic.AddDofs(model_part)
    solver = structural_solver_dynamic.DynamicStructuralSolver(model_part,domain_size)

#creating a fluid solver object
model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
model_part.Properties[2].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
print "Linear elastic model selected"

solver.Initialize()
(solver.solver).SetEchoLevel(2);

Dt = fsi_structure_var.Dt
nsteps = 10

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())

for step in range(0,nsteps):
    print "line49"

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        solver.Solve()

        #print the results
        gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
#        gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time,domain_size)

gid_io.FinalizeResults()

          
        

