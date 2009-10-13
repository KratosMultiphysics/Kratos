##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

print "before importing kratos"

#importing Kratos main library
from Kratos import *
print "Kratos library imported"
kernel = Kernel()   #defining kernel
print "kernel created"
#importing applications
import applications_interface
applications_interface.Import_ElectrostaticApplication = True

applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosR1ElectrostaticApplication import *

##x = raw_input("stopped to allow debug: set breakpoints and press enter to continue");

#defining a model part
model_part = ModelPart("ASABERQUEElectrostaticPart");  
#model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);

#adding of Variables to Model Part
#import nonlinear_convection_diffusion_solver
import static_poisson_solver
static_poisson_solver.AddVariables(model_part)

model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
model_part.AddNodalSolutionStepVariable(ELECTRICAL_PERMITTIVITY);
model_part.AddNodalSolutionStepVariable(ELECTRIC_FIELD);
model_part.AddNodalSolutionStepVariable(ELECTRIC_DISPLACEMENT_FIELD);

model_part.AddNodalSolutionStepVariable(ELECTROSTATIC_POTENTIAL);
model_part.AddNodalSolutionStepVariable(ELECTROSTATIC_POINT_CHARGE);

#reading a model
#gid_mode = GiDPostMode.GiD_PostAscii
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("cubito2",gid_mode,multifile,deformed_mesh_flag,write_conditions)
gid_io.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#importing the solver files
#nonlinear_convection_diffusion_solver.AddDofs(model_part)
static_poisson_solver.AddDofs(model_part)
for node in model_part.Nodes:
    #adding dofs
    node.AddDof(ELECTROSTATIC_POTENTIAL);
print "variables for the electromagnetic solver added correctly"
    
#creating a fluid solver object
solver = static_poisson_solver.StaticPoissonSolver(model_part,domain_size)
solver.time_order = 1
solver.linear_solver = SkylineLUFactorizationSolver();
solver.echo_level = 0
solver.Initialize()

print "static poisson solver created"

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())

print "GiD IO initialized"

solver.Solve()

print "solver solved"

gid_io.WriteNodalResults(ELECTROSTATIC_POTENTIAL,model_part.Nodes,0,0)
gid_io.PrintOnGaussPoints(ELECTRIC_FIELD, model_part, 0)
gid_io.PrintOnGaussPoints(ELECTRIC_DISPLACEMENT_FIELD, model_part, 0)
gid_io.PrintOnGaussPoints(ELECTRICAL_PERMITTIVITY, model_part, 0)
gid_io.FinalizeResults()

