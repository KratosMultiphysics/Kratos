import kElectrostatic

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = kElectrostatic.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = kElectrostatic.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path = kElectrostatic.kratos_path + '/applications' ##kratos_root/applications
#kratos_benchmarking_path = kElectrostatic.kratos_path +'/benchmarking' ##kratos_root/benchmarking
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
#sys.path.append(kratos_benchmarking_path)

#import benchmarking

print "Before importing kratos..."

#importing Kratos main library
from Kratos import **
print "Kratos library imported"
kernel = Kernel()   #defining kernel
print "kernel created"
#importing applications
import applications_interface
applications_interface.Import_ElectrostaticApplication = True
#applications_interface.Import_ExternalSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosR1ElectrostaticApplication import **
#from KratosExternalSolversApplication import **

##x = raw_input("stopped to allow debug: set breakpoints and press enter to continue");

#defining a model part for the electrostatic
electrostatic_model_part = ModelPart("ASABERQUEElectrostaticPart");  
#model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);

SolverType = kElectrostatic.SolverType

#adding of Variables to Model Part
#import nonlinear_convection_diffusion_solver
import static_poisson_solver
static_poisson_solver.AddVariables(electrostatic_model_part)

#introducing input file name
input_file_name = kElectrostatic.problem_name

print "Before reading from GiD... " + input_file_name

#reading the electrostatic part
#gid_mode = GiDPostMode.GiD_PostAscii
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
#gid_io = ElectrostaticGidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag,write_conditions)
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag,write_conditions)
#gid_io.ReadModelPart(model_part)
model_part_io_electrostatic = ModelPartIO(input_file_name)
model_part_io_electrostatic.ReadModelPart(electrostatic_model_part)

print "After reading from GiD..."

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
electrostatic_model_part.SetBufferSize(3)

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((electrostatic_model_part).GetMesh());
gid_io.FinalizeMesh()
print electrostatic_model_part

#importing the solver files
#nonlinear_convection_diffusion_solver.AddDofs(model_part)
static_poisson_solver.AddDofs(electrostatic_model_part)
for node in electrostatic_model_part.Nodes:
    #adding dofs
    node.AddDof(ELECTROSTATIC_POTENTIAL);
print "variables for the electromagnetic solver added correctly"
    
#creating a electrostatic solver object
electrostatic_solver = static_poisson_solver.StaticPoissonSolver(electrostatic_model_part,domain_size)
electrostatic_solver.time_order = 1
electrostatic_solver.linear_solver = SkylineLUFactorizationSolver();
electrostatic_solver.echo_level = 0
electrostatic_solver.Initialize()

print "static poisson solver created"

gid_io.InitializeResults(mesh_name,(electrostatic_model_part).GetMesh())

print "GiD IO initialized"

electrostatic_solver.Solve()

print "solver solved"

gid_io.WriteNodalResults(ELECTROSTATIC_POTENTIAL,electrostatic_model_part.Nodes,0,0)
gid_io.PrintOnGaussPoints(ELECTRIC_FIELD, electrostatic_model_part, 0)
gid_io.PrintOnGaussPoints(ELECTRIC_DISPLACEMENT_FIELD, electrostatic_model_part, 0)
gid_io.PrintOnGaussPoints(ELECTRICAL_PERMITTIVITY, electrostatic_model_part, 0)
gid_io.FinalizeResults()

