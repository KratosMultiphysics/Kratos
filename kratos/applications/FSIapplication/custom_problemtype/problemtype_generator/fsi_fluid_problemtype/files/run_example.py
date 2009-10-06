import fsi_fluid_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = fsi_fluid_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = fsi_fluid_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = fsi_fluid_var.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_ALEApplication = True
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_StructuralApplication = True
applications_interface.Import_FSIApplication = True
applications_interface.Import_ExternalSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosStructuralApplication import *
from KratosIncompressibleFluidApplication import *
from KratosExternalSolversApplication import *

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  
structure_model_part = ModelPart("StructurePart");  

#############################################
##importing the solvers needed
import mesh_solver
import incompressible_fluid_solver
import NonConformant_OneSideMap
import Conformant_OneSideMap
import structural_solver_dynamic
import FractionalStepCouplingAitken

incompressible_fluid_solver.AddVariables(fluid_model_part)
mesh_solver.AddVariables(fluid_model_part)
if(domain_size == 3):
    NonConformant_OneSideMap.AddVariables(fluid_model_part,structure_model_part)
else:
    Conformant_OneSideMap.AddVariables(fluid_model_part,structure_model_part)
structural_solver_dynamic.AddVariables(structure_model_part)
FractionalStepCouplingAitken.AddVariables(fluid_model_part,structure_model_part)
fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
fluid_model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
structure_model_part.AddNodalSolutionStepVariable(FORCE)
structure_model_part.AddNodalSolutionStepVariable(ACCELERATION)
structure_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
structure_model_part.AddNodalSolutionStepVariable(NORMAL)

#introducing input file name
input_file_name = fsi_fluid_var.problem_name
input_file_name_structure = fsi_fluid_var.structure_file

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#reading the structural part
model_part_io_structure = ModelPartIO(input_file_name_structure)
model_part_io_structure.ReadModelPart(structure_model_part)
print structure_model_part
print "structural model read correctly"

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)
structure_model_part.SetBufferSize(2)

##adding dofs
incompressible_fluid_solver.AddDofs(fluid_model_part)
mesh_solver.AddDofs(fluid_model_part)
structural_solver_dynamic.AddDofs(structure_model_part)


#########select here the laplacian form!!!!!!!!!!!!!!!!!
laplacian_form = fsi_fluid_var.laplacian_form 
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

##check to ensure that no node has zero density or pressure
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise 'node with zero density found'

#creating the solvers
#fluid solver
fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
fluid_solver.laplacian_form = laplacian_form; #standard laplacian form
fluid_solver.predictor_corrector = False
fluid_solver.max_press_its = 10
fluid_solver.Initialize()
print "fluid solver created"

#mesh solver
reform_dofs_at_each_step = False
mesh_solver = mesh_solver.MeshSolver(fluid_model_part,domain_size,reform_dofs_at_each_step)
pDiagPrecond = DiagonalPreconditioner()
mesh_solver.linear_solver = CGSolver(1e-9, 1000, pDiagPrecond)
mesh_solver.time_order = 2
mesh_solver.Initialize()
mesh_solver.solver.SetEchoLevel(0);
print "mesh solver created"

#structure solver
structure_solver = structural_solver_dynamic.DynamicStructuralSolver(structure_model_part,domain_size)
structure_solver.echo_level = 0
structure_solver.toll = 1e-12
structure_solver.absolute_tol = 1e-15
structure_solver.structure_linear_solver = SuperLUSolver()
structure_solver.Initialize()
if(domain_size == 3):
    structure_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
else:
    structure_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
print "structural solver created"


#mapper
#non conformant mapper
print fluid_model_part
print structure_model_part
if(domain_size == 3):
    mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(fluid_model_part,structure_model_part)
else:
    #conformant point to point
    utilities = VariableUtils()
    interface_fluid_nodes = (utilities).SelectNodeList(IS_INTERFACE,1.0,fluid_model_part.Nodes)
    interface_structure_nodes = (utilities).SelectNodeList(IS_INTERFACE,1.0,structure_model_part.Nodes)
    print "interface fluid nodes = ",len(interface_fluid_nodes)
    print "interface structure nodes = ",len(interface_structure_nodes)
    mapper = Conformant_OneSideMap.Conformant_OneSideMap(interface_fluid_nodes,interface_structure_nodes)

print "mapper created"



#creating the coupled solver
coupled_solver = FractionalStepCouplingAitken.FractionalStepCoupling(fluid_model_part,structure_model_part,structure_solver,mesh_solver,mapper,domain_size)
coupled_solver.pressure_linear_solver = SuperLUSolver()
coupled_solver.fsi_absolute_toll = 1e-9
coupled_solver.fsi_convergence_toll = 0.0001; #0.001; #pressure tolerance
coupled_solver.max_coupled_its = 50
coupled_solver.laplacian_form = laplacian_form
coupled_solver.complete_mesh_move_at_iterations = False
coupled_solver.force_prediction_order = 0 #if 0 no pressure prediction ... using the value available
coupled_solver.incremental_structural_solution = True
coupled_solver.switch_off_accelerator = True
coupled_solver.Initialize()


print "coupled solver created"


it_count = open("iteration_count.csv", 'w')
it_count.write( "count of iterations at the various time steps" )
 


#settings to be changed
nsteps = fsi_fluid_var.nsteps 
output_step = 1

Dt = fsi_fluid_var.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01

out = 0

#mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( fluid_model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())


time = 0.0
for step in range(0,nsteps):

    if(step < 5):
        Dt = initial_Dt
    else:
        Dt = full_Dt
        
    time = time + Dt
    fluid_model_part.CloneTimeStep(time)
    structure_model_part.CloneTimeStep(time)

    if(step >= 3):
        its = coupled_solver.Solve()

        it_count.write(str(time) + " " + str(its) + "\n")
        it_count.flush()

    if(out == output_step):
        gid_io.WriteNodalResults(MESH_VELOCITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISPLACEMENT,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_INTERFACE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(EXTERNAL_PRESSURE,fluid_model_part.Nodes,time,0)

        out = 0

    out = out + 1

gid_io.FinalizeResults()
          
        

