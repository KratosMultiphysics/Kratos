##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs/' ##kratos_root/libs
kratos_applications_path = '../../../../applications/' ##kratos_root/applications
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
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosStructuralApplication import *



#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  
structure_model_part = ModelPart("StructurePart");  

#adding variables for all of the tools needed
import incompressible_fluid_solver
incompressible_fluid_solver.AddVariables(fluid_model_part)
import mesh_solver
mesh_solver.AddVariables(fluid_model_part)
import Conformant_OneSideMap
Conformant_OneSideMap.AddVariables(fluid_model_part,structure_model_part)
import structural_solver_dynamic
structural_solver_dynamic.AddVariables(structure_model_part)
import ExplicitCoupling
#NonConformant_OneSideMap.AddVariables(fluid_model_part,structure_model_part)


#introducing input file name
input_file_name = "flag_test"

#reading the fluid part
write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
write_elements = WriteConditionsFlag.WriteElementsOnly
post_mode = GiDPostMode.GiD_PostBinary
multi_file_flag = MultiFileFlag.MultipleFiles
gid_io = GidIO( input_file_name+str("_fluid"), post_mode, multi_file_flag, write_deformed_flag, write_elements )
gid_io.ReadModelPart(fluid_model_part)
print fluid_model_part
print "fluid model read correctly"

#reading the structural part
data_io = DatafileIO(input_file_name+str("_structure"))
data_io.ReadModelPart(structure_model_part)
print structure_model_part
print "structural model read correctly"

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)
structure_model_part.SetBufferSize(2)

##adding dofs
incompressible_fluid_solver.AddDofs(fluid_model_part)
mesh_solver.AddDofs(fluid_model_part)
structural_solver_dynamic.AddDofs(structure_model_part)

#assigning the fluid properties
density = 1.18;
viscosity = 1.82e-5 / density
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,viscosity);
    node.SetSolutionStepValue(DENSITY,0,density);


#creating the solvers
#fluid solver
fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
fluid_solver.laplacian_form = 1; #standard laplacian form
fluid_solver.predictor_corrector = False
fluid_solver.max_press_its = 10
fluid_solver.Initialize()
print "fluid solver created"

#mesh solver
reform_dofs_at_each_step = False
mesh_solver = mesh_solver.MeshSolver(fluid_model_part,domain_size,reform_dofs_at_each_step)
pDiagPrecond = DiagonalPreconditioner()
mesh_solver.linear_solver = CGSolver(1e-3, 300, pDiagPrecond)
mesh_solver.time_order = 2
mesh_solver.Initialize()
mesh_solver.solver.SetEchoLevel(0);
print "mesh solver created"

#structure solver
structure_solver = structural_solver_dynamic.DynamicStructuralSolver(structure_model_part,domain_size)
#structure_solver.echo_level = 0
#pILUPrecond = ILU0Preconditioner()
#structure_solver.structure_linear_solver = BICGSTABSolver(1e-8, 5000,pILUPrecond)
structure_solver.Initialize()
structure_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
print "structural solver created"

#mapper
#non conformant mapper
#mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(fluid_model_part,structure_model_part)
#conformant point to point
utilities = VariableUtils()
interface_fluid_nodes = (utilities).SelectNodeList(IS_INTERFACE,1.0,fluid_model_part.Nodes)
interface_structure_nodes = (utilities).SelectNodeList(IS_INTERFACE,1.0,structure_model_part.Nodes)
print "interface fluid nodes = ",len(interface_fluid_nodes)
print "interface structure nodes = ",len(interface_structure_nodes)
mapper = Conformant_OneSideMap.Conformant_OneSideMap(interface_fluid_nodes,interface_structure_nodes)

print "mapper created"

#creating the coupled solver
coupled_solver = ExplicitCoupling.ExplicitCoupling(fluid_model_part,structure_model_part,structure_solver,fluid_solver,mesh_solver,mapper,domain_size)
print "coupled solver created"


    



 


#settings to be changed
nsteps = 5000
output_step = 10

Dt = 0.005

out = 0

#mesh to be printed
gid_io.InitializeMesh(0.0)
gid_io.WriteMesh((fluid_model_part).GetMesh());
##gid_io.WriteMesh((structure_model_part).GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);
gid_io.FinalizeMesh()

for step in range(0,nsteps):

    time = Dt*step
    fluid_model_part.CloneTimeStep(time)
    structure_model_part.CloneTimeStep(time)

    print "time = " , time

    #solving the fluid problem
    if(step < 5 ):
        if(step > 3):
            print "solving only the fluid - starting procedure"
            fluid_solver.Solve()
            print "solution -only fluid- complete"
    else:
        print "coupled solve"
        coupled_solver.Solve()

    #print the results
##    gid_io.WriteNodalResults(DISPLACEMENT,structure_model_part.Nodes,time,0)
##    gid_io.WriteNodalResults(POSITIVE_FACE_PRESSURE,structure_model_part.Nodes,time,0)

    if(out == output_step):
        gid_io.InitializeResults( time, (fluid_model_part).GetMesh() )
#        gid_io.WriteNodalResults(MESH_VELOCITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISPLACEMENT,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
#        gid_io.WriteNodalResults(IS_INTERFACE,fluid_model_part.Nodes,time,0)

##    gid_io.WriteNodalResults(IS_INTERFACE,structure_model_part.Nodes,time,0)
##    gid_io.WriteNodalResults(POSITIVE_FACE_PRESSURE,structure_model_part.Nodes,time,0)
##    gid_io.WriteNodalResults(NEGATIVE_FACE_PRESSURE,structure_model_part.Nodes,time,0)
##    gid_io.WriteNodalResults(DISPLACEMENT,structure_model_part.Nodes,time,0)

        #gid_io.Flush();
        out = 0

    out = out + 1

          
        

