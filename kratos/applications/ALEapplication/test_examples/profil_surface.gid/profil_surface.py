#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *


#############################################
#defining the domain size
domain_size = 3
 
#############################################
#############################################
#attention!! order is crucial here!!!
#importing the applications library and adding them to the kernel
ale_app = KratosALEApplication()
kernel.AddApplication(ale_app)

#dynamic renumbering of variables to ensure consistency 
kernel.Initialize()
kernel.InitializeApplication(ale_app);
#############################################
#############################################
 
#defining a model part
fluid_model_part = ModelPart("FluidPart");  
structure_model_part = ModelPart("StructurePart");  

#importing the solver files
import incompressible_fluid_solver
import mesh_solver

#adding variables
incompressible_fluid_solver.AddVariables(model_part)
mesh_solver.AddVariables(model_part)

#reading fluid
data_io = DatafileIO("profil_surface_fluid")
data_io.ReadMesh(fluid_model_part.GetMesh())

#reading structure 
gid_io = GidIO("profil_surface_structure",GiDPostMode.GiD_PostBinary)
gid_io.ReadMesh(structure_model_part.GetMesh())
print fluid_model_part

#the buffer size should be set up here after the mesh is read for the first time
fluid_model_part.SetBufferSize(2)
structure_model_part.SetBufferSize(2)

#adding variables
fluid_model_part.AddNodalSolutionStepVariable(NODAL_AREA)
fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)
fluid_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
fluid_model_part.AddNodalSolutionStepVariable(AUX)
fluid_model_part.AddNodalSolutionStepVariable(VAUX)
fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)

for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(IS_INTERFACE,1.0)

structure_model_part.AddNodalSolutionStepVariable(NODAL_AREA)
structure_model_part.AddNodalSolutionStepVariable(PRESSURE)
structure_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
structure_model_part.AddNodalSolutionStepVariable(AUX)
structure_model_part.AddNodalSolutionStepVariable(VAUX)
structure_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)

for node in structure_model_part.Nodes:
    node.SetSolutionStepValue(IS_INTERFACE,1.0)
    
#applying a function to the fluid
for node in fluid_model_part.Nodes:
    value = node.X**2 + node.Y**2 + node.Z**2
    value = value/1000.0
    node.SetSolutionStepValue(PRESSURE,0,value)

for node in structure_model_part.Nodes:
    disp = Vector(3)
    disp[0] = node.X
    disp[1] = node.Y
    disp[2] = node.Z
    node.SetSolutionStepValue(DISPLACEMENT,0,disp)

#creating the mapper
#mapper = NMPointsMapper(fluid_model_part,structure_model_part)
fluid_to_structure_mapper = AdvancedNMPointsMapper(fluid_model_part,structure_model_part)
structure_to_fluid_mapper = AdvancedNMPointsMapper(structure_model_part,fluid_model_part)


print structure_model_part
print fluid_model_part

print "before neighbour search"
time = 0.0                      
search_radius_factor = 3.0
fluid_to_structure_mapper.FindNeighbours(search_radius_factor)
structure_to_fluid_mapper.FindNeighbours(search_radius_factor)

print "before mapping"
it_max = 10
tol = 1e-3
#fluid_to_structure_mapper.ScalarMap(PRESSURE,PRESSURE,it_max,tol)
#structure_to_fluid_mapper.VectorMap(DISPLACEMENT,DISPLACEMENT,it_max,tol)

nsteps = 10
output_step = 1

Dt = 0.1

out = 0

output_model_part = structure_model_part
gid_io.WriteConditions(output_model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);

for step in range(0,nsteps):
    print "line49"

    time = Dt*step
    output_model_part.CloneTimeStep(time)
    fluid_model_part.CloneTimeStep(time)

    for node in fluid_model_part.Nodes:
        value = node.X**2 + node.Y**2 + node.Z**2
        value = value * time
        node.SetSolutionStepValue(PRESSURE,0,value)

    fluid_to_structure_mapper.ScalarMap(PRESSURE,PRESSURE,it_max,tol)

    print " output of structure        "
    gid_io.WriteNodalResults(PRESSURE, output_model_part.Nodes, time, 0);
    gid_io.WriteNodalResults(IS_INTERFACE, output_model_part.Nodes, time, 0);
    gid_io.Flush()



##print " output of structure        "
##output_model_part = structure_model_part
##gid_io.WriteConditions(output_model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);
##gid_io.WriteNodalResults(PRESSURE, output_model_part.Nodes, time, 0);
##gid_io.WriteNodalResults(IS_INTERFACE, output_model_part.Nodes, time, 0);
##gid_io.Flush()
##gid_io.CloseResultFile();
##
##print " output of fluid "
##file_name = "fluid"
##
##gid_io.ChangeOutputName(file_name,GiDPostMode.GiD_PostBinary);
##gid_io.WriteConditions(fluid_model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);
##gid_io.WriteNodalResults(DISPLACEMENT, fluid_model_part.Nodes, time, 0);
##gid_io.Flush()
##gid_io.CloseResultFile();
##
          
        

