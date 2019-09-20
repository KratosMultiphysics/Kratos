from __future__ import print_function, absolute_import, division #makes kratos backward compatible with python 2.6 and 2.7 - breaks compatibility with 2.5
import ProjectParameters
domain_size = 2
#
#
import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import * #need this for the cube mesher
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *


# defining a model part for the fluid
base_model_part = ModelPart("BasePart")
base_model_part.AddNodalSolutionStepVariable(VELOCITY) ##needed to compute the normal
base_model_part.AddNodalSolutionStepVariable(TEMPERATURE) ##needed to compute the normal


#reading the base model part
model_part_io_fluid = ModelPartIO("base")
model_part_io_fluid.ReadModelPart(base_model_part)



# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
base_model_part.SetBufferSize(2)




#refine the mesh once
#for elem in base_model_part.Elements:
    #elem.SetValue(SPLIT_ELEMENT,True)
 
####compute the nodal neighbours on the initial mesh
#number_of_avg_elems = 20
#number_of_avg_nodes = 20
#nodal_neighbour_search = FindNodalNeighboursProcess(base_model_part,number_of_avg_elems,number_of_avg_nodes)
#nodal_neighbour_search.Execute()

####perform the refinement
#Refine = LocalRefineTriangleMesh(base_model_part)
#refine_on_reference = False;
#interpolate_internal_variables = False;
#Refine.LocalRefineMesh(refine_on_reference,interpolate_internal_variables) 
 
#assign velocity field
import math
xc = 0.7 - 0.5
yc = 0.7 - 0.5
for node in base_model_part.Nodes:
    xlocal = node.X - 0.5
    ylocal = node.Y - 0.5
    r = math.sqrt( (xlocal)**2 + (ylocal)**2 )
    if( r < 0.45):
        node.SetSolutionStepValue(VELOCITY_X,0,-ylocal)
        node.SetSolutionStepValue(VELOCITY_Y,0,xlocal)
    
        d = math.sqrt( (xlocal - xc)**2 + (ylocal-yc)**2 ) - 0.1
        
        node.SetSolutionStepValue(TEMPERATURE,0,d)
        #if(d <= 0.0):
            #node.SetSolutionStepValue(TEMPERATURE,0,-1.0)

print("base_model_part =",base_model_part)

# initialize GiD  I/O
from gid_output import GiDOutput
gid_io_base = GiDOutput("base",
                   ProjectParameters.VolumeOutput,
                   ProjectParameters.GiDPostMode,
                   ProjectParameters.GiDMultiFileFlag,
                   ProjectParameters.GiDWriteMeshFlag,
                   ProjectParameters.GiDWriteConditionsFlag)


#printing the mesh of the base
gid_io_base.initialize_results(base_model_part)


#mount the search structure
locator = BinBasedFastPointLocator2D(base_model_part)
locator.UpdateSearchDatabase()
#locator.UpdateSearchDatabaseAssignedSize(0.01)

#construct the utility to move the points
bfecc_utility = BFECCConvection2D(locator)

base_model_part.CloneTimeStep(0.00)

substepping  = 10.0
final_time = 10.0# math.pi
t=0.0
dt = 0.05 #1
step = 0
while( t < final_time):
    t = t+dt
    base_model_part.CloneTimeStep(t)
    
    print("time = ",t)
    bfecc_utility.BFECCconvect(base_model_part,TEMPERATURE,VELOCITY,substepping)
        
    gid_io_base.write_results(
            t,
            base_model_part,
            ProjectParameters.nodal_results,
            ProjectParameters.gauss_points_results)


    step +=1

gid_io_base.finalize_results()
