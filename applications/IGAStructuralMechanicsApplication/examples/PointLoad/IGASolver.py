from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#### TIME MONITORING START ####

# time control starts
import time as timer
print(timer.ctime())
# measure process time
t0p = timer.clock()
# measure wall time
t0w = timer.time()

#
def StartTimeMeasuring():
    # measure process time
    time_ip = timer.clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    # measure process time
    time_fp = timer.clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[IGA Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

#### TIME MONITORING END ####

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.IGAStructuralMechanicsApplication import *
from KratosMultiphysics.NurbsBrepApplication import *
import KratosMultiphysics.MappingApplication as KratosMapping
from KratosMultiphysics.ExternalSolversApplication import *


from KratosExternalSolversApplication import *

import iga_structural_mechanics_solver

#region Set Up Model Part / Parsing the Parameters
#### PARSING THE PARAMETERS ####

#import define_output
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

#set echo level
echo_level = ProjectParameters["solver_settings"]["echo_level"].GetInt()

iga_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
# Not yet implemented in the parser
iga_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())


###TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : iga_model_part}


#construct the iga_solver (main setting methods are located in the solver_module)
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
iga_solver = solver_module.CreateSolver(iga_model_part, ProjectParameters["solver_settings"])

iga_solver.AddVariables()

##NURBSBREPAPPLICATION
nurbs_brep_time = StartTimeMeasuring()
import nurbs_brep_process
NurbsBrepProcess = nurbs_brep_process.Factory(ProjectParameters["nurbs_brep_configuration"], iga_model_part)

NurbsBrepProcess.ExecuteInitialize()
NurbsBrepProcess.ExecuteFinalize()
StopTimeMeasuring(nurbs_brep_time, "NurbsBrepApplication time", True)

import read_materials_process
read_materials_process.Factory(ProjectParameters,Model)

for properties in iga_model_part.Properties:
    print(properties)

nurbs_brep_time = StartTimeMeasuring()
iga_solver.ImportModelPartNurbsBrep(NurbsBrepProcess.model_part_integration_domain, ProjectParameters)
StopTimeMeasuring(nurbs_brep_time, "Import model part time", True)
iga_solver.AddDofs()


#build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
##TODO: replace MODEL for the Kratos one ASAP
##get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name: iga_model_part.GetSubModelPart(part_name)})

#print iga_model_part
if(echo_level>1):
    print("")
    print(NurbsBrepProcess.model_part_integration_domain)
    print(iga_model_part)
    #for properties in iga_model_part.Properties:
        #print(properties)
#### iga_model_part settings end ####
#endregion
#region Process Settings
#--- processes settings start ####
import process_factory
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )

##print list of constructed processes
if(echo_level>4):
    for process in list_of_processes:
        print(process)

##TODO: decide which is the correct place to initialize the processes
for process in list_of_processes:
    process.ExecuteInitialize()
#### processes settings end ####
#endregion

#region Start Solution / Initialize / GiD I/O
#### START SOLUTION ####

#TODO: think if there is a better way to do this
iga_computing_model_part = iga_solver.GetComputeModelPart()

print(iga_model_part)


#### output settings start ####

problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()


#### output settings end ####
# initialize IGA_IO
import iga_io_process
IGA_IO = iga_io_process.Factory(ProjectParameters["output_configuration"],Model)

IGA_IO.ExecuteInitialize()

#iga_solver.solver.SetEchoLevel(5)

# initialize iga_solver and processes
iga_solver.Initialize()

for process in list_of_processes:
	process.ExecuteBeforeSolutionLoop()

#endregion

#region Time Integration
## Stepping and time settings (get from process info or solving info)
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
step = 0
time = ProjectParameters["problem_data"]["start_time"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

mapper_params = Parameters("""
{
    "mapper_type" : "iga_dem",
    "interface_submodel_part_origin": "STRUCTURAL_ANALYSIS_2",
    "search_radius": 1.0,
    "search_iterations": 1
}
""")

#mapper = KratosMapping.MapperFactory.CreateMapper(iga_model_part, spheres_mp, mapper_params)

#condition_model_part = iga_model_part.CreateSubModelPart("ConditionModelPart")
iga_solver.solver.SetEchoLevel(5)
# solving the problem (time integration)
while(time <= end_time):

    #TODO: this must be done by a solving_info utility in the iga_solver
    # store previous time step
    #~ iga_computing_model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = delta_time
    # set new time step ( it can change when solve is called )
    #~ delta_time = iga_computing_model_part.ProcessInfo[DELTA_TIME]

    time = time + delta_time
    step = step + 1
    iga_model_part.ProcessInfo[TIME_STEPS] = step
    iga_model_part.CloneTimeStep(time)

    #NurbsBrepProcess.modeler.GetInterfaceConditions(dem_analysis.spheres_model_part, condition_model_part, NurbsBrepProcess.model_part_integration_domain)
    #print(condition_model_part)
    #mapper.UpdateInterface()

	
    #dem_analysis.RunSingleTemporalLoop()
	
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()
    nurbs_brep_time = StartTimeMeasuring()
    iga_solver.Solve()
    #print("ckeck")
    StopTimeMeasuring(nurbs_brep_time, "Solving step time", True)

	
	
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    IGA_IO.ExecuteFinalizeSolutionStep()


for process in list_of_processes:
	process.ExecuteFinalize()

# ending the problem (time integration finished)
# check solving information for any problem
#~ iga_solver.InfoCheck() # InfoCheck not implemented yet.
#endregion

#### END SOLUTION ####
# measure process time
tfp = timer.clock()
# measure wall time
tfw = timer.time()
print("::[IGA Simulation]:: [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")
print(timer.ctime())

# Write results in proper format for Rhino

