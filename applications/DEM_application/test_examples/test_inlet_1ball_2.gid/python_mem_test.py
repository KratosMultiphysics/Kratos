from KratosMultiphysics import *
import KratosMultiphysics.DEMApplication
import os

import main_script as Main

for x in range(0, 100000):
    print(x)
    solution = Main.Solution()   
    solution.AddVariables()

    solution.ReadModelParts()

    solution.FillAnalyticSubModelParts()

    # Setting up the buffer size
    solution.procedures.SetUpBufferSizeInAllModelParts(solution.spheres_model_part, 1, solution.cluster_model_part, 1, solution.DEM_inlet_model_part, 1, solution.rigid_face_model_part, 1)

    # Adding dofs
    solution.solver.AddDofs(solution.spheres_model_part)
    solution.solver.AddDofs(solution.cluster_model_part)
    solution.solver.AddDofs(solution.DEM_inlet_model_part)

    os.chdir(solution.main_path)
    solution.KRATOSprint("\nInitializing Problem...")    
    
    solution.GraphicalOutputInitialize()        

# Perform a partition to balance the problem
solution.solver.search_strategy = solution.parallelutils.GetSearchStrategy(solution.solver, solution.spheres_model_part)
solution.solver.BeforeInitialize()
solution.parallelutils.Repart(solution.spheres_model_part)

#Setting up the BoundingBox
solution.bounding_box_time_limits = solution.procedures.SetBoundingBoxLimits(solution.all_model_parts, solution.creator_destructor)

#Finding the max id of the nodes... (it is necessary for anything that will add spheres to the solution.spheres_model_part, for instance, the INLETS and the CLUSTERS read from mdpa file.z
max_Id = solution.procedures.FindMaxNodeIdAccrossModelParts(solution.creator_destructor, solution.all_model_parts)

solution.creator_destructor.SetMaxNodeId(solution.all_model_parts.MaxNodeId)

#Strategy Initialization
os.chdir(solution.main_path)
solution.solver.Initialize() # Possible modifications of number of elements and number of nodes

#Constructing a model part for the DEM inlet. It contains the DEM elements to be released during the simulation
#Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itsolution according to some options of the DEM model part
solution.SetInlet()

solution.SetInitialNodalValues()

solution.DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(solution.DEM_parameters, solution.graphs_path, solution.spheres_model_part, solution.rigid_face_model_part)

os.chdir(solution.graphs_path)
solution.DEMEnergyCalculator = DEM_procedures.DEMEnergyCalculator(solution.DEM_parameters, solution.spheres_model_part, solution.cluster_model_part, "EnergyPlot.grf")

solution.materialTest.Initialize(solution.DEM_parameters, solution.procedures, solution.solver, solution.graphs_path, solution.post_path, solution.spheres_model_part, solution.rigid_face_model_part)

solution.KRATOSprint("Initialization Complete" + "\n")

solution.report.Prepare(timer, solution.DEM_parameters["ControlTime"].GetDouble())

#solution.procedures.ModelData(solution.spheres_model_part, solution.solver) #check link with ModelDataInfo = "OFF"

solution.materialTest.PrintChart()
solution.materialTest.PrepareDataForGraph()

solution.post_utils = DEM_procedures.PostUtils(solution.DEM_parameters, solution.spheres_model_part)
solution.report.total_steps_expected = int(solution.final_time / solution.dt)
solution.KRATOSprint(solution.report.BeginReport(timer))
os.chdir(solution.main_path)


solution.Initialize()
del solution
#solution.Run()

