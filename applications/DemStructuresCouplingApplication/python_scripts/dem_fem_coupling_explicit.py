##  Main authors: Klaus B. Sautter

import KratosMultiphysics
from KratosMultiphysics import StructuralMechanicsApplication
import structural_mechanics_analysis

from KratosMultiphysics import DEMApplication
import main_script

import KratosMultiphysics.MappingApplication as KratosMapping


def RunCoupledSystem():
    ## create structural analysis
    model = KratosMultiphysics.Model()
    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    structural_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model,parameters)

    structural_analysis.Initialize()

    print("-----> Initialized fem part")

    ## create dem analysis
    model_dem = KratosMultiphysics.Model()
    dem_analysis = main_script.Solution(model_dem)
    dem_analysis.Initialize()
    dem_analysis.InitializeTime()

    print("-----> Initialized dem part")

    ## create model parts
    mp_struct = model["Structure.computing_domain"]

    mp_dem = dem_analysis.rigid_face_model_part.GetSubModelPart('1')
    mp_dem_particle = dem_analysis.spheres_model_part

    ####################################################################### WALL
    ## create point load condition on all fem nodes
    cond_counter = 0
    node_id_list = []
    for node_i in mp_struct.Nodes:
        node_id = node_i.Id
        node_id_list.append(node_id)
        cond_counter+=1
        mp_struct.CreateNewCondition("PointLoadCondition3D1N",cond_counter,[node_id],mp_struct.GetProperties()[0])

    print("-----> Created fem PointLoadConditions")

    ## create structural submodal part used for mapping
    struct_smp = mp_struct.CreateSubModelPart("struct_sub");
    struct_smp.AddNodes(node_id_list)
    ## set up mapper
    print("-----> Added node list for mapping part")

    mapper_settings = KratosMultiphysics.Parameters("""{"mapper_settings":
                    {"mapper_type": "nearest_neighbor",
                    "interface_submodel_part_destination": "struct_sub",
                    "echo_level":1}}""")

    print("-----> Wrote mapper settings")
    mapper = KratosMapping.MapperFactory.CreateMapper(
        mp_dem,
        mp_struct,
        mapper_settings["mapper_settings"])


    print("-----> Initialized mapper")


    dem_mesh_moving_utility = DEMApplication.MoveMeshUtility()

    print("-----> Starting time loop:")
    ## run solving loop
    while dem_analysis.time < dem_analysis.final_time:
        print('current t: ',dem_analysis.time)

        ## call dem functions
        dem_analysis.UpdateTimeParameters()

        dem_analysis.BeforeSolveOperations(dem_analysis.time)
        dem_analysis.SolverSolve()
        dem_analysis.AfterSolveOperations()

        dem_analysis.FinalizeSingleTimeStep()

        ## map loads
        mapper.Map(DEMApplication.CONTACT_FORCES, StructuralMechanicsApplication.POINT_LOAD)

        ## call fem functions
        structural_analysis.time = structural_analysis._GetSolver().AdvanceInTime(structural_analysis.time)
        structural_analysis.InitializeSolutionStep()
        structural_analysis._GetSolver().Predict()
        structural_analysis._GetSolver().SolveSolutionStep()
        structural_analysis.FinalizeSolutionStep()
        structural_analysis.OutputSolutionStep()


        ## map-1 vel,disp
        mapper.InverseMap(KratosMultiphysics.VELOCITY, KratosMultiphysics.VELOCITY)
        mapper.InverseMap(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.DISPLACEMENT)

        ## move dem wall mesh
        dem_mesh_moving_utility.MoveDemMesh(mp_dem.Nodes,True)

        dem_analysis.OutputSingleTimeLoop()   ## do this at the end of the loop to see deformed net


    ## finalize dem
    dem_analysis.Finalize()
    dem_analysis.CleanUpOperations()

    ## finalize fem
    structural_analysis.Finalize()



if __name__ == "__main__":
    RunCoupledSystem()








