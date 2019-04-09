'''
Main Script for FSI with Kratos Mutliphysics

This script is intended to be modified. Each solver can be imported and used as "BlackBox"

Chair of Structural Analysis, Technical University of Munich
All rights reserved
'''
'''
This example is based on the dissertation of Daniel Mok
"Partitionierte Lösungsansätze in der Strukturdynamik und der Fluid-Struktur-Interaktion"
Chapter 7.3 "Flexible Klappe in Kanalströmung mit Einschnürung"
'''
# ----- Importing the modules -----
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.FluidDynamicsApplication as KratosFluidDynamics
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructuralMechanics


# Import the "BlackBox" Solvers
from structural_mechanics_analysis import StructuralMechanicsAnalysis
from fluid_dynamics_analysis import FluidDynamicsAnalysis

import fsi_utilities # here auxiliary functions e.g. for relaxation are declared

fluid_model = KratosMultiphysics.Model()
structural_model = KratosMultiphysics.Model()

fluid_project_params_file_name = "ProjectParametersCFD.json"
with open(fluid_project_params_file_name,'r') as parameter_file:
    parameters_fluid = KratosMultiphysics.Parameters(parameter_file.read())
structural_project_params_file_name = "ProjectParametersCSM.json"
with open(structural_project_params_file_name,'r') as parameter_file:
    parameters_structure = KratosMultiphysics.Parameters(parameter_file.read())

'''
# --------------------------------------------------------
# ----- Setting up and initializing the Fluid Solver -----
# --------------------------------------------------------
'''
fluid_solver = FluidDynamicsAnalysis(fluid_model, parameters_fluid)

fluid_solver.Initialize()

fluid_model_part = fluid_model["MainModelPart"]

print("======================================================================")
print("||||||||||||||||||||||| SETTING UP FLUID DONE ||||||||||||||||||||||||")
print("======================================================================")

'''
# -------------------------------------------------------------
# ----- Setting up and initializing the Structural Solver -----
# -------------------------------------------------------------
'''
structural_solver = StructuralMechanicsAnalysis(structural_model, parameters_structure)

structural_solver.Initialize()

structural_model_part = structural_model["Structure"]

print("======================================================================")
print("||||||||||||||||| SETTING UP STRUCTURAL DYNAMICS DONE ||||||||||||||||")
print("======================================================================")

'''
# ------------------------------------------------------
# ----- Setting up the FSI-related functionalities -----
# ------------------------------------------------------
# '''
# ----- Setting up the time parameters -----
start_time = 0.0
end_time   = 15.0
delta_time = 0.001

num_steps = int((end_time - start_time) / delta_time)
round_val = fsi_utilities.TimeRoundValue(delta_time)

time = start_time
step = 0

# ----- Setting up the FSI Parameters -----
# FSI parameters
max_iter = 10    # number of inner iterations (set to 1 for explicit coupling)
interface_epsilon = 1e-5  # interface residual (only needed for implicit coupling)
relaxation_coefficient = 0.125 # initial value

# ---------------
# ----- ALE -----
# ---------------
fluid_solver._GetSolver().GetMeshMotionSolver().SetEchoLevel(0) # Fix until the source of the prints is found
# -------------------
# ----- Mapping -----
# -------------------
parameter_file = open(fluid_project_params_file_name,'r') # Here the time information of the fluid solver is used (can be changed if desired)
mapper_params = KratosMultiphysics.Parameters( parameter_file.read())

project_parameters_mapper_1 = mapper_params["mapper_settings"][0]
project_parameters_mapper_2 = mapper_params["mapper_settings"][1]
project_parameters_mapper_3 = mapper_params["mapper_settings"][2]

mapper_1 = KratosMapping.MapperFactory.CreateMapper(structural_model_part,fluid_model_part, project_parameters_mapper_1)
mapper_2 = KratosMapping.MapperFactory.CreateMapper(structural_model_part,fluid_model_part, project_parameters_mapper_2)
mapper_3 = KratosMapping.MapperFactory.CreateMapper(structural_model_part,fluid_model_part, project_parameters_mapper_3)

def NeumannToStructure(mapper, flag):
    mapper.InverseMap(KratosStructuralMechanics.POINT_LOAD, KratosMultiphysics.REACTION, flag)

def DisplacementToMesh(mapper):
    mapper.Map(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.MESH_DISPLACEMENT)

print("======================================================================")
print("|||||||||||||||||||||||| SETTING UP FSI DONE |||||||||||||||||||||||||")
print("======================================================================")

file_writer = fsi_utilities.FileWriter("Mok_Results.dat", ["Time", "Disp_X", "Disp_Y", "Disp_Z", "Coupling_Iterations"])
tip_node = structural_model_part.GetNode(1)

# ----- Solving the problem (time integration) -----
while(time <= end_time):
    new_time_fluid = fluid_solver._GetSolver().AdvanceInTime(time)
    new_time_structure = structural_solver._GetSolver().AdvanceInTime(time)

    fluid_solver._GetSolver().Predict()
    structural_solver._GetSolver().Predict()

    fluid_solver.InitializeSolutionStep()
    structural_solver.InitializeSolutionStep()

    time = time + delta_time
    if abs(time-new_time_fluid) > 1e-12:
        raise Exception("Fluid has wrong time!")
    if abs(time-new_time_structure) > 1e-12:
        raise Exception("Structure has wrong time!")
    step += 1

    print("\n--- Step =", step, "/", num_steps, "---")
    print("--- Time =", round(time, round_val), "/", end_time, "---")

    residual = 1
    old_displacements = fsi_utilities.GetDisplacements(structural_model_part.GetSubModelPart("GENERIC_Beam").Nodes, 2)

    num_inner_iter = 1
    ### Inner FSI Loop (executed once in case of explicit coupling)
    for k in range(max_iter):

        # Apply Dirichlet B.C.'s from structural solver to mesh solver
        DisplacementToMesh(mapper_1)
        DisplacementToMesh(mapper_2)
        DisplacementToMesh(mapper_3)

        # Mesh and Fluid are currently solved independently, since the ALE solver does not copy the mesh velocity
        # Solve Mesh
        fluid_solver._GetSolver().SolveSolutionStep()

        # Apply Neumann B.C.'s from fluid solver to structural solver
        NeumannToStructure(mapper_1, KratosMapping.Mapper.SWAP_SIGN | KratosMapping.Mapper.CONSERVATIVE)
        NeumannToStructure(mapper_2, KratosMapping.Mapper.SWAP_SIGN | KratosMapping.Mapper.ADD_VALUES | KratosMapping.Mapper.CONSERVATIVE)

        # # Solver Structure
        structural_solver._GetSolver().SolveSolutionStep()

        # Convergence Checking (only for implicit coupling)
        if max_iter > 1:
            displacements = fsi_utilities.GetDisplacements(structural_model_part.GetSubModelPart("GENERIC_Beam").Nodes, 2)

            # Compute Residual
            old_residual = residual
            residual = fsi_utilities.CalculateResidual(displacements,old_displacements)

            if (fsi_utilities.Norm(residual) <= interface_epsilon):
                fsi_utilities.SetDisplacements(displacements,structural_model_part.GetSubModelPart("GENERIC_Beam").Nodes, 2)
                print("******************************************************")
                print("************ CONVERGENCE AT INTERFACE ACHIEVED *******")
                print("******************************************************")
                break # TODO check if this works bcs it is nested
            else:
                relaxation_coefficient = fsi_utilities.ComputeAitkenRelaxation(relaxation_coefficient, residual, old_residual, k)
                relaxed_displacements = fsi_utilities.CalculateRelaxation(relaxation_coefficient, old_displacements, residual)
                old_displacements = relaxed_displacements
                fsi_utilities.SetDisplacements(relaxed_displacements, structural_model_part.GetSubModelPart("GENERIC_Beam").Nodes, 2)
                num_inner_iter += 1

            if (k+1 >= max_iter):
                print("######################################################")
                print("##### CONVERGENCE AT INTERFACE WAS NOT ACHIEVED ######")
                print("######################################################")

            print("==========================================================")
            print("COUPLING RESIDUAL = ", fsi_utilities.Norm(residual))
            print("COUPLING ITERATION = ", k+1, "/", max_iter)
            print("RELAXATION COEFFICIENT = ",relaxation_coefficient)
            print("==========================================================")

    fluid_solver.FinalizeSolutionStep()
    structural_solver.FinalizeSolutionStep()

    fluid_solver.OutputSolutionStep()
    structural_solver.OutputSolutionStep()

    disp = tip_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
    file_writer.WriteToFile([time, disp[0], disp[1], disp[2], num_inner_iter])

# TIME LOOP END
fluid_solver.Finalize()
structural_solver.Finalize()

file_writer.CloseFile()