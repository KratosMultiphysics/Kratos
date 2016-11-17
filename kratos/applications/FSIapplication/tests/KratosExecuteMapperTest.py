from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.FSIApplication import *

import os
import process_factory
import NonConformant_OneSideMap

class KratosExecuteMapperTest:

    def __init__(self, ProjectParameters):
        
        # Json format solvers settings
        ProjectParametersFluid = ProjectParameters["fluid_solver_settings"]
        ProjectParametersSolid = ProjectParameters["structure_solver_settings"]
        
        # Defining a model part for the fluid and one for the structure
        self.structure_main_model_part = ModelPart("structure_part")
        self.fluid_main_model_part = ModelPart("fluid_part")

        # Set the domain size (2D test)
        self.structure_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, 2)
        self.fluid_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, 2)

        # Create the fluid and solid solvers
        fluid_solver_module = __import__(ProjectParametersFluid["solver_settings"]["solver_type"].GetString())
        fluid_solver = fluid_solver_module.CreateSolver(self.fluid_main_model_part,ProjectParametersFluid["solver_settings"])
        solid_solver_module = __import__(ProjectParametersSolid["solver_settings"]["solver_type"].GetString())
        solid_solver = solid_solver_module.CreateSolver(self.structure_main_model_part,ProjectParametersSolid["solver_settings"])

        FluidModel = {ProjectParameters["structure_solver_settings"]["problem_data"]["model_part_name"].GetString() : self.structure_main_model_part}
        SolidModel = {ProjectParameters["fluid_solver_settings"]["problem_data"]["model_part_name"].GetString() : self.fluid_main_model_part}

        # Variables addition
        fluid_solver.AddVariables()
        solid_solver.AddVariables()
        NonConformant_OneSideMap.AddVariables(self.fluid_main_model_part, self.structure_main_model_part)

        # Import model_parts
        fluid_solver.ImportModelPart()
        solid_solver.ImportModelPart()

        # Add DOFs
        fluid_solver.AddDofs()
        solid_solver.AddDofs()
                
        # Get the list of the skin submodel parts where the fluid interface submodel part is stored
        for i in range(ProjectParameters["fluid_solver_settings"]["solver_settings"]["skin_parts"].size()):
            skin_part_name = ProjectParameters["fluid_solver_settings"]["solver_settings"]["skin_parts"][i].GetString()
            FluidModel.update({skin_part_name: self.fluid_main_model_part.GetSubModelPart(skin_part_name)})

        # Get the list of the submodel parts where the structure interface is stored
        for i in range(ProjectParameters["structure_solver_settings"]["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = ProjectParameters["structure_solver_settings"]["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            SolidModel.update({part_name: self.structure_main_model_part.GetSubModelPart(part_name)})

        # Fluid domain processes
        self.list_of_processes = process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["fluid_solver_settings"]["boundary_conditions_process_list"] )

        # Solid domain processes
        self.list_of_processes += process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( ProjectParameters["structure_solver_settings"]["constraints_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( ProjectParameters["structure_solver_settings"]["loads_process_list"] )
            
        # Initialize solvers
        fluid_solver.Initialize()
        solid_solver.Initialize()
        
        for process in self.list_of_processes:
            process.ExecuteInitialize()
        
        # Mapper construction
        search_radius_factor = 1.0
        mapper_max_iteration = 25
        self.mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(self.fluid_main_model_part, 
                                                                        self.structure_main_model_part, 
                                                                        search_radius_factor,
                                                                        mapper_max_iteration)

        # Output settings 
        #~ self.problem_path = os.getcwd()
        
        #~ output_post_solid = ProjectParametersSolid.Has("output_configuration")
        #~ output_post_fluid = ProjectParametersFluid.Has("output_configuration")

        #~ if (output_post_solid == True) and (output_post_fluid == True):
            #~ self.output_post = True
        
        #~ if (self.output_post == True):
            #~ from gid_output_process import GiDOutputProcess
            
            #~ self.gid_output_structure = GiDOutputProcess(solid_solver.GetComputingModelPart(),
                                        #~ ProjectParameters["structure_solver_settings"]["problem_data"]["problem_name"].GetString()+"_structure",
                                        #~ ProjectParameters["structure_solver_settings"]["output_configuration"])

            #~ self.gid_output_fluid = GiDOutputProcess(fluid_solver.GetComputingModelPart(),
                                    #~ ProjectParameters["fluid_solver_settings"]["problem_data"]["problem_name"].GetString()+"_fluid",
                                    #~ ProjectParameters["fluid_solver_settings"]["output_configuration"])

            #~ self.gid_output_structure.ExecuteInitialize()
            #~ self.gid_output_fluid.ExecuteInitialize()
        
        
        #~ if (self.output_post == True):
            #~ self.gid_output_structure.ExecuteBeforeSolutionLoop()
            #~ self.gid_output_fluid.ExecuteBeforeSolutionLoop()


    def Solve(self):
        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()
            
        self.structure_main_model_part.ProcessInfo[TIME_STEPS] = 1
        self.structure_main_model_part.CloneTimeStep(1.0)
        self.fluid_main_model_part.ProcessInfo[TIME_STEPS] = 1
        self.fluid_main_model_part.CloneTimeStep(1.0)

        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()
            
        #~ if (self.output_post == True):
            #~ self.gid_output_structure.ExecuteInitializeSolutionStep()
            #~ self.gid_output_fluid.ExecuteInitializeSolutionStep()
                    
        # Solve: Information transfer
        self.mapper.FluidToStructure_ScalarMap(PRESSURE, PRESSURE,   True)
        self.mapper.FluidToStructure_VectorMap(REACTION, POINT_LOAD, True, True)
        self.mapper.StructureToFluid_VectorMap(VELOCITY, VELOCITY,   True, False)
        
        #~ if (self.output_post == True):
            #~ self.gid_output_structure.ExecuteFinalizeSolutionStep()
            #~ self.gid_output_fluid.ExecuteFinalizeSolutionStep()

        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        for process in self.list_of_processes:
            process.ExecuteBeforeOutputStep()

        for process in self.list_of_processes:
            process.ExecuteAfterOutputStep()

        #~ if (self.output_post == True):
            #~ self.gid_output_structure.PrintOutput()
            #~ self.gid_output_fluid.PrintOutput()

        #~ if (self.output_post == True):
            #~ self.gid_output_structure.ExecuteFinalize()
            #~ self.gid_output_fluid.ExecuteFinalize()

        for process in self.list_of_processes:
            process.ExecuteFinalize()
