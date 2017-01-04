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
        self.structure_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParametersFluid["problem_data"]["domain_size"].GetInt())
        self.fluid_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParametersSolid["problem_data"]["domain_size"].GetInt())

        # Set the fluid and solid models
        FluidModel = {ProjectParametersFluid["problem_data"]["model_part_name"].GetString() : self.structure_main_model_part}
        SolidModel = {ProjectParametersSolid["problem_data"]["model_part_name"].GetString() : self.fluid_main_model_part}

        # Fluid model part variables addition
        self.fluid_main_model_part.AddNodalSolutionStepVariable(VELOCITY)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(PRESSURE)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(REACTION)

        # Structure model part variables addition
        self.structure_main_model_part.AddNodalSolutionStepVariable(VELOCITY)
        self.structure_main_model_part.AddNodalSolutionStepVariable(PRESSURE)
        self.structure_main_model_part.AddNodalSolutionStepVariable(POINT_LOAD)

        # Mapper variables addition
        NonConformant_OneSideMap.AddVariables(self.fluid_main_model_part, self.structure_main_model_part)

        # Fluid domain model reading
        ModelPartIO(ProjectParametersFluid["solver_settings"]["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.fluid_main_model_part)
        prepare_model_part_settings_fluid = Parameters("{}")
        prepare_model_part_settings_fluid.AddValue("volume_model_part_name",ProjectParametersFluid["solver_settings"]["volume_model_part_name"])
        prepare_model_part_settings_fluid.AddValue("skin_parts",ProjectParametersFluid["solver_settings"]["skin_parts"])
        import check_and_prepare_model_process_fluid
        check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.fluid_main_model_part, prepare_model_part_settings_fluid).Execute()

        # Solid domain model reading
        computing_model_part_name = "computing_domain"
        ModelPartIO(ProjectParametersSolid["solver_settings"]["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.structure_main_model_part)
        prepare_model_part_settings_structure = Parameters("{}")
        prepare_model_part_settings_structure.AddEmptyValue("computing_model_part_name").SetString(computing_model_part_name)
        prepare_model_part_settings_structure.AddValue("problem_domain_sub_model_part_list",ProjectParametersSolid["solver_settings"]["problem_domain_sub_model_part_list"])
        prepare_model_part_settings_structure.AddValue("processes_sub_model_part_list",ProjectParametersSolid["solver_settings"]["processes_sub_model_part_list"])
        import check_and_prepare_model_process_solid
        check_and_prepare_model_process_solid.CheckAndPrepareModelProcess(self.structure_main_model_part, prepare_model_part_settings_structure).Execute()

        # Get the list of the skin submodel parts where the fluid interface submodel part is stored
        for i in range(ProjectParametersFluid["solver_settings"]["skin_parts"].size()):
            skin_part_name = ProjectParametersFluid["solver_settings"]["skin_parts"][i].GetString()
            FluidModel.update({skin_part_name: self.fluid_main_model_part.GetSubModelPart(skin_part_name)})

        # Get the list of the submodel parts where the structure interface is stored
        for i in range(ProjectParametersSolid["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = ProjectParametersSolid["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            SolidModel.update({part_name: self.structure_main_model_part.GetSubModelPart(part_name)})

        # Construct processes
        self.list_of_processes = process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParametersFluid["boundary_conditions_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( ProjectParametersSolid["constraints_process_list"] )

        # Set fluid and structure interfaces
        for node in self.fluid_main_model_part.GetSubModelPart("Fluid_interface").Nodes:
            node.Set(INTERFACE, True)
        for node in self.structure_main_model_part.GetSubModelPart("Structure_interface").Nodes:
            node.Set(INTERFACE, True)

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
        self.output_post = False # Set this variable to True if it is need to print the results for debugging purposes
        self.problem_path = os.getcwd()

        if (self.output_post == True):
            from gid_output_process import GiDOutputProcess

            self.gid_output_structure = GiDOutputProcess(self.structure_main_model_part,
                                                         ProjectParametersSolid["problem_data"]["problem_name"].GetString()+"_structure",
                                                         ProjectParametersSolid["output_configuration"])

            self.gid_output_fluid = GiDOutputProcess(self.fluid_main_model_part,
                                                     ProjectParametersFluid["problem_data"]["problem_name"].GetString()+"_fluid",
                                                     ProjectParametersFluid["output_configuration"])

            self.gid_output_structure.ExecuteInitialize()
            self.gid_output_fluid.ExecuteInitialize()

    def Solve(self):

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        if (self.output_post == True):
            self.gid_output_structure.ExecuteBeforeSolutionLoop()
            self.gid_output_fluid.ExecuteBeforeSolutionLoop()

        self.structure_main_model_part.ProcessInfo[TIME_STEPS] = 1
        self.structure_main_model_part.CloneTimeStep(1.0)
        self.fluid_main_model_part.ProcessInfo[TIME_STEPS] = 1
        self.fluid_main_model_part.CloneTimeStep(1.0)

        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()

        if (self.output_post == True):
            self.gid_output_structure.ExecuteInitializeSolutionStep()
            self.gid_output_fluid.ExecuteInitializeSolutionStep()

        # Solve: Information transfer
        self.mapper.FluidToStructure_ScalarMap(PRESSURE, PRESSURE,   True)
        self.mapper.FluidToStructure_VectorMap(REACTION, POINT_LOAD, True, True)
        self.mapper.StructureToFluid_VectorMap(VELOCITY, VELOCITY,   True, False)

        if (self.output_post == True):
            self.gid_output_structure.ExecuteFinalizeSolutionStep()
            self.gid_output_fluid.ExecuteFinalizeSolutionStep()

        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        for process in self.list_of_processes:
            process.ExecuteBeforeOutputStep()

        for process in self.list_of_processes:
            process.ExecuteAfterOutputStep()

        if (self.output_post == True):
            self.gid_output_structure.PrintOutput()
            self.gid_output_fluid.PrintOutput()

        if (self.output_post == True):
            self.gid_output_structure.ExecuteFinalize()
            self.gid_output_fluid.ExecuteFinalize()

        for process in self.list_of_processes:
            process.ExecuteFinalize()
