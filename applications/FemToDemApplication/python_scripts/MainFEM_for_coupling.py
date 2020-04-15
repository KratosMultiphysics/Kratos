from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication.MainFemDem as MainFemDem
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.DemStructuresCouplingApplication as DEM_Structures

# Python script created to modify the existing one due to the coupling of the DEM app in 2D

class FEM_for_coupling_Solution(MainFemDem.FEM_Solution):

    def Info(self):
        print("FEM part of the FEMDEM application")


    def Initialize(self):

        #### INITIALIZE ####

        # Add variables (always before importing the model part)
        self.solver.AddVariables()

        # For remeshing purposes
        self.main_model_part.AddNodalSolutionStepVariable(KratosFemDem.NODAL_STRESS_VECTOR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFemDem.EQUIVALENT_NODAL_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFemDem.EQUIVALENT_NODAL_STRESS_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFemDem.NODAL_DAMAGE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFemDem.EQUIVALENT_STRESS_VM)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFemDem.DISPLACEMENT_INCREMENT)
        self.main_model_part.AddNodalSolutionStepVariable(DEM.DEM_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TOTAL_FORCES)
        self.main_model_part.AddNodalSolutionStepVariable(DEM.DELTA_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(DEM.CONTACT_FORCES)
        self.main_model_part.AddNodalSolutionStepVariable(DEM.ELASTIC_FORCES)
        self.main_model_part.AddNodalSolutionStepVariable(DEM.TANGENTIAL_ELASTIC_FORCES)
        self.main_model_part.AddNodalSolutionStepVariable(DEM.SHEAR_STRESS)

        # For the Substepping
        self.main_model_part.AddNodalSolutionStepVariable(DEM_Structures.BACKUP_LAST_STRUCTURAL_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(DEM_Structures.BACKUP_LAST_STRUCTURAL_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(DEM_Structures.SMOOTHED_STRUCTURAL_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(DEM.CONTACT_IMPULSE)


        # Read model_part (note: the buffer_size is set here) (restart is read here)
        self.solver.ImportModelPart()

        # Add dofs (always after importing the model part)
        if((self.main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
                self.solver.AddDofs()
        else:
            self.solver.AddDofs()

        # Add materials (assign material to model_parts if Materials.json exists)
        self.AddMaterials()

        # Add processes
        self.model_processes = self.AddProcesses()
        self.model_processes.ExecuteInitialize()

        # Print model_part and properties
        if(self.echo_level > 1):
            print("")
            print(self.main_model_part)
            for properties in self.main_model_part.Properties:
                print(properties)

        #### START SOLUTION ####
        self.computing_model_part = self.solver.GetComputingModelPart()

        ## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
        self.solver.Initialize()
        #self.solver.InitializeStrategy()
        self.solver.SetEchoLevel(self.echo_level)

        # Initialize GiD  I/O (gid outputs, file_lists)
        self.SetGraphicalOutput()

        self.GraphicalOutputExecuteInitialize()

        print(" ")
        print("=================================================")
        print(" - Kratos FemDem Application Calculation Start - ")
        print("=================================================")

        self.model_processes.ExecuteBeforeSolutionLoop()

        self.GraphicalOutputExecuteBeforeSolutionLoop()

        # Set time settings
        self.step       = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        self.time       = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.end_time   = self.ProjectParameters["problem_data"]["end_time"].GetDouble()
        self.delta_time = self.ComputeDeltaTime()


#============================================================================================================================

    def ComputeDeltaTime(self):

        if self.ProjectParameters["problem_data"].Has("time_step"):
            return self.ProjectParameters["problem_data"]["time_step"].GetDouble()

        elif self.ProjectParameters["problem_data"].Has("variable_time_steps"):

            current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
            for key in self.ProjectParameters["problem_data"]["variable_time_steps"].keys():
                interval_settings = self.ProjectParameters["problem_data"]["variable_time_steps"][key]
                interval = KratosMultiphysics.IntervalUtility(interval_settings)
                # Getting the time step of the interval
                if interval.IsInInterval(current_time):
                    return interval_settings["time_step"].GetDouble()
                # If we arrive here we raise an error because the intervals are not well defined
                raise Exception("::[MechanicalSolver]:: Time stepping not well defined!")
        else:
            raise Exception("::[MechanicalSolver]:: Time stepping not defined!")