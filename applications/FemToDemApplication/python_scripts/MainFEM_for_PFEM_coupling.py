from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication.MainFEM_for_coupling as MainFEM_for_coupling
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

# Python script created to modify the existing one due to the coupling of the DEM app in 2D

class FEM_for_PFEM_coupling_Solution(MainFEM_for_coupling.FEM_for_coupling_Solution):

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

        # Adding PFEM Variables TODO put in another place
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POISSON_RATIO)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YOUNG_MODULUS)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_ERROR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)


        #VARIABLES FOR PAPANASTASIOU MODEL
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FLOW_INDEX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELD_SHEAR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT)

        #VARIABLES FOR MU-I RHEOLOGY MODEL
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.STATIC_FRICTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.DYNAMIC_FRICTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ZERO)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DIAMETER)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INFINITE_FRICTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ONE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ALPHA_PARAMETER)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION_OLD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION_RATE)

        # PFEM fluid variables
        # self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NORMVELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELDED)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FREESURFACE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_ACCELERATION)

        
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

