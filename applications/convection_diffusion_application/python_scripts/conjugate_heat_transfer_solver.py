from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")
KratosMultiphysics.CheckRegisteredApplications("ConvectionDiffusionApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
# TODO: USE THE MAPPING APPLICATION INSTEAD OF THE FSI MAPPER
import KratosMultiphysics.FSIApplication as FSIApplication
import NonConformant_OneSideMap as ncosm

# Importing the base class
from python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    return ConjugateHeatTransferSolver(main_model_part, custom_settings)

class ConjugateHeatTransferSolver(PythonSolver):
    
    def __init__(self, model, custom_settings):

        super(ConjugateHeatTransferSolver, self).__init__(model, custom_settings)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "conjugate_heat_transfer",
            "domain_size": -1,
            "echo_level": 0,
            "fluid_domain_solver_settings": {
                "fluid_solver_settings": {
                    "solver_type": "navier_stokes_solver_vmsmonolithic",
                    "model_import_settings": {
                        "input_type": "mdpa",
                        "input_filename": "unknown_name"
                    }
                },
                "thermal_solver_settings":{
                    "model_part_name": "FluidThermalModelPart",
                    "solver_type": "Transient",
                    "analysis_type": "linear",
                    "model_import_settings": {
                        "input_type": "use_input_model_part"
                    },
                    "material_import_settings": {
                        "materials_filename": "ThermicMaterialsFluid.json"
                    }
                }
            },
            "solid_domain_solver_settings":{
                "solid_solver_settings": {
                },
                "thermal_solver_settings": {
                    "model_part_name": "SolidThermalModelPart",
                    "solver_type": "Transient",
                    "analysis_type": "linear",
                    "model_import_settings": {
                        "input_type": "mdpa",
                        "input_filename": "unknown_name"
                    },
                    "material_import_settings": {
                        "materials_filename": "ThermicMaterialsSolid.json"
                    }
                }
            },
            "coupling_settings":{
                "max_iteration": 10,
                "relaxation_factor": 0.7,
                "temperature_relative_tolerance": 1e-5,
                "fluid_interfaces_list": [],
                "solid_interfaces_list": []
            }
        }
        """)
        
        ## Overwrite the default settings with user-provided parameters
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        ## Set the fluid dynamics solver
        import python_solvers_wrapper_fluid
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_domain_solver_settings"]["fluid_solver_settings"], "OpenMP")

        # Set the fluid and solid heat solvers
        import python_solvers_wrapper_convection_diffusion
        self.fluid_thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model, self.settings["fluid_domain_solver_settings"]["thermal_solver_settings"], "OpenMP")
        self.solid_thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model, self.settings["solid_domain_solver_settings"]["thermal_solver_settings"], "OpenMP")
        
    def AddVariables(self):
        self.fluid_solver.AddVariables()
        self.fluid_thermal_solver.AddVariables() 
        
        #TODO: WHY ARE WE USING NODAL_PAUX?
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.AUX_FLUX)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_PAUX)
        self.fluid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_PAUX)
        
        # Temporary container for un-relaxed temperature
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part,
                                                             self.fluid_thermal_solver.main_model_part)
        
        self.solid_thermal_solver.AddVariables()
        self.solid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.AUX_FLUX)
        self.solid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_PAUX)
        self.solid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.AUX_TEMPERATURE)

    def ImportModelPart(self):
        # Check that both thermal solvers have a different model part name. If
        # both model part names coincide the solver will fail to acces them. This
        # is the case if the default one in the convection diffusion is taken.
        fluid_thermal_model_part_name = self.settings["fluid_domain_solver_settings"]["thermal_solver_settings"]["model_part_name"].GetString()
        solid_thermal_model_part_name = self.settings["solid_domain_solver_settings"]["thermal_solver_settings"]["model_part_name"].GetString()
        if fluid_thermal_model_part_name == solid_thermal_model_part_name:
            err_msg = "\nFluid thermal solver settings model_part_name and solid thermal solver settings model_part_name can not coincide.\n"
            err_msg += "- fluid model_part_name: " + fluid_thermal_model_part_name + "\n"
            err_msg += "- solid model_part_name: " + solid_thermal_model_part_name + "\n"
            err_msg += "Provide different model_part_names in the JSON settings file."
            raise Exception(err_msg)

        # Import the fluid domain in the fluid dynamics solver
        self.fluid_solver.ImportModelPart()
        
        # In order to consider the buoyancy effects, the nodes in the fluid model part must
        # be shared with the nodes in the fluid thermal model part. To do that, we use the modeler
        # Save the convection diffusion settings
        convection_diffusion_settings = self.fluid_thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        
        # Here the fluid model part is cloned to be thermal model part so that the nodes are shared
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        if(self.domain_size == 2):
            modeler.GenerateModelPart(self.fluid_solver.main_model_part, 
                                      self.fluid_thermal_solver.main_model_part, 
                                      "Element2D3N",
                                      "Condition2D2N")
        else:
            modeler.GenerateModelPart(self.fluid_solver.main_model_part,
                                      self.fluid_thermal_solver.main_model_part,
                                      "Element3D4N",
                                      "Condition3D3N")

        # Set the saved convection diffusion settings to the new thermal model part
        self.fluid_thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)

        # Import the solid domain 
        self.solid_thermal_solver.ImportModelPart()

    def AddDofs(self):
        (self.fluid_solver).AddDofs()
        (self.fluid_thermal_solver).AddDofs()
        (self.solid_thermal_solver).AddDofs()

    def GetComputingModelPart(self):        
        return self.fluid_solver.GetComputingModelPart()

    def ComputeDeltaTime(self):                
        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.fluid_thermal_solver.GetMinimumBufferSize()

        return max(buffer_size_fluid, buffer_size_thermal)

    def Initialize(self):
        # Initialize the fluid and solid solvers
        (self.fluid_solver).Initialize()
        (self.fluid_thermal_solver).Initialize()
        (self.solid_thermal_solver).Initialize()

        # Create the fluid and solid interface mapper
        self._SetUpMapper()

    def Clear(self):
        (self.fluid_solver).Clear()
        (self.fluid_thermal_solver).Clear()
        (self.solid_thermal_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.fluid_thermal_solver).Check()
        (self.solid_thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.fluid_thermal_solver).SetEchoLevel(level)
        (self.solid_thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared between the fluid and thermal solvers
        new_time = self.fluid_solver.AdvanceInTime(current_time)
       
        # Do the time advance in the solid thermal solver
        self.solid_thermal_solver.main_model_part.CloneTimeStep(new_time)
        self.solid_thermal_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def PrepareModelPart(self):  
        self.fluid_solver.PrepareModelPart()
        # TODO: CHECK THIS (if we switch the order in thermal solvers the solver breaks)
        self.solid_thermal_solver.PrepareModelPart()
        self.fluid_thermal_solver.PrepareModelPart()

        self._SetUpDirichletCouplingBoundary(self.solid_thermal_solver.GetComputingModelPart())

    def InitializeSolutionStep(self):
        self.fluid_solver.InitializeSolutionStep()
        self.fluid_thermal_solver.InitializeSolutionStep()
        self.solid_thermal_solver.InitializeSolutionStep()

    def Predict(self):
        self.fluid_solver.Predict()
        self.fluid_thermal_solver.Predict()
        self.solid_thermal_solver.Predict()
        print("finished predict")

    def SolveSolutionStep(self):
        # Solve the buoyancy solver
        self.fluid_solver.SolveSolutionStep()

        max_iteration = self.settings["coupling_settings"]["max_iteration"].GetInt()
        relaxation_factor = self.settings["coupling_settings"]["relaxation_factor"].GetDouble()
        temp_rel_tol = self.settings["coupling_settings"]["temperature_relative_tolerance"].GetDouble()

        # Couple the solid and fluid thermal problems
        iteration = 0
        while iteration < max_iteration:
            # Solve Dirichlet side to get reactions from solid domain
            self.solid_thermal_solver.Solve()

            # Map reactions to the fluid interface. Note that we first call the redistribution utility to convert the point values to distributed ones
            redistribution_tolerance = 1e-5
            redistribution_max_iteration = 50
            KratosMultiphysics.VariableRedistributionUtility.DistributePointValues(self.mapper.str_interface,
                                                                                   KratosMultiphysics.REACTION_FLUX,
                                                                                   KratosConvDiff.AUX_FLUX,
                                                                                   redistribution_tolerance,
                                                                                   redistribution_max_iteration)
            
            self.mapper.StructureToFluid_ScalarMap(KratosConvDiff.AUX_FLUX,
                                                   KratosMultiphysics.FACE_HEAT_FLUX,
                                                   False) # Sign is not kept (inverted) when mapping

            # Solve Neumann side to get temperature values from fluid domain
            self.fluid_thermal_solver.SolveSolutionStep()

            # Map back the Neumann domain obtained temperature
            self.mapper.FluidToStructure_ScalarMap(KratosMultiphysics.TEMPERATURE,
                                                   KratosConvDiff.AUX_TEMPERATURE,
                                                   True) # Sign is kept when mapping

            temperature_difference = 0.0
            for node in self.mapper.str_interface.Nodes:
                old_temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                new_temperature = node.GetSolutionStepValue(KratosConvDiff.AUX_TEMPERATURE)

                # Accumulate the squared nodal residal
                temperature_difference += (old_temperature - new_temperature)**2
                
                # Update the Dirichlet side temperature values by performing a relaxation
                interpolated_temperature = (1.0-relaxation_factor)*old_temperature + relaxation_factor*new_temperature
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, interpolated_temperature )

            iteration += 1
            rel_res_norm = (temperature_difference**0.5) / len(self.mapper.str_interface.Nodes)
            self.print_on_rank_zero("::[ConjugateHeatTransferSolver]::", "Iteration: " + str(iteration) + " Relative residual: " + str(rel_res_norm))

            # Check convergence
            if rel_res_norm <= temp_rel_tol:
                self.print_on_rank_zero("::[ConjugateHeatTransferSolver]::", "Converged in " + str(iteration) + " iterations.")
                break
            elif iteration == max_iteration:
                self.print_on_rank_zero("::[ConjugateHeatTransferSolver]::", "Did not converge in " + str(iteration) + " iterations.")

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.fluid_thermal_solver.FinalizeSolutionStep()

    def Solve(self):
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

    def _SetUpDirichletCouplingBoundary(self, model_part):
        # Run the solid interfaces list to fix the temperature DOFs.
        for i_int in range(self.settings["coupling_settings"]["solid_interfaces_list"].size()):
            solid_int_name = self.settings["coupling_settings"]["solid_interfaces_list"][i_int].GetString()
            for node in self.model.GetModelPart(solid_int_name).Nodes:
                node.Fix(KratosMultiphysics.TEMPERATURE)

    def _SetUpMapper(self):
        for i_int in range(self.settings["coupling_settings"]["fluid_interfaces_list"].size()):
            fluid_int_name = self.settings["coupling_settings"]["fluid_interfaces_list"][i_int].GetString()
            for node in self.model.GetModelPart(fluid_int_name).Nodes:
                node.Set(KratosMultiphysics.INTERFACE, True)
        
        for i_int in range(self.settings["coupling_settings"]["solid_interfaces_list"].size()):
            solid_int_name = self.settings["coupling_settings"]["solid_interfaces_list"][i_int].GetString()
            for node in self.model.GetModelPart(solid_int_name).Nodes:
                node.Set(KratosMultiphysics.INTERFACE, True)

        mapper_tolerance = 1e-5
        mapper_max_iteration = 50
        mapper_search_radius_factor = 1.0
        self.mapper = ncosm.NonConformant_OneSideMap(self.fluid_solver.GetComputingModelPart(),
                                                     self.solid_thermal_solver.GetComputingModelPart(), 
                                                     mapper_search_radius_factor, 
                                                     mapper_max_iteration,
                                                     mapper_tolerance)
