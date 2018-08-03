from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")
KratosMultiphysics.CheckRegisteredApplications("ConvectionDiffusionApplication")


# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff
import KratosMultiphysics.FSIApplication as FSIApplication
import NonConformant_OneSideMap as ncosm

def CreateSolver(main_model_part, custom_settings):
    
    return ConjugateHeatTransferSolver(main_model_part, custom_settings)

class ConjugateHeatTransferSolver(object):
    
    def __init__(self, model, custom_settings):
        self.model = model
        
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "ConjugateHeatTransfer",
            "domain_size": -1,
            "ThermallyCoupled": {
                "fluid_solver_settings": {
                    "solver_type": "Monolithic",
                    "domain_size": 3,
                    "model_import_settings": {
                        "input_type": "mdpa",
                        "input_filename": "Test"
                    },
                    "echo_level": 0,
                    "compute_reactions": false,
                    "dynamic_tau": 0.0,
                    "oss_switch": 0,
                    "maximum_iterations": 10,
                    "relative_velocity_tolerance": 0.001,
                    "absolute_velocity_tolerance": 0.00001,
                    "relative_pressure_tolerance": 0.001,
                    "absolute_pressure_tolerance": 0.00001,
                    "volume_model_part_name": "Parts_Parts_Auto1",
                    "skin_parts": [
                        "Outlet3D_Outlet_pressure_Auto1",
                        "NoSlip3D_No_Slip_Auto1"
                    ],
                    "no_skin_parts": [],
                    "time_stepping": {
                        "automatic_time_step": false,
                        "time_step": 0.01
                    }
                },
                "thermal_solver_settings": {
                    "domain_size": 3,
                    "echo_level": 0,
                    "solver_type": "Transient",
                    "analysis_type": "linear",
                    "model_import_settings": {
                        "input_type": "use_input_model_part",
                        "input_filename": "unknown_name"
                    },
                    "transient_parameters": {
                        "dynamic_tau": 1.0,
                        "theta": 0.5
                    },
                    "computing_model_part_name": "Thermal",
                    "material_import_settings": {
                        "materials_filename": "ThermicMaterials.json"
                    },
                    "convection_diffusion_variables": {
                        "density_variable": "DENSITY",
                        "diffusion_variable": "CONDUCTIVITY",
                        "unknown_variable": "TEMPERATURE",
                        "volume_source_variable": "HEAT_FLUX",
                        "surface_source_variable": "FACE_HEAT_FLUX",
                        "projection_variable": "PROJECTED_SCALAR1",
                        "convection_variable": "CONVECTION_VELOCITY",
                        "mesh_velocity_variable": "MESH_VELOCITY",
                        "transfer_coefficient_variable": "",
                        "velocity_variable": "VELOCITY",
                        "specific_heat_variable": "SPECIFIC_HEAT",
                        "reaction_variable": "REACTION_FLUX"
                    }
                }
            },
            "HeatTransfer": {
                "solid_solver_settings": {

                 },
                "thermal_solver_settings": {
                    "domain_size": 3,
                    "echo_level": 0,
                    "solver_type": "Transient",
                    "analysis_type": "linear",
                    "model_import_settings": {
                        "input_type": "use_input_model_part",
                        "input_filename": "unknown_name"
                    },
                    "transient_parameters": {
                        "dynamic_tau": 1.0,
                        "theta": 0.5
                    },
                    "computing_model_part_name": "Thermal",
                    "material_import_settings": {
                        "materials_filename": "ThermicMaterials.json"
                    },
                    "convection_diffusion_variables": {
                        "density_variable": "DENSITY",
                        "diffusion_variable": "CONDUCTIVITY",
                        "unknown_variable": "TEMPERATURE",
                        "volume_source_variable": "HEAT_FLUX",
                        "surface_source_variable": "FACE_HEAT_FLUX",
                        "projection_variable": "PROJECTED_SCALAR1",
                        "convection_variable": "CONVECTION_VELOCITY",
                        "mesh_velocity_variable": "MESH_VELOCITY",
                        "transfer_coefficient_variable": "",
                        "velocity_variable": "VELOCITY",
                        "specific_heat_variable": "SPECIFIC_HEAT",
                        "reaction_variable": "REACTION_FLUX"
                    }
                }
            }
        }

        """)
        
        ## Overwrite the default settings with user-provided parameters

        self.settings = custom_settings

        
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        domain_size = self.settings["domain_size"].GetInt()
        
        if(self.settings["ThermallyCoupled"]["fluid_solver_settings"]["domain_size"].GetInt() != domain_size):
            raise Exception("domain size for the fluid solver is not consistent")
        
        if(self.settings["ThermallyCoupled"]["thermal_solver_settings"]["domain_size"].GetInt() != domain_size):
            raise Exception("domain size for the fluid solver is not consistent")
        self.domain_size = domain_size
         
        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        
        import python_solvers_wrapper_fluid

        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["ThermallyCoupled"]["fluid_solver_settings"],"OpenMP")
        
        self.main_model_part = self.fluid_solver.main_model_part




        
        
        import python_solvers_wrapper_convection_diffusion
        
        self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model,custom_settings["ThermallyCoupled"]["thermal_solver_settings"],"OpenMP")
           
                

        self.solid_thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model,custom_settings["HeatTransfer"]["thermal_solver_settings"],"OpenMP")

        settings_aux=custom_settings["HeatTransfer"]["thermal_solver_settings"]

             
        #adding the model part name for the solid part        
        #settings_aux.AddEmptyValue("model_part_name")
        #settings_aux["model_part_name"].SetString("SolidModelPart")
        #settings_aux.AddEmptyValue("volume_model_part_name")
        #settings_aux["volume_model_part_name"].SetString("Parts_solids")
        
    def AddVariables(self):
        self.fluid_solver.AddVariables()
        self.thermal_solver.AddVariables() 
        
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_PAUX)
        self.thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_PAUX)
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_PAUX)
        
        # Temporary container for un-relaxed temperature
        
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part)
        
        self.solid_thermal_solver.AddVariables()
        self.solid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_PAUX)

    def ImportModelPart(self):
        
        self.fluid_solver.ImportModelPart()
        
        self.solid_thermal_solver.ImportModelPart()
        
        if not self.model.HasModelPart(self.fluid_solver.main_model_part.Name):
            self.model.AddModelPart(self.fluid_solver.main_model_part)
        if not self.model.HasModelPart(self.solid_thermal_solver.main_model_part.Name):
            self.model.AddModelPart(self.solid_thermal_solver.main_model_part)

        #here cloning the fluid modelpart to thermal_model_part so that the nodes are shared
        convection_diffusion_settings = self.thermal_solver.GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()


        if(self.domain_size == 2):
            modeler.GenerateModelPart(self.main_model_part, self.thermal_solver.GetComputingModelPart(), "Element2D3N", "LineCondition2D2N")
        else:
            modeler.GenerateModelPart(self.main_model_part, self.thermal_solver.GetComputingModelPart(), "Element3D4N", "SurfaceCondition3D3N")

        self.thermal_solver.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)

        self.thermal_solver.ImportModelPart()

        count_nodes = 0
        for node in self.main_model_part.Nodes:
            count_nodes += 1
            node.Id = count_nodes
        count_elems = 0
        for elem in self.main_model_part.Elements:
            count_elems += 1
            elem.Id = count_elems
        count_conds = 0
        for cond in self.main_model_part.Conditions:
            count_conds += 1
            cond.Id = count_conds

        for node in self.solid_thermal_solver.main_model_part.Nodes:
            count_nodes += 1
            node.Id = count_nodes
        for elem in self.solid_thermal_solver.main_model_part.Elements:
            count_elems += 1
            elem.Id = count_elems
        for cond in self.solid_thermal_solver.main_model_part.Conditions:
            count_conds += 1
            cond.Id = count_conds

    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.thermal_solver.AddDofs()
        self.solid_thermal_solver.AddDofs()
        print(self.solid_thermal_solver)
        

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        
        return self.fluid_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
                
        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()

        return max(buffer_size_fluid, buffer_size_thermal)

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.thermal_solver.Initialize()
        self.solid_thermal_solver.Initialize()

    def SaveRestart(self):
        pass #one should write the restart file here

    def Clear(self):
        (self.fluid_solver).Clear()
        (self.thermal_solver).Clear()
        (self.solid_thermal_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.thermal_solver).Check()
        (self.solid_thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)
        (self.solid_thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt

        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

       
        self.solid_thermal_solver.main_model_part.CloneTimeStep(new_time)
        self.solid_thermal_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        return new_time

    def PrepareModelPart(self):
          
        self.fluid_solver.PrepareModelPart()
        self.thermal_solver.PrepareModelPart()
        self.solid_thermal_solver.PrepareModelPart()

    def InitializeSolutionStep(self):
        self.fluid_solver.InitializeSolutionStep()
        
        self.thermal_solver.InitializeSolutionStep()

        self.solid_thermal_solver.InitializeSolutionStep()

        
    def Predict(self):
        self.fluid_solver.Predict()
        self.thermal_solver.Predict()
        self.solid_thermal_solver.Predict()
        


        for node in self.fluid_solver.GetComputingModelPart().Nodes:
            temperature=node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)	
            gravity_y=-10.0*(1-0.001*(temperature-0.0))
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y,gravity_y)	
 
    
    def SolveSolutionStep(self):
        self.fluid_solver.SolveSolutionStep()
        
        for node in self.solid_thermal_solver.GetComputingModelPart().Nodes:
            if(node.X<-0.249999999999):
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,300.15)	
                node.Fix(KratosMultiphysics.TEMPERATURE)

        for node in self.fluid_solver.GetComputingModelPart().Nodes:
            if(node.X>.999999):
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,595.15)	
                node.Fix(KratosMultiphysics.TEMPERATURE)


        num_coupling_iterations = 10
        temperature_relaxation_factor = 0.7
        coupling_relative_tolerance = 1e-5

        self.setUpMapper()
        

        self.setUpDirichletCouplingBoundary(self.solid_thermal_solver.GetComputingModelPart())
        

        iter = 0
        while iter < num_coupling_iterations:

            # Solve Dirichlet side -> Get reactions
            
            self.solid_thermal_solver.Solve()
            # Map reactions

            FSIApplication.VariableRedistributionUtility.DistributePointValues( self.mapper.str_interface, KratosMultiphysics.REACTION_FLUX, KratosMultiphysics.NODAL_PAUX, 1e-5, 50)
                
            self.mapper.StructureToFluid_ScalarMap(KratosMultiphysics.NODAL_PAUX,KratosMultiphysics.FACE_HEAT_FLUX,False)
            
            # Solve Neumann side
            self.thermal_solver.SolveSolutionStep()
        #   # Get updated temperature


            self.mapper.FluidToStructure_ScalarMap(KratosMultiphysics.TEMPERATURE,KratosMultiphysics.NODAL_PAUX,True)
            temperature_difference = 0.0
            for node in self.mapper.str_interface.Nodes:
                #self.fl_interface
                old_temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                new_temperature = node.GetSolutionStepValue(KratosMultiphysics.NODAL_PAUX)
                interpolated_temperature = (1.0-temperature_relaxation_factor)*old_temperature + temperature_relaxation_factor*new_temperature
                temperature_difference += (old_temperature-new_temperature)**2
                
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, interpolated_temperature )

            iter += 1
            if (temperature_difference**0.5)/len(self.mapper.str_interface.Nodes) <= coupling_relative_tolerance:
                
                break




    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

    def Solve(self):
        
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

    def setUpDirichletCouplingBoundary(self,model_part):
        for cond in model_part.Conditions:
            
            for node in cond.GetNodes():
                node.Fix(KratosMultiphysics.TEMPERATURE)

    def setUpMapper(self):

        for cond in self.fluid_solver.GetComputingModelPart().Conditions:
            for node in cond.GetNodes():
                node.Set(KratosMultiphysics.INTERFACE,True)
                

        for cond in self.solid_thermal_solver.GetComputingModelPart().Conditions:
            
            for node in cond.GetNodes():
                node.Set(KratosMultiphysics.INTERFACE,True)

        #for node in self.fluid_solver.GetComputingModelPart().Nodes:
        #    if(node.X>0.0001):
        #        node.Set(KratosMultiphysics.INTERFACE,False)	

        self.mapper = ncosm.NonConformant_OneSideMap(self.fluid_solver.GetComputingModelPart(),self.solid_thermal_solver.GetComputingModelPart(), search_radius_factor=0.0, it_max=50, tol=1e-5)
        