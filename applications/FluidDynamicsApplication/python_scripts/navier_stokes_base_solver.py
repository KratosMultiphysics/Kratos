from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

#import KratosMultiphysics.MeshingApplication as KratosMeshing

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesBaseSolver(main_model_part, custom_settings)

class NavierStokesBaseSolver:

    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        base_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "navier_stokes_base_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "maximum_iterations": 10,
            "dynamic_tau": 0.0,
            "oss_switch": 0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": false,
            "divergence_clearance_steps": 0,
            "reform_dofs_at_each_step": true,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "alpha":-0.1,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "regularization_coef": 1000,
            "MoveMeshFlag": false,
            "use_slip_conditions": false,
            "turbulence_model": "None",
            "use_spalart_allmaras": false
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(base_settings)

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        ## Set the element replace settings
        self.settings.AddEmptyValue("element_replace_settings")
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                "element_name":"Element3D4N",
                "condition_name": "SurfaceCondition3D3N"
                }
                """)
        elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                "element_name":"Element2D3N",
                "condition_name": "SurfaceCondition2D2N"
                }
                """)
        else:
            raise Exception("domain size is not 2 or 3")

        print("Construction of NavierStokesBaseSolver finished")

    def AddVariables(self):
        ## Add base class variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_TAU) # Variable stored in cfd_variables.h
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.OSS_SWITCH)  # Variable stored in cfd_variables.h
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.M)           # Variable stored in cfd_variables.h
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)          # PATCH_INDEX belongs to FluidDynamicsApp.

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMeshing.AUXILIAR_GRADIENT)  
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMeshing.AUXILIAR_HESSIAN)  
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMeshing.ANISOTROPIC_RATIO)  
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMeshing.MMG_METRIC)  

        # TODO: TURBULENCE MODELS ARE NOT ADDED YET
        #~ if config is not None:
            #~ if hasattr(config, "TurbulenceModel"):
                #~ if config.TurbulenceModel == "Spalart-Allmaras":
                    #~ model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY)
                    #~ model_part.AddNodalSolutionStepVariable(MOLECULAR_VISCOSITY)
                    #~ model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
                    #~ model_part.AddNodalSolutionStepVariable(DISTANCE)

        print("Base class fluid solver variables added correctly")

    def ImportModelPart(self):
        ## Read model part
        self._ModelPartReading()
        ## Replace elements and conditions, check the input reading and set KINEMATIC_VISCOSITY and DENSITY
        self._ExecuteAfterReading()
        ## Set buffer size
        self._SetBufferSize()

        # Adding C_SMAGORINSKY
        for elem in self.main_model_part.Elements:
            elem.SetValue(KratosMultiphysics.C_SMAGORINSKY, 0.0)
    
        print ("Base class model reading finished.")

    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE)
            node.AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z)

        print("Base class fluid solver DOFs added correctly.")

        # TODO: TURBULENCE MODELS ARE NOT ADDED YET
        #~ if config is not None:
            #~ if hasattr(config, "TurbulenceModel"):
                #~ if config.TurbulenceModel == "Spalart-Allmaras":
                    #~ for node in model_part.Nodes:
                        #~ node.AddDof(TURBULENT_VISCOSITY)

    def GetMinimumBufferSize(self):
        return 3;

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here

    def SolverInitialize(self):
        (self.solver).Initialize()

    def SolverInitializeSolutionStep(self):
        (self.solver).InitializeSolutionStep()

    def SolverPredict(self):
        (self.solver).Predict()

    def SolverSolveSolutionStep(self):
        (self.solver).SolveSolutionStep()

    def SolverFinalizeSolutionStep(self):
        (self.solver).FinalizeSolutionStep()

    def Solve(self):
        self.SolverInitialize()
        self.SolverInitializeSolutionStep()
        self.SolverPredict()
        self.SolverSolveSolutionStep()
        self.SolverFinalizeSolutionStep()

    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    def Clear(self):
        (self.solver).Clear()

    def Check(self):
        (self.solver).Check()

    def _ModelPartReading(self):
        ## Model part reading
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            ## Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
        else:
            raise Exception("Other input options are not yet implemented.")

    def _ExecuteAfterReading(self):
        ## Replace element and conditions
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])

        import check_and_prepare_model_process_fluid
        check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()

        # Read the KINEMATIC VISCOSITY and DENSITY and we apply it to the nodes
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            kin_viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
            break

        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)

    def _SetBufferSize(self):
        ## Set the buffer size
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )
