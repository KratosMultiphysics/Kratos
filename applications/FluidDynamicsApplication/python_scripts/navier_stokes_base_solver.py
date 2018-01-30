from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

#import KratosMultiphysics.MeshingApplication as KratosMeshing

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesBaseSolver(main_model_part, custom_settings)

class NavierStokesBaseSolver(object):

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
            "time_stepping": {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "alpha":-0.1,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "regularization_coef": 1000,
            "MoveMeshFlag": false,
            "use_slip_conditions": false,
            "turbulence_model": "None",
            "use_spalart_allmaras": false,
            "reorder": false
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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)          # Q_VALUE belongs to FluidDynamicsApp.


        print("Base class fluid solver variables added correctly")

    def ImportModelPart(self):
        ## Read model part
        self._ModelPartReading()
        ## Replace elements and conditions, check the input reading and set KINEMATIC_VISCOSITY and DENSITY
        self._ExecuteAfterReading()
        ## Set buffer size
        self._SetBufferSize()

        print ("Base class model reading finished.")

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()

        ## Model part writing
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,self.main_model_part)

        #print("main",self.main_model_part.Nodes[12091])
        #print("compute",self.main_model_part.GetSubModelPart("fluid_computational_model_part").Nodes[12091])
        #err

        print("Base class fluid solver DOFs added correctly.")

    def AdaptMesh(self):
        pass

    def GetMinimumBufferSize(self):
        return 3

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            delta_time = self.EstimateDeltaTimeUtility.EstimateDt()
        # User-defined delta time
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble()

        return delta_time

    def Initialize(self):
        raise Exception("Calling the Navier-Stokes base solver. Please implement the custom Initialize() method of your solver.")

    def SaveRestart(self):
        pass #one should write the restart file here

    def Clear(self):
        (self.solver).Clear()

    def Check(self):
        (self.solver).Check()

    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    def InitializeSolutionStep(self):
        (self.solver).InitializeSolutionStep()

    def Predict(self):
        (self.solver).Predict()

    def SolveSolutionStep(self):
        is_converged = (self.solver).SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        (self.solver).FinalizeSolutionStep()

    def Solve(self):
        (self.solver).Solve()

    def _ModelPartReading(self):
        ## Model part reading
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            ## Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            if(self.settings["reorder"].GetBool()):
                print("******************************************************* REORDERING ********************************************************")
                tmp = KratosMultiphysics.Parameters("{}")
                KratosMultiphysics.ReorderAndOptimizeModelPartProcess(self.main_model_part,tmp).Execute()
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
            # kin_viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
            break

        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        # KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)



    def _SetBufferSize(self):
        ## Set the buffer size
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize(self.GetMinimumBufferSize())

    def _GetAutomaticTimeSteppingUtility(self):
        if (self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility2D(self.computing_model_part,
                                                                     self.settings["time_stepping"])
        else:
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility3D(self.computing_model_part,
                                                                     self.settings["time_stepping"])

        return EstimateDeltaTimeUtility
