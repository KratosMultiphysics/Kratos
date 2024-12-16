
# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.DelaunayMeshingApplication  as KratosDelaunay


# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    return CoupledPfemFluidThermalSolver(main_model_part, custom_settings)

class CoupledPfemFluidThermalSolver(PythonSolver):

    def __init__(self, model, custom_settings):

        self._validate_settings_in_baseclass = True

        super(CoupledPfemFluidThermalSolver, self).__init__(model, custom_settings)

        # Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        from KratosMultiphysics.PfemFluidDynamicsApplication import python_solvers_wrapper_pfem_fluid
        self.fluid_solver = python_solvers_wrapper_pfem_fluid.CreateSolverByParameters(self.model,self.settings["fluid_solver_settings"],"OpenMP")

        from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion
        self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model,self.settings["thermal_solver_settings"],"OpenMP")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver_type": "coupled_pfem_fluid_thermal_solver",
            "model_part_name": "PfemFluidModelPart",
            "time_stepping"               : {
                    "automatic_time_step" : false,
                    "time_step"           : 0.001
                },
                "domain_size": 2,
            "echo_level"                         : 1,
            "fluid_solver_settings":{
                "physics_type"   : "fluid",
                "model_import_settings":{
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
                },
                "buffer_size": 3,
                "echo_level": 1,
                "reform_dofs_at_each_step": false,
                "clear_storage": false,
                "compute_reactions": true,
                "move_mesh_flag": true,
                "dofs"                : [],
                "stabilization_factor": 1.0,
                "line_search": false,
                "compute_contact_forces": false,
                "block_builder": false,
                "component_wise": false,
                "predictor_corrector": true,
                "time_order": 2,
                "maximum_velocity_iterations": 1,
                "maximum_pressure_iterations": 7,
                "velocity_tolerance": 1e-5,
                "pressure_tolerance": 1e-5,
                "pressure_linear_solver_settings":  {
                    "solver_type"                    : "amgcl",
                    "max_iteration"                  : 5000,
                    "tolerance"                      : 1e-9,
                    "provide_coordinates"            : false,
                    "scaling"                        : false,
                    "smoother_type"                  : "damped_jacobi",
                    "krylov_type"                    : "cg",
                    "coarsening_type"                : "aggregation",
                    "verbosity"                      : 0
                },
                "velocity_linear_solver_settings": {
                    "solver_type"                    : "bicgstab",
                    "max_iteration"                  : 5000,
                    "tolerance"                      : 1e-9,
                    "preconditioner_type"            : "none",
                    "scaling"                        : false
                },
                "solving_strategy_settings":{
                   "time_step_prediction_level": 0,
                   "max_delta_time": 1.0e-5,
                   "fraction_delta_time": 0.9,
                   "rayleigh_damping": false,
                   "rayleigh_alpha": 0.0,
                   "rayleigh_beta" : 0.0
                },
                "bodies_list": [],
                "problem_domain_sub_model_part_list": [],
                "constitutive_laws_list": [],
                "processes_sub_model_part_list": [],
                "constraints_process_list": [],
                "loads_process_list"       : [],
                "output_process_list"      : [],
                "output_configuration"     : {},
                "problem_process_list"     : [],
                "processes"                : {},
                "output_processes"         : {},
                "check_process_list": []
            },
            "thermal_solver_settings": {
                "solver_type": "Transient",
                "analysis_type": "linear",
                "computing_model_part_name": "thermal_computing_domain",
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                        "materials_filename": "ThermalMaterials.json"
                },
                "reform_dofs_at_each_step": true
            },
            "coupling_settings": {}
        }""")

        this_defaults.AddMissingParameters(super(CoupledPfemFluidThermalSolver, cls).GetDefaultParameters())

        return this_defaults

    def AddVariables(self):
        # Import the fluid and thermal solver variables. Then merge them to have them in both pfem and thermal solvers.
        self.fluid_solver.AddVariables()
        self.AddMaterialVariables()
        self.AddPfemVariables()
        self.thermal_solver.AddVariables()
        KM.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part)
        print("::[PfemFluidThermallyCoupledSolver]:: Variables MERGED")

    def ImportModelPart(self):
        # Call the fluid solver to import the model part from the mdpa
        self.fluid_solver._ImportModelPart(self.fluid_solver.main_model_part,self.settings["fluid_solver_settings"]["model_import_settings"])

    def CloneThermalModelPart(self):
        # Save the convection diffusion settings
        convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KM.CONVECTION_DIFFUSION_SETTINGS)

        # Here the pfem model part is cloned to be the thermal model part so that the nodes are shared
        modeler = KM.ConnectivityPreserveModeler()
        if self.domain_size == 2:
            modeler.GenerateModelPart(self.fluid_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "EulerianConvDiff2D",
                                      "ThermalFace2D2N")
        else:
            modeler.GenerateModelPart(self.fluid_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "EulerianConvDiff3D",
                                      "ThermalFace3D3N")

        # Set the saved convection diffusion settings to the new thermal model part
        self.thermal_solver.main_model_part.ProcessInfo.SetValue(KM.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)

    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.thermal_solver.AddDofs()

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()

    def ComputeDeltaTime(self):
        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal)

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.thermal_solver.Initialize()

    def Clear(self):
        (self.fluid_solver).Clear()
        (self.thermal_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.fluid_solver.AdvanceInTime(current_time)
        return new_time

    def InitializeSolutionStep(self):
        KratosPfemFluid.UpdateThermalModelPartProcess(
            self.fluid_solver.main_model_part, \
            self.thermal_solver.main_model_part, \
            self.thermal_solver.GetComputingModelPart(), \
            self.domain_size).Execute()
        self.fluid_solver.InitializeSolutionStep()
        self.thermal_solver.InitializeSolutionStep()

    def Predict(self):
        self.fluid_solver.Predict()
        self.thermal_solver.Predict()

    def SolveSolutionStep(self):
        pfem_is_converged = self.fluid_solver.SolveSolutionStep()
        KratosPfemFluid.SetMeshVelocityForThermalCouplingProcess(self.fluid_solver.main_model_part).Execute()
        KratosPfemFluid.SetMaterialPropertiesForThermalCouplingProcess(self.fluid_solver.main_model_part,self.thermal_solver.main_model_part).Execute()
        KM.Logger.Print("\nSolution of convection-diffusion at t=" + "{:.3f}".format(self.fluid_solver.main_model_part.ProcessInfo[KM.TIME]) + "s", label="")
        KM.Logger.Flush()
        thermal_is_converged = self.thermal_solver.SolveSolutionStep()
        return (pfem_is_converged and thermal_is_converged)

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

    def PrepareModelPart(self):
        self.CloneThermalModelPart()
        self.PrepareThermalModelPart()
        self.fluid_solver.PrepareModelPart()

    def PrepareThermalModelPart(self):
        # Thermal model part is being prepared here instead of calling its application
        # because there is a small change needed in that function (self.thermal_solver.PrepareModelPart()).
        # ATTENTION: in the future, better call self.thermal_solver.PrepareModelPart(), but sending a flag
        # to not execute TetrahedralMeshOrientationCheck inside that function.
        if not self.thermal_solver.is_restarted():
            materials_imported = self.thermal_solver.import_materials()
            if materials_imported:
                KM.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Materials were successfully imported.")
            else:
                KM.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Materials were not imported.")

            KM.ReplaceElementsAndConditionsProcess(self.thermal_solver.main_model_part,self.thermal_solver._get_element_condition_replace_settings()).Execute()
            self.thermal_solver._set_and_fill_buffer()

        if (self.thermal_solver.settings["echo_level"].GetInt() > 0):
            KM.Logger.PrintInfo(self.thermal_solver.model)

        KM.Logger.PrintInfo("::[ConvectionDiffusionSolver]::", "ModelPart prepared for Solver.")

    def AddMaterialVariables(self):

        if self.settings.Has("fluid_solver_settings"):
            if self.settings["fluid_solver_settings"].Has("constitutive_laws_list"):
                self.constitutive_laws_names     = self.settings["fluid_solver_settings"]["constitutive_laws_list"]

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.DENSITY):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.DENSITY)
        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.BULK_MODULUS):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.BULK_MODULUS)
        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.DYNAMIC_VISCOSITY):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.DYNAMIC_VISCOSITY)

        for i in range(self.constitutive_laws_names.size()):
            if (self.constitutive_laws_names[i].GetString()=="FrictionalViscoplasticTemperatureDependent2DLaw" or self.constitutive_laws_names[i].GetString()=="FrictionalViscoplasticTemperatureDependent3DLaw"):
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.INTERNAL_FRICTION_ANGLE):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.INTERNAL_FRICTION_ANGLE)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.COHESION):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.COHESION)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT)
            elif (self.constitutive_laws_names[i].GetString()=="HypoelasticTemperatureDependent2DLaw" or self.constitutive_laws_names[i].GetString()=="HypoelasticTemperatureDependent3DLaw"):
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.POISSON_RATIO):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.POISSON_RATIO)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.YOUNG_MODULUS):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.YOUNG_MODULUS)
            elif (self.constitutive_laws_names[i].GetString()=="BinghamTemperatureDependent2DLaw" or self.constitutive_laws_names[i].GetString()=="BinghamTemperatureDependent3DLaw"):
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.FLOW_INDEX):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FLOW_INDEX)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.YIELD_SHEAR):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELD_SHEAR)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.YIELDED):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELDED)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT)
            elif (self.constitutive_laws_names[i].GetString()=="TemperatureDependentMuIRheology2DLaw" or self.constitutive_laws_names[i].GetString()=="TemperatureDependentMuIRheology3DLaw"):
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.STATIC_FRICTION):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.STATIC_FRICTION)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.DYNAMIC_FRICTION):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.DYNAMIC_FRICTION)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ZERO):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ZERO)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DIAMETER):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DIAMETER)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DENSITY):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DENSITY)
                if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT):
                    self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT)
            elif (self.constitutive_laws_names[i].GetString()!="None" and self.constitutive_laws_names[i].GetString()!="NewtonianTemperatureDependent2DLaw" and self.constitutive_laws_names[i].GetString()!="NewtonianTemperatureDependent3DLaw"):
                print("ERROR: THE CONSTITUTIVE LAW PROVIDED FOR THIS SUBMODEL PART IS NOT IN THE PFEM FLUID DATABASE")


    def AddPfemVariables(self):
        print("Add Pfem Variables in pfem_fluid_thermally_coupled_analysis")
        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.MESH_VELOCITY):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.MESH_VELOCITY)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.BODY_FORCE):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.BODY_FORCE)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.NODAL_MASS):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_MASS)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.REACTION):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.REACTION)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.NORMAL):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.NORMAL)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.VOLUME_ACCELERATION):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_VELOCITY):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_VELOCITY)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_ACCELERATION):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_ACCELERATION)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.ISOLATED_NODE):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ISOLATED_NODE)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.NODAL_H_WALL):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_H_WALL)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.NODAL_H):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_H)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosDelaunay.SHRINK_FACTOR):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.SHRINK_FACTOR)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KratosDelaunay.PROPERTY_ID):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.PROPERTY_ID)

        if not self.fluid_solver.main_model_part.HasNodalSolutionStepVariable(KM.HEAT_FLUX):
            self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KM.HEAT_FLUX)


