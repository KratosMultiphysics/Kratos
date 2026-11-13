# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================


# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# ==============================================================================
def CreateController(model_settings, model):
    return ModelPartController(model_settings, model)

# ==============================================================================
class ModelPartController:
    # --------------------------------------------------------------------------
    def __init__(self, model_settings, model):
        self.model_settings = model_settings

        default_settings = KM.Parameters("""
        {
            "domain_size"           : 3,
            "model_part_name"       : "OPTIMIZATION_MODEL_PART_NAME",
            "model_import_settings"              : {
                "input_type"     : "mdpa",
                "input_filename" : "OPTIMIZATION_MODEL_PART_FILENAME"
            },
            "design_surface_sub_model_part_name" : "DESIGN_SURFACE_NAME",
            "damping" : {
                "apply_damping"      : false,
                "recalculate_damping": false,
                "max_neighbor_nodes" : 10000,
                "damping_regions"    : []
            },
            "direction_damping" : {
                "recalculate_damping": true,
                "max_neighbor_nodes" : 10000,
                "damping_regions"    : []
            },
            "mesh_motion" : {
                "apply_mesh_solver" : false
            },
            "write_iteration_restart_files": false            
        }""")

        self.model_settings.ValidateAndAssignDefaults(default_settings)
        self.model_settings["model_import_settings"].ValidateAndAssignDefaults(default_settings["model_import_settings"])
        self.model_settings["damping"].ValidateAndAssignDefaults(default_settings["damping"])
        self.model_settings["direction_damping"].ValidateAndAssignDefaults(default_settings["direction_damping"])

        for direction_damping_settings in self.model_settings["direction_damping"]["damping_regions"].values():
            if not direction_damping_settings.Has("max_neighbor_nodes"):
                max_neighbors = self.model_settings["direction_damping"]["max_neighbor_nodes"].GetInt()
                direction_damping_settings.AddEmptyValue("max_neighbor_nodes").SetInt(max_neighbors)

        self.model = model

        model_part_name = self.model_settings["model_part_name"].GetString()
        self.optimization_model_part = model.CreateModelPart(model_part_name)
        self.optimization_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, self.model_settings["domain_size"].GetInt())

        if self.model_settings["mesh_motion"]["apply_mesh_solver"].GetBool():
            from KratosMultiphysics.ShapeOptimizationApplication.mesh_controllers.mesh_controller_with_solver import MeshControllerWithSolver
            self.mesh_controller = MeshControllerWithSolver(self.model_settings["mesh_motion"], self.model)
        else:
            from KratosMultiphysics.ShapeOptimizationApplication.mesh_controllers.mesh_controller_basic_updating import MeshControllerBasicUpdating
            self.mesh_controller = MeshControllerBasicUpdating(self.optimization_model_part)

        self.design_surface = None
        self.damping_utility = None
        self.direction_dampings = []
        self.is_iteration_restart_files_written = self.model_settings["write_iteration_restart_files"].GetBool()

    # --------------------------------------------------------------------------
    def Initialize(self):
        self.__ImportOptimizationModelPart()
        self.__IdentifyDesignSurface()

        self.mesh_controller.Initialize()

    def InitializeDamping(self):
        """Initialize damping utilities, should be called after mapper is initialized"""
        if self.model_settings["damping"]["apply_damping"].GetBool():
            self.damping_utility = KSO.DampingUtilities(
                self.design_surface, self.model_settings["damping"]
            )

        for direction_damping_settings in self.model_settings["direction_damping"]["damping_regions"].values():
            self.direction_dampings.append(
                KSO.DirectionDampingUtilities(
                    self.design_surface, direction_damping_settings
                )
            )

    # --------------------------------------------------------------------------
    def SetMinimalBufferSize(self, buffer_size):
        if self.optimization_model_part.GetBufferSize() < buffer_size:
            self.optimization_model_part.SetBufferSize(buffer_size)

    # --------------------------------------------------------------------------
    def UpdateTimeStep(self, step):
        self.optimization_model_part.CloneTimeStep(step)
        self.optimization_model_part.ProcessInfo.SetValue(KM.STEP, step)

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable(self, InputVariable):
        self.mesh_controller.UpdateMeshAccordingInputVariable(InputVariable)

        if self.model_settings["damping"]["recalculate_damping"].GetBool():
            self.damping_utility = KSO.DampingUtilities(
                self.design_surface, self.model_settings["damping"]
            )

        if self.model_settings["direction_damping"]["recalculate_damping"].GetBool():
            self.direction_dampings = []
            for direction_damping_settings in self.model_settings["direction_damping"]["damping_regions"].values():
                self.direction_dampings.append(
                    KSO.DirectionDampingUtilities(
                        self.design_surface, direction_damping_settings
                    )
                )

    # --------------------------------------------------------------------------
    def SetMeshToReferenceMesh(self):
        KSO.MeshControllerUtilities(self.optimization_model_part).SetMeshToReferenceMesh()

    # --------------------------------------------------------------------------
    def SetReferenceMeshToMesh(self):
        KSO.MeshControllerUtilities(self.optimization_model_part).SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def SetDeformationVariablesToZero(self):
        KSO.MeshControllerUtilities(self.optimization_model_part).SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def GetOptimizationModelPart(self):
        return self.optimization_model_part

    # --------------------------------------------------------------------------
    def GetModel(self):
        return self.model

    # --------------------------------------------------------------------------
    def GetDesignSurface(self):
        return self.design_surface

    # --------------------------------------------------------------------------
    def IsIterationRestartFilesWritten(self):
        return self.is_iteration_restart_files_written        

    # --------------------------------------------------------------------------
    def DampNodalSensitivityVariableIfSpecified(self, variable):
        if self.model_settings["damping"]["apply_damping"].GetBool():
            self.damping_utility.DampNodalVariable(variable)

        for direction_damping in reversed(self.direction_dampings):
            direction_damping.DampNodalVariable(variable)

    # --------------------------------------------------------------------------
    def DampNodalUpdateVariableIfSpecified(self, variable):
        for direction_damping in self.direction_dampings:
            direction_damping.DampNodalVariable(variable)

        if self.model_settings["damping"]["apply_damping"].GetBool():
            self.damping_utility.DampNodalVariable(variable)

    # --------------------------------------------------------------------------
    def ComputeUnitSurfaceNormals(self):
        KSO.GeometryUtilities(self.GetDesignSurface()).ComputeUnitSurfaceNormals()

    # --------------------------------------------------------------------------
    def ProjectNodalVariableOnUnitSurfaceNormals(self, variable):
        KSO.GeometryUtilities(self.GetDesignSurface()).ProjectNodalVariableOnUnitSurfaceNormals(variable)

    # --------------------------------------------------------------------------
    def __ImportOptimizationModelPart(self):
        input_type = self.model_settings["model_import_settings"]["input_type"].GetString()
        if input_type != "mdpa":
            raise RuntimeError("The model part for the optimization has to be read from the mdpa file!")
        input_filename = self.model_settings["model_import_settings"]["input_filename"].GetString()

        model_part_io = KM.ModelPartIO(input_filename)
        model_part_io.ReadModelPart(self.optimization_model_part)

        self.SetMinimalBufferSize(1)

    # --------------------------------------------------------------------------
    def __IdentifyDesignSurface(self):
        nameOfDesignSurface = self.model_settings["design_surface_sub_model_part_name"].GetString()
        if self.optimization_model_part.HasSubModelPart(nameOfDesignSurface):
            self.design_surface = self.optimization_model_part.GetSubModelPart(nameOfDesignSurface)
            KM.Logger.Print("")
            KM.Logger.PrintInfo("ShapeOpt", "The following design surface was defined:\n\n",self.design_surface)
        else:
            raise ValueError("The following sub-model part (design surface) specified for shape optimization does not exist: ",nameOfDesignSurface)

# ==============================================================================
