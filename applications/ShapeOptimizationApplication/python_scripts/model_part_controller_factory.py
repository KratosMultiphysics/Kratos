# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

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
                "max_neighbor_nodes" : 10000,
                "damping_regions"    : []
            },
            "mesh_motion" : {
                "apply_mesh_solver" : false
            }
        }""")

        self.model_settings.ValidateAndAssignDefaults(default_settings)
        self.model_settings["model_import_settings"].ValidateAndAssignDefaults(default_settings["model_import_settings"])
        self.model_settings["damping"].ValidateAndAssignDefaults(default_settings["damping"])

        self.model = model

        model_part_name = self.model_settings["model_part_name"].GetString()
        self.optimization_model_part = model.CreateModelPart(model_part_name)
        self.optimization_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, self.model_settings["domain_size"].GetInt())

        if self.model_settings["mesh_motion"]["apply_mesh_solver"].GetBool():
            from .mesh_controller_with_solver import MeshControllerWithSolver
            self.mesh_controller = MeshControllerWithSolver(self.model_settings["mesh_motion"], self.model)
        else:
            from .mesh_controller_basic_updating import MeshControllerBasicUpdating
            self.mesh_controller = MeshControllerBasicUpdating(self.optimization_model_part)

        self.design_surface = None
        self.damping_regions = {}
        self.damping_utility = None

    # --------------------------------------------------------------------------
    def Initialize(self):
        self.__ImportOptimizationModelPart()
        self.__IdentifyDesignSurface()

        self.mesh_controller.Initialize()

        if self.model_settings["damping"]["apply_damping"].GetBool():
            self.__IdentifyDampingRegions()
            self.damping_utility = KSO.DampingUtilities(self.design_surface, self.damping_regions, self.model_settings["damping"])

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
    def DampNodalVariableIfSpecified(self, variable):
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

    # --------------------------------------------------------------------------
    def __IdentifyDampingRegions(self):
        KM.Logger.Print("")
        KM.Logger.PrintInfo("ShapeOpt", "The following damping regions are defined: \n")
        if self.model_settings["damping"]["apply_damping"].GetBool():
            if self.model_settings["damping"].Has("damping_regions"):
                for regionNumber in range(self.model_settings["damping"]["damping_regions"].size()):
                    regionName = self.model_settings["damping"]["damping_regions"][regionNumber]["sub_model_part_name"].GetString()
                    if self.optimization_model_part.HasSubModelPart(regionName):
                        KM.Logger.Print(regionName)
                        self.damping_regions[regionName] = self.optimization_model_part.GetSubModelPart(regionName)
                    else:
                        raise ValueError("The following sub-model part specified for damping does not exist: ",regionName)
            else:
                raise ValueError("Definition of damping regions required but not availabe!")
        KM.Logger.Print("")

# ==============================================================================
