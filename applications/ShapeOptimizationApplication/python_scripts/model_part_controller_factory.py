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
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# ==============================================================================
def CreateController(model_settings, model):
    return ModelPartController(model_settings, model)

# ==============================================================================
class ModelPartController:
    # --------------------------------------------------------------------------
    def __init__(self, model_settings, model):
        self.model_settings = model_settings

        default_settings = Parameters("""
        {
            "domain_size"           : 3,
            "model_part_name"       : "OPTIMIZATION_MODEL_PART_NAME",
            "model_import_settings"              : {
                "input_type"     : "mdpa",
                "input_filename" : "OPTIMIZATION_MODEL_PART_FILENAME"
            },
            "design_surface_sub_model_part_name" : "DESIGN_SURFACE_NAME",
            "damping" : {
                "apply_damping"   : false,
                "damping_regions" : []
            },
            "mesh_motion" : {
                "apply_mesh_solver" : false,
                "solver_settings" : { },
                "boundary_conditions_process_list" : []
            }
        }""")

        self.model_settings.ValidateAndAssignDefaults(default_settings)
        self.model_settings["model_import_settings"].ValidateAndAssignDefaults(default_settings["model_import_settings"])

        self.model = model

        model_part_name = self.model_settings["model_part_name"].GetString()

        self.optimization_model_part = ModelPart(model_part_name)
        self.model.AddModelPart(self.optimization_model_part)
        # TODO use this line after model_v3 is merged:
        # self.optimization_model_part = model.CreateModelPart(model_part_name)

        self.optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.model_settings["domain_size"].GetInt())

        mesh_motion_settings = self.model_settings["mesh_motion"]

        if mesh_motion_settings["apply_mesh_solver"].GetBool():
            from mesh_controller_with_solver import MeshControllerWithSolver
            self.mesh_controller = MeshControllerWithSolver(mesh_motion_settings, model)
        else:
            from mesh_controller_basic_updating import MeshControllerBasicUpdating
            self.mesh_controller = MeshControllerBasicUpdating(self.optimization_model_part)

        self._design_surface = None
        self._damping_utility = None

    # --------------------------------------------------------------------------
    def ImportOptimizationModelPart(self):
        input_type = self.model_settings["model_import_settings"]["input_type"].GetString()
        if input_type != "mdpa":
            raise RuntimeError("The model part for the optimization has to be read from the mdpa file!")
        input_filename = self.model_settings["model_import_settings"]["input_filename"].GetString()

        model_part_io = ModelPartIO(input_filename)
        model_part_io.ReadModelPart(self.optimization_model_part)

    # --------------------------------------------------------------------------
    def InitializeMeshController(self):
        self.mesh_controller.Initialize()

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable(self, InputVariable):
        self.mesh_controller.UpdateMeshAccordingInputVariable(InputVariable)

    # --------------------------------------------------------------------------
    def SetMeshToReferenceMesh(self):
        MeshControllerUtilities(self.optimization_model_part).SetMeshToReferenceMesh()

    # --------------------------------------------------------------------------
    def SetReferenceMeshToMesh(self):
        MeshControllerUtilities(self.optimization_model_part).SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def SetDeformationVariablesToZero(self):
        MeshControllerUtilities(self.optimization_model_part).SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def GetOptimizationModelPart(self):
        return self.optimization_model_part

    # --------------------------------------------------------------------------
    def GetModel(self):
        return self.model

    # --------------------------------------------------------------------------
    def GetDesignSurface(self):
        if self._design_surface is None:
            self.__IdentifyDesignSurface()
        return self._design_surface

    # --------------------------------------------------------------------------
    def DampNodalVariableIfSpecified(self, variable):
        if self.model_settings["damping"]["apply_damping"].GetBool():
            self.__GetDampingUtility().DampNodalVariable(variable)

    # --------------------------------------------------------------------------
    def ComputeUnitSurfaceNormals(self):
        GeometryUtilities(self.GetDesignSurface()).ComputeUnitSurfaceNormals()

    # --------------------------------------------------------------------------
    def ProjectNodalVariableOnUnitSurfaceNormals(self, variable):
        GeometryUtilities(self.GetDesignSurface()).ProjectNodalVariableOnUnitSurfaceNormals(variable)

    # --------------------------------------------------------------------------
    def __IdentifyDesignSurface(self):
        nameOfDesignSurface = self.model_settings["design_surface_sub_model_part_name"].GetString()
        if self.optimization_model_part.HasSubModelPart(nameOfDesignSurface):
            self._design_surface = self.optimization_model_part.GetSubModelPart(nameOfDesignSurface)
            print("> The following design surface was defined:\n\n",self._design_surface)
        else:
            raise ValueError("The following sub-model part (design surface) specified for shape optimization does not exist: ",nameOfDesingSurface)

    # --------------------------------------------------------------------------
    def __GetDampingUtility(self):
        if self._damping_utility == None:
            self._damping_utility = DampingUtilities(self.GetDesignSurface(), self.__IdentifyDampingRegions(), self.model_settings["damping"])
        return self._damping_utility

    # --------------------------------------------------------------------------
    def __IdentifyDampingRegions(self):
        print("> The following damping regions are defined: \n")
        damping_regions = {}
        if self.model_settings["damping"]["apply_damping"].GetBool():
            if self.model_settings["damping"].Has("damping_regions"):
                for regionNumber in range(self.model_settings["damping"]["damping_regions"].size()):
                    regionName = self.model_settings["damping"]["damping_regions"][regionNumber]["sub_model_part_name"].GetString()
                    if self.optimization_model_part.HasSubModelPart(regionName):
                        print(regionName)
                        damping_regions[regionName] = self.optimization_model_part.GetSubModelPart(regionName)
                    else:
                        raise ValueError("The following sub-model part specified for damping does not exist: ",regionName)
            else:
                raise ValueError("Definition of damping regions required but not availabe!")
        print("")
        return damping_regions

# ==============================================================================