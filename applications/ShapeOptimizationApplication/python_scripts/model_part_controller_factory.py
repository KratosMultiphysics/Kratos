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
def CreateController(optimization_settings, optimization_mdpa):
    return ModelPartController(optimization_settings, optimization_mdpa)

# ==============================================================================
class ModelPartController:
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, optimization_mdpa):
        self.optimization_settings = optimization_settings
        self.optimization_mdpa = optimization_mdpa

        self.design_surface = None
        self.damping_regions = None

        mesh_motion_settings = self.optimization_settings["design_variables"]["mesh_motion"]

        if mesh_motion_settings["apply_mesh_solver"].GetBool():
            from mesh_controller_with_solver import MeshControllerWithSolver
            self.mesh_controller = MeshControllerWithSolver(mesh_motion_settings, optimization_mdpa)
        else:
            from mesh_controller_basic_updating import MeshControllerBasicUpdating
            self.mesh_controller = MeshControllerBasicUpdating(optimization_mdpa)

    # --------------------------------------------------------------------------
    def ImportOptimizationModelPart(self):
        if self.__IsOptimizationModelPartAlreadyImported():
            print("> Skipping import of optimization model part as already done by another application. ")
        else:
            model_part_io = ModelPartIO(self.optimization_settings["design_variables"]["optimization_model_part_name"].GetString())
            model_part_io.ReadModelPart(self.optimization_mdpa)

    # --------------------------------------------------------------------------
    def InitializeMeshController(self):
        self.mesh_controller.Initialize()

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable(self, InputVariable):
        self.mesh_controller.UpdateMeshAccordingInputVariable(InputVariable)

    # --------------------------------------------------------------------------
    def SetMeshToReferenceMesh(self):
        MeshControllerUtilities(self.optimization_mdpa).SetMeshToReferenceMesh()

    # --------------------------------------------------------------------------
    def SetReferenceMeshToMesh(self):
        MeshControllerUtilities(self.optimization_mdpa).SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def SetDeformationVariablesToZero(self):
        MeshControllerUtilities(self.optimization_mdpa).SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def GetOptimizationModelPart(self):
        return self.optimization_mdpa

    # --------------------------------------------------------------------------
    def GetDesignSurface(self):
        if self.design_surface is None:
            self.__IdentifyDesignSurface()
        return self.design_surface

    # --------------------------------------------------------------------------
    def GetDampingRegions(self):
        if self.damping_regions is None:
            self.__IdentifyDampingRegions()
        return self.damping_regions

    # --------------------------------------------------------------------------
    def __IsOptimizationModelPartAlreadyImported(self):
        if self.optimization_mdpa.NumberOfNodes()>0:
            self.__CheckIfDomainSizeIsSet()
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __CheckIfDomainSizeIsSet(self):
        if self.optimization_mdpa.ProcessInfo.GetValue(DOMAIN_SIZE) == 0:
            raise ValueError("DOMAIN_SIZE not specified for given optimization model part!")

    # --------------------------------------------------------------------------
    def __IdentifyDesignSurface(self):
        nameOfDesingSurface = self.optimization_settings["design_variables"]["design_surface_sub_model_part_name"].GetString()
        if self.optimization_mdpa.HasSubModelPart(nameOfDesingSurface):
            self.design_surface = self.optimization_mdpa.GetSubModelPart(nameOfDesingSurface)
            print("> The following design surface was defined:\n\n",self.design_surface)
        else:
            raise ValueError("The following sub-model part (design surface) specified for shape optimization does not exist: ",nameOfDesingSurface)

    # --------------------------------------------------------------------------
    def __IdentifyDampingRegions(self):
        print("> The following damping regions are defined: \n")
        self.damping_regions = {}
        if self.optimization_settings["design_variables"]["damping"]["perform_damping"].GetBool():
            if self.optimization_settings["design_variables"]["damping"].Has("damping_regions"):
                for regionNumber in range(self.optimization_settings["design_variables"]["damping"]["damping_regions"].size()):
                    regionName = self.optimization_settings["design_variables"]["damping"]["damping_regions"][regionNumber]["sub_model_part_name"].GetString()
                    if self.optimization_mdpa.HasSubModelPart(regionName):
                        print(regionName)
                        self.damping_regions[regionName] = self.optimization_mdpa.GetSubModelPart(regionName)
                    else:
                        raise ValueError("The following sub-model part specified for damping does not exist: ",regionName)
            else:
                raise ValueError("Definition of damping regions required but not availabe!")
        print("")

# ==============================================================================