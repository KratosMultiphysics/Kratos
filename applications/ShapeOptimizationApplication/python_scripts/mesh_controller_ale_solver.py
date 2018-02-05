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
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import time as timer

from mesh_controller_base import MeshController

# # ==============================================================================
class MeshControllerUsingALESolver( MeshController) :
    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart, MeshMotionSettings ):
        self.OptimizationModelPart = OptimizationModelPart
        self.MeshMotionSettings = MeshMotionSettings

        mesh_solver_module = __import__(self.MeshMotionSettings["mesh_solver_settings"]["solver_type"].GetString())
        self.mesh_solver = mesh_solver_module.CreateSolver(self.OptimizationModelPart, self.MeshMotionSettings["mesh_solver_settings"])
        
        self.mesh_solver.AddVariables()

    # --------------------------------------------------------------------------    
    def Initialize( self ):
        self.mesh_solver.AddDofs()
        self.mesh_solver.Initialize()
        self.mesh_solver.SetEchoLevel(0)    

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable( self, InputVariable ):
        print("\n> Starting to update the mesh...")
        startTime = timer.time()

        # Extract surface nodes
        sub_model_part_of_surface = "surface_nodes"     
        GeometryUtilities(self.OptimizationModelPart).ExtractSurfaceNodes(sub_model_part_of_surface)

        # Apply shape update as boundary condition for computation of mesh displacement 
        for node in self.OptimizationModelPart.GetSubModelPart(sub_model_part_of_surface).Nodes:
            node.Fix(MESH_DISPLACEMENT_X)
            node.Fix(MESH_DISPLACEMENT_Y)
            node.Fix(MESH_DISPLACEMENT_Z)              
            disp = Vector(3)
            disp[0] = node.GetSolutionStepValue(SHAPE_UPDATE_X)
            disp[1] = node.GetSolutionStepValue(SHAPE_UPDATE_Y)
            disp[2] = node.GetSolutionStepValue(SHAPE_UPDATE_Z)
            node.SetSolutionStepValue(MESH_DISPLACEMENT,0,disp)

        # Solve for mesh-update
        self.mesh_solver.Solve()

        # Update reference mesh (Since shape updates are imposed as incremental quantities)
        self.mesh_solver.get_mesh_motion_solver().UpdateReferenceMesh()

        # Log absolute mesh displacement
        for node in self.OptimizationModelPart.Nodes:
            mesh_change = Vector(3)
            mesh_change[0] = node.GetSolutionStepValue(MESH_CHANGE_X) + node.GetSolutionStepValue(MESH_DISPLACEMENT_X)
            mesh_change[1] = node.GetSolutionStepValue(MESH_CHANGE_Y) + node.GetSolutionStepValue(MESH_DISPLACEMENT_Y)
            mesh_change[2] = node.GetSolutionStepValue(MESH_CHANGE_Z) + node.GetSolutionStepValue(MESH_DISPLACEMENT_Z)
            node.SetSolutionStepValue(MESH_CHANGE,0,mesh_change)     

        print("> Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")    

    # # --------------------------------------------------------------------------
    # def AddMeshDerivativeToSurfaceGradient( self, DesignSurface, GradientVariable ):
    
    #     print("\n> Starting to compute contribution of mesh-motion to gradient...")
    #     startTime = timer.time()

    #     # Here we solve the pseudo-elastic mesh-motion system again using modified BCs
    #     # The contributions from the mesh derivatives appear as reaction forces

    #     # Apply rhs (not that setting rhs for nodes on the design surface is actually not necessary but still done for performance reasons)
    #     for node in self.OptimizationModelPart.Nodes:
    #         rhs = node.GetSolutionStepValue(GradientVariable)
    #         print(rhs)
    #         node.SetSolutionStepValue(MESH_RHS,0,rhs)

    #     # Apply dirichlet conditions
    #     for node in DesignSurface.Nodes:
    #         node.Fix(MESH_DISPLACEMENT_X)
    #         node.Fix(MESH_DISPLACEMENT_Y)
    #         node.Fix(MESH_DISPLACEMENT_Z)              
    #         xs = Vector(3)
    #         xs[0] = 0.0
    #         xs[1] = 0.0
    #         xs[2] = 0.0
    #         node.SetSolutionStepValue(MESH_DISPLACEMENT,0,xs)            

    #     # Solve mesh-motion problem with previously modified BCs
    #     self.mesh_solver.Solve()

    #     # Compute and add gradient contribution from mesh motion
    #     for node in DesignSurface.Nodes:
    #         gradient = node.GetSolutionStepValue(GradientVariable)
    #         gradient_contribution = node.GetSolutionStepValue(MESH_REACTION)
    #         print(gradient)
    #         print(gradient_contribution)
    #         modified_gradient = gradient + gradient_contribution 
    #         node.SetSolutionStepValue(GradientVariable,modified_gradient)            

    #     # Reset mesh displacement to zero
    #     for node in self.OptimizationModelPart.Nodes:            
    #         zero_vector = Vector(3)
    #         zero_vector[0] = 0.0
    #         zero_vector[1] = 0.0
    #         zero_vector[2] = 0.0
    #         node.SetSolutionStepValue(MESH_DISPLACEMENT,zero_vector)    

    #     print("> Time needed for computing mesh-motion contribution to gradient = ",round(timer.time() - startTime,2),"s")

# ==============================================================================