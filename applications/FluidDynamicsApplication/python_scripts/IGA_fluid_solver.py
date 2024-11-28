# Importing the Kratos Library
import KratosMultiphysics

import numpy as np
# from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import json
import os

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication import navier_stokes_solver_vmsmonolithic
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def CreateSolver(main_model_part, custom_settings):
    return IGAFluidSolver(main_model_part, custom_settings)


class IGAFluidSolver(navier_stokes_solver_vmsmonolithic.NavierStokesSolverMonolithic):

    # Write all the points of the skin boundary in an external file
    directory = "txt_files"
    file_name = os.path.join(directory, "true_points.txt")
    if os.path.exists(file_name):
        os.remove(file_name)
    name_mdpa_true_boundary = "mdpa_files/Weird_shape3..."
    # name_mdpa_true_boundary = "mdpa_files/sphere" 
    file_mdpa_exists = os.path.isfile(os.path.join(name_mdpa_true_boundary + ".mdpa"))
    
    if (file_mdpa_exists) :
        print('mdpa file of embedded boundary exists!')
        current_model = KratosMultiphysics.Model()
        skin_model_part2 = current_model.CreateModelPart("skin_model_part2")
        skin_model_part2.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        KratosMultiphysics.ModelPartIO(name_mdpa_true_boundary).ReadModelPart(skin_model_part2)

        with open(file_name, 'w') as file:
            for condition in skin_model_part2.Conditions :
                if (condition.GetGeometry().WorkingSpaceDimension() <= 2) :
                    file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y}\n")
                    file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y}\n")
                else :
                    file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y} {condition.GetNodes()[0].Z} \n")
                    file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y} {condition.GetNodes()[1].Z} \n")
                    file.write(f"{condition.GetNodes()[2].X} {condition.GetNodes()[2].Y} {condition.GetNodes()[2].Z} \n")
                    
    

    # Write all the points of the skin boundary in an external file
    directory = "txt_files"
    file_name = os.path.join(directory, "true_points_outer.txt")
    if os.path.exists(file_name):
        os.remove(file_name)
    name_mdpa_true_boundary = "mdpa_files/external_bunny.."
    file_mdpa_exists = os.path.isfile(os.path.join(name_mdpa_true_boundary + ".mdpa"))
    if (file_mdpa_exists) :
        current_model = KratosMultiphysics.Model()
        skin_model_part2_outer = current_model.CreateModelPart("skin_model_part2_outer")
        skin_model_part2_outer.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        KratosMultiphysics.ModelPartIO(name_mdpa_true_boundary).ReadModelPart(skin_model_part2_outer)
        
        with open(file_name, 'w') as file:
            for condition in skin_model_part2_outer.Conditions :
                file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y}\n")
                file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y}\n")


    @classmethod
    def GetDefaultParameters(cls):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "monolithic",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "formulation": {
                "element_type": "vms"
            },
            "maximum_iterations": 10,
            "use_old_stiffness_in_first_iteration" : false,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": true,
            "assign_neighbour_elements_to_conditions": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"        : {
                "solver_type" : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01,
                "time_step"           : 0.0
            },
            "time_scheme":"bossak",
            "alpha":-0.3,
            "velocity_relaxation":0.9,
            "pressure_relaxation":0.9,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "move_mesh_flag": false,
            "line_search_settings": {
                "max_line_search_iterations" : 5,
                "first_alpha_value"          : 0.5,
                "second_alpha_value"         : 1.0,
                "min_alpha"                  : 0.1,
                "max_alpha"                  : 2.0,
                "line_search_tolerance"      : 0.5
            }
        }""")

        default_settings.AddMissingParameters(super(IGAFluidSolver, cls).GetDefaultParameters())
        return default_settings
    
    def __init__(self, main_model_part, custom_settings):
        super().__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    def AddVariables(self):
        super().AddVariables()

    def AddDofs(self):
        super().AddDofs()
        # KratosMultiphysics.VariableUtils().AddDof(KratosConvDiff.SCALAR_LAGRANGE_MULTIPLIER, self.main_model_part)
    
    def _CreateScheme(self):
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # Cases in which the element manages the time integration
        self.element_integrates_in_time = True
        if self.element_integrates_in_time:
            # "Fake" scheme for those cases in where the element manages the time integration
            # It is required to perform the nodal update once the current time step is solved
            scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(
                domain_size,
                domain_size + 1)
            # In case the BDF2 scheme is used inside the element, the BDF time discretization utility is required to update the BDF coefficients
            if (self.settings["time_scheme"].GetString() == "bdf2"):
                time_order = 2
                self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
            else:
                if  (self.settings["time_scheme"].GetString()!= "crank_nicolson"):
                    err_msg = "Requested elemental time scheme \"" + self.settings["time_scheme"].GetString()+ "\" is not available.\n"
                    err_msg += "Available options are: \"bdf2\" and \"crank_nicolson\""
                    raise Exception(err_msg)
        # Cases in which a time scheme manages the time integration
        else:
            # Bossak time integration scheme
            if self.settings["time_scheme"].GetString() == "bossak":
                if self.settings["consider_periodic_conditions"].GetBool() == True:
                    scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                        self.settings["alpha"].GetDouble(),
                        domain_size,
                        KratosCFD.PATCH_INDEX)
                else:
                    scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                        self.settings["alpha"].GetDouble(),
                        self.settings["move_mesh_strategy"].GetInt(),
                        domain_size)
            # BDF2 time integration scheme
            elif self.settings["time_scheme"].GetString() == "bdf2":
                scheme = KratosCFD.BDF2TurbulentScheme()
            # Time scheme for steady state fluid solver
            elif self.settings["time_scheme"].GetString() == "steady":
                scheme = KratosCFD.ResidualBasedSimpleSteadyScheme(
                        self.settings["velocity_relaxation"].GetDouble(),
                        self.settings["pressure_relaxation"].GetDouble(),
                        domain_size)
            else:
                if  (self.settings["time_scheme"].GetString()!= "crank_nicolson"):
                    err_msg = "Requested time scheme " + self.settings["time_scheme"].GetString() + " is not available.\n"
                    err_msg += "Available options are: \"bossak\", \"bdf2\" ,\"steady\" and \"crank_nicolson\""
                    raise Exception(err_msg)

        return scheme

    def Initialize(self):
        super().Initialize()
        # Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME_INTEGRATION_THETA] = 1.0

        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.USE_CONSTITUTIVE_LAW] = False
        

        
    def InitializeSolutionStep(self):
        
        current_time = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]

        main_model_part = self.GetComputingModelPart()
        for node in main_model_part.Nodes :
            # if (node.X == 0):

            # if (node.X == 0.0 and node.Y == 2.0) or (node.X == 2.0 and node.Y == 2.0):
            #     node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X , 0 , 1.0)
            #     node.Fix(KratosMultiphysics.VELOCITY_X)

            if (node.X == 0.0 and node.Y == 0.0) or (node.X == 2.0 and node.Y == 0.0):
                # node.SetSolutionStepValue(KratosMultiphysics.PRESSURE ,0 , 2*np.pi*(np.cos(2*np.pi*node.Y)-np.cos(2*np.pi*node.X)))
                # node.SetSolutionStepValue(KratosMultiphysics.PRESSURE , 0 , node.X*node.X + node.Y*node.Y)
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE , 0 , 0.0)

                # node.SetSolutionStepValue(KratosMultiphysics.PRESSURE , 0 , 2*np.pi*(np.cos(2*np.pi*node.Y)-np.cos(2*np.pi*node.X)))
                # node.SetSolutionStepValue(KratosMultiphysics.PRESSURE , 0 , ((node.X)**2 + (node.Y)**2)*np.cos(current_time))
                # node.SetSolutionStepValue(KratosMultiphysics.PRESSURE , 0 , ((node.X)**2 + (node.Y)**2) * (current_time))
                # node.SetSolutionStepValue(KratosMultiphysics.PRESSURE , 0 , (node.X + node.Y) * (current_time))
                # node.SetSolutionStepValue(KratosMultiphysics.PRESSURE , 0 , (node.X + node.Y))
                node.Fix(KratosMultiphysics.PRESSURE)
        
        super().InitializeSolutionStep()

        
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        # self.printDofsAndCPs() # Print control points and dofs

    def Finalize(self):
        super().Finalize()

    def _CreateSolutionStrategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            solution_strategy = self._CreateLinearStrategy()
        elif analysis_type == "non_linear":
            solution_strategy = self._CreateNewtonRaphsonStrategy()
        elif analysis_type == "line_search":
            solution_strategy = self._CreateLineSearchStrategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return solution_strategy

    def _CreateLineSearchStrategy(self):
        linesearch_settings = self.settings["line_search_settings"]
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        linesearch_settings.AddValue("max_iteration", self.settings["maximum_iterations"])
        linesearch_settings.AddValue("compute_reactions", self.settings["compute_reactions"])
        linesearch_settings.AddValue("reform_dofs_at_each_step", self.settings["reform_dofs_at_each_step"])
        linesearch_settings.AddValue("move_mesh_flag", self.settings["move_mesh_flag"])
        return KratosMultiphysics.LineSearchStrategy(computing_model_part,
            time_scheme,
            convergence_criterion,
            builder_and_solver,
            linesearch_settings)



















    def printDofsAndCPs(self) :
        # fig, ax = plt.subplots()

        free_node_x = []
        free_node_y = []
        free_node_z = []
        fixed_node_x = []
        fixed_node_y = []
        fixed_node_z = []
        dof = []

        z_ref = 1.0
        # Set Free the active ones
        if os.path.exists("txt_files/Id_active_control_points.txt"):
            with open('txt_files/Id_active_control_points.txt', 'r') as file:
                lines = file.readlines()
            for line in lines:
                numbers = line.split()
                node = self.main_model_part.GetNode(int(numbers[0]))
                node.Free(KratosMultiphysics.TEMPERATURE)
                node.Set(KratosMultiphysics.VISITED, False)


                if (np.abs(node.Y - z_ref) > 1e-1): continue 
                free_node_x.append(node.X)
                free_node_y.append(node.Y)
                free_node_z.append(node.Z)
                dof.append(numbers[1])
        
        dof2 = []
        free_node_x2 = []
        free_node_y2= []
        free_node_z2= []
        # Set Free the active ones
        if os.path.exists("txt_files/Id_active_control_points_condition.txt"):
            with open('txt_files/Id_active_control_points_condition.txt', 'r') as file:
                lines = file.readlines()
            for line in lines:
                numbers = line.split()
                node = self.main_model_part.GetNode(int(numbers[0]))

                if (np.abs(node.Y - z_ref) > 1e-1): continue 
                free_node_x2.append(node.X)
                free_node_y2.append(node.Y)
                free_node_z2.append(node.Z)
                dof2.append(numbers[1])
        
        for node in self.main_model_part.GetNodes() :

            if (np.abs(node.Y - z_ref) > 1e-1): continue 
            fixed_node_x.append(node.X)
            fixed_node_y.append(node.Y)
            fixed_node_z.append(node.Z)

        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # ax.scatter(fixed_node_x, fixed_node_y, fixed_node_z, marker='x', color='black', s=12,  label='Fixed Nodes')
        ax.scatter(free_node_x, free_node_y, free_node_z, marker='d', color='red',s=8, label='Free Nodes')

        ax.scatter(free_node_x2, free_node_y2, free_node_z2, marker='x', alpha=0.3, color='green',s=12, label='Free Nodes')
        
        # count = 0
        # for i in range(len(dof)) :
        #     ax.annotate(str(dof[i]), (free_node_x[i], free_node_y[i], free_node_z[i]), textcoords="offset points", xytext=(0,10), ha='center')

        # for i in range(len(dof2)) :
        #     ax.annotate(str(dof2[i]), (free_node_x2[i], free_node_y2[i], free_node_z2[i]), textcoords="offset points", color='green', xytext=(0,8), ha='center')
     
        # Read the refinements.iga.json
        # with open('refinements.iga.json', 'r') as file:
        #     refinements_parameters = json.load(file)
        # insert_nb_per_span_u = refinements_parameters['refinements'][0]['parameters']['insert_nb_per_span_u']
        # insert_nb_per_span_v = refinements_parameters['refinements'][0]['parameters']['insert_nb_per_span_v']

        # knots_u = [0.0]
        # knots_v = [0.0]
        # initial = 0.0
        # total = 2.0
        # for j in range(1, insert_nb_per_span_u + 1):
        #     knots_u.append(initial + total / (insert_nb_per_span_u + 1) * j)
        # for j in range(1, insert_nb_per_span_v + 1):
        #     knots_v.append(initial + total / (insert_nb_per_span_v + 1) * j)
        # knots_u.append(2.0)
        # knots_v.append(2.0)

        # ax.set_xlim(min(knots_u)-0.2, max(knots_u)+0.2)
        # ax.set_ylim(min(knots_v)-0.2, max(knots_v)+0.2)

        # for u in knots_u:
        #     ax.axvline(x=u, color='red', linestyle='--', linewidth=1.0)
        # for v in knots_v:
        #     ax.axhline(y=v, color='red', linestyle='--', linewidth=1.0)
        # ax.grid(True, linestyle='--', linewidth=0.5)
        # ax.gca().set_aspect('equal', adjustable='box')

        plt.show()


