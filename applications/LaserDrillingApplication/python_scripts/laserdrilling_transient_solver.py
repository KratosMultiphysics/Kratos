import os
import numpy as np
from scipy.interpolate import interp1d
import h5py
import time as timer
from sympy import *

import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
import KratosMultiphysics.LaserDrillingApplication as LaserDrillingApplication
if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_transient_solver

def CreateSolver(model, custom_settings):
    return LaserDrillingTransientSolver(model, custom_settings)

class SurfaceFEMProjector:
    def __init__(self, n_elements, R_far, H, density):
        self.n = n_elements # number of elements
        n = self.n
        self.R_far = R_far
        self.lumped_projection_option = True
        self.X = np.linspace(0.0, R_far, n+1)
        self.h = self.R_far / n
        self.A = np.zeros((n+1,n+1))
        self.b = np.zeros((n+1,1))        
        self.evaporated_volume_lists = [[] for i in range(n)]
        self.H = H
        self.density = density

    def NAsc(self, i, x):
        X = self.X
        return (x - X[i - 1]) / (X[i] - X[i - 1])

    def NDesc(self, i, x):
        X = self.X
        return (X[i + 1] - x) / (X[i + 1] - X[i])

    def N(self, i, x):
        NDesc = self.NDesc
        NAsc = self.NAsc        
        n = self.n
        X = self.X

        if i == 0:
            if np.any(x >= X[i]) and np.any(x <= X[i + 1]):
                return NDesc(i, x)
            else:
                return 0.0
        if i == n:
            if np.any(x >= X[i - 1]) and np.any(x <= X[i]):
                return NAsc(i, x)
            else:
                return 0.0
        if np.any(x >= X[i - 1]) and np.any(x <= X[i]):
            return NAsc(i, x)
        elif np.any(x >= X[i]) and np.any(x <= X[i + 1]):
            return NDesc(i,x)
        else:
            return 0.0

    def EvaluateFEMFunction(self, nodal_values, x):
        result = 0.0
        for i, U in enumerate(nodal_values):
            result += U * self.N(i, x)
        return result        

    def FillUpMassMatrix(self):
        X = self.X
        n = self.n
        A = self.A
        NDesc = self.NDesc
        NAsc = self.NAsc

        x = symbols('x')
        for i in range(0, n+1):
            for j in range(0, n+1):
                if i == j:
                    if i == 0:
                        integrand = NDesc(i, x) * NDesc(j, x) * x
                        A[i, j] = integrate(integrand, (x, X[0], X[1]))
                    elif i == n:
                        integrand = NAsc(i, x) * NAsc(j, x) * x
                        A[i, j] = integrate(integrand, (x, X[n-1], X[n]))
                    else:
                        integrand1 = NAsc(i, x)**2 * x
                        integrand2 = NDesc(i, x)**2 * x
                        A[i, j] = integrate(integrand1, (x, X[i-1], X[i])) + integrate(integrand2, (x, X[i], X[i+1]))
                elif j - i == 1:
                    integrand = NDesc(i, x) * NAsc(j, x) * x
                    A[i, j] = integrate(integrand, (x, X[i], X[j]))
        for i in range(0, n+1):
            for j in range(0, n+1):
                if i - j == 1:
                    A[i, j] = A[j, i]

        if self.lumped_projection_option:
            for i in range(0, n+1):
                sum_i = 0.0
                for j in range(0, n+1):
                    sum_i += A[i, j]
                A[i, i] = sum_i
                for j in range(0, n+1):
                    if i!=j:
                        A[i, j] = 0.0

    def FillUpDeltasRHS(self, evap_element_centers, evap_volumes):
        n = self.n
        N = self.N        
        b = self.b

        for i in range(n + 1):
            for r, vol in zip(evap_element_centers, evap_volumes):
                b[i] += self.H * self.density * vol * N(i, r)

    def FillUpFEMRHS(self, j):
        X = self.X
        n = self.n
        NDesc = self.NDesc
        NAsc = self.NAsc
        b = self.b

        x = symbols('x')
        for i in range(0, n+1):
            if i == j and i == 0:
                integrand = NDesc(i, x) * NDesc(j, x) #* x
                b[i] = integrate(integrand, (x, X[0], X[1]))
            elif i == j and i == n:
                integrand = NAsc(i, x) * NAsc(j, x) #* x
                b[i] = integrate(integrand, (x, X[n-1], X[n]))
            elif i == j:
                integrand1 = NAsc(i, x)**2 #* x
                integrand2 = NDesc(i, x)**2 #* x
                b[i] = integrate(integrand1, (x, X[i-1], X[i])) + integrate(integrand2, (x, X[i], X[i+1]))
            elif j - i == 1:
                integrand = NDesc(i, x) * NAsc(j, x) #* x
                b[i] = integrate(integrand, (x, X[i], X[j]))
            elif i - j == 1:
                integrand = NDesc(j, x) * NAsc(i, x) #* x
                b[i] = integrate(integrand, (x, X[j], X[i]))                

    def CalculateEnergyOfFEMFunction(self, nodal_values):
        X = self.X
        n = self.n
        NDesc = self.NDesc
        NAsc = self.NAsc

        total_energy = 0.0

        for i, u_value in enumerate(nodal_values):
            x = symbols('x')
            U = u_value[0]

            if i == 0:
                integrand = U * NDesc(i, x) * x
                total_energy += integrate(integrand, (x, X[i], X[i+1]))
            elif i == n:
                integrand = U * NAsc(i, x) * x
                total_energy += integrate(integrand, (x, X[i-1], X[i]))
            else:
                integrand1 = U * NAsc(i, x) * x
                integrand2 = U * NDesc(i, x) * x
                total_energy += integrate(integrand1, (x, X[i-1], X[i])) + integrate(integrand2, (x, X[i], X[i+1]))
        return float(total_energy)

    def InterpolateFunctionAndNormalize(self, f): #, normalization_value):
        X = self.X
        n = self.n

        nodal_values = np.zeros((n+1,1))
        for i, x in enumerate(X):
            nodal_values[[i]] = f(x)

        #obtained_energy = self.CalculateEnergyOfFEMFunction(nodal_values)
        #nodal_values *= normalization_value / obtained_energy
        return nodal_values       

    def AssignDeltasToTestFunctionSupports(self, radii):    
        X = self.X
        n = self.n        
        for i_rad, r in enumerate(radii):
            for j, list_of_centers in enumerate(self.evaporated_volume_lists):
                if j == 0 and r <= X[1]:
                    list_of_centers.append(i_rad)
                    continue
                elif j == n and r > X[n-1]:
                    list_of_centers.append(i_rad)
                    continue
                elif r > X[j-1] and r <= X[j+1]:
                    list_of_centers.append(i_rad)        

    def Project(self):
        self.solution = np.linalg.solve(self.A, self.b)
        return self.solution

    def q(self, r):
        return 80 * np.exp(-r)

class LaserDrillingTransientSolver(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver):
    def __init__(self, model, custom_settings):
        # Construct the base solver and validate the settings in base class
        super().__init__(model, custom_settings)   

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        self.jump_between_pulses_counter += self.delta_time
        if self.jump_between_pulses_counter >= self.time_jump_between_pulses:
            self.jump_between_pulses_counter = 0
            self.pulse_number += 1
            self.RemoveElementsByAblation()
            self.ResidualHeatStage()

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters(r"""{
            "time_integration_method" : "implicit",
            "transient_parameters" : {
                "dynamic_tau": 1.0,
                "theta"    : 0.5
            },
            "ambient_temperature" : 0.0
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AllocateKratosMemory(self):
        # Set element counter variable to zero
        for elem in self.main_model_part.Elements:
            elem.SetValue(LaserDrillingApplication.THERMAL_COUNTER, 0)
            elem.SetValue(KratosMultiphysics.TEMPERATURE, 0.0)
            elem.SetValue(LaserDrillingApplication.THERMAL_DECOMPOSITION, 0.01)
            elem.SetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, 0.0)
            elem.Set(KratosMultiphysics.ACTIVE, True)

        for node in self.main_model_part.Nodes:
            node.SetValue(LaserDrillingApplication.DECOMPOSED_NODE, 0.0)

    def SetParameters(self):
        self.element_id_to_study = 1578

        self.some_elements_are_above_the_evap_temp = True
        self.plot_decomposed_volume_graph = False
        self.jump_between_pulses_counter = 0
        self.number_of_decomposed_elements = 0
        self.pulse_number = 1

        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()

        with open(materials_filename, 'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())

        self.material_settings = materials["properties"][0]["Material"]

        with open("ProjectParameters.json", 'r') as project_parameters_file:
            project_parameters = KratosMultiphysics.Parameters(project_parameters_file.read())

        if not project_parameters["problem_data"].Has("energy"):
            self.Q = 25e-6
        else:
            self.Q = project_parameters["problem_data"]["energy"].GetDouble()
        if not project_parameters["problem_data"].Has("residual_heat_fraction"):
            self.residual_heat_fraction = 0.05
        else:
            self.residual_heat_fraction = project_parameters["problem_data"]["residual_heat_fraction"].GetDouble()
        if not project_parameters["problem_data"].Has("mask_aperture_diameter"):
            mask_aperture_diameter = 0.025
        else:
            mask_aperture_diameter = project_parameters["problem_data"]["mask_aperture_diameter"].GetDouble()
        if not project_parameters["problem_data"].Has("vaporisation_temperature"):
            self.T_e = 1000.0
        else:
            self.T_e = project_parameters["problem_data"]["vaporisation_temperature"].GetDouble()
        if not project_parameters["problem_data"].Has("time_jump_between_pulses"):
            self.time_jump_between_pulses = 1e6
        else:
            self.time_jump_between_pulses = project_parameters["problem_data"]["time_jump_between_pulses"].GetDouble()
        if not project_parameters["problem_data"].Has("compute_vaporisation"):
            self.compute_vaporisation = False
        else:
            self.compute_vaporisation = project_parameters["problem_data"]["compute_vaporisation"].GetBool()
        if not self.material_settings["Variables"].Has("ENTHALPY"):
            self.evaporation_enthalpy = 4e5 # J/Kg. Value found on the internet for a given epoxy resin.
        else:
            self.evaporation_enthalpy = self.material_settings['Variables']['ENTHALPY'].GetDouble()
        if self.compute_vaporisation:
            self.decomposed_nodes_coords_filename = "list_of_decomposed_nodes_coords_with_evap.txt"
        else:
            self.decomposed_nodes_coords_filename = "list_of_decomposed_nodes_coords_no_evap.txt"

        self.R_far = mask_aperture_diameter * 0.5

        self.cp = self.material_settings['Variables']['SPECIFIC_HEAT'].GetDouble()
        self.conductivity = self.material_settings['Variables']['CONDUCTIVITY'].GetDouble()
        self.rho = self.material_settings['Variables']['DENSITY'].GetDouble()
        self.T0 = self.settings['ambient_temperature'].GetDouble()
        self.kappa = self.conductivity / (self.rho * self.cp)
        self.ablation_energy_fraction = 1 - self.residual_heat_fraction
        self.sigma = 0.5 * self.R_far
        self.K = 1 / (2 * self.sigma**2)
        self.C = self.ablation_energy_fraction * self.Q * self.K / (np.pi * (1 - np.exp(-self.K * self.R_far**2)))
        self.irradiated_surface_area = np.pi * self.R_far**2

        # Problem calibrated for Q = 15e-6 (3W power). For different powers, the ablated volume should change linearly as experiments show
        Q_reference = 15e-6
        self.V = 4.72e-7 # mm3. Approximate ablated volume for 1 pulses (experimental). For 5 pulses it should be around 2.36e-6
        self.V *= self.Q / Q_reference

        if self.ablation_energy_fraction:
            # Find F_th_fraction multiplying F_th so Radius_th = R_far
            # This gives a F_th_fraction of 0.313, a little too flat (maximum not captured)
            # F_th_fraction = self.C * np.exp(-self.K * self.R_far**2) / (self.ablation_energy_fraction * self.Q / A)

            ####################################
            # Calibrated for:                  #
            # residual heat fraction   = 0.05  #
            # Vaporisation temperature = 1000K #
            ####################################
            if not self.compute_vaporisation:
                F_th_fraction = 0.333 # 5 pulses, no evaporation
            else:
                F_th_fraction = 0.333 #0.43 # 5 pulses, with evaporation

            self.F_th = F_th_fraction * self.ablation_energy_fraction * self.Q / self.irradiated_surface_area
            self.radius_th = np.sqrt(np.log(self.C / self.F_th) / self.K)
            
            if not self.compute_vaporisation:
                self.l_s = self.PenetrationDepth() # It gives self.l_s = 0.002mm # 5 pulses, no evaporation
            else:
                self.l_s = self.PenetrationDepth() #0.001 # 5 pulses, with evaporation     

        self.l_th = 0.05 # micrometers # Thermal depth, assumed

        l_th_in_meters = self.l_th * 1e-6
        kappa_in_square_meters = self.kappa * 1e-6
        self.thermal_penetration_time = l_th_in_meters**2 / kappa_in_square_meters

        # Finite Elements
        self.n_surface_elements = 10 # number of elements
        self.lumped_projection_option = True
        self.surface_nodes_Y_values = np.linspace(0.0, self.R_far, self.n_surface_elements + 1)
        self.surface_element_size = self.R_far / self.n_surface_elements

    def SetUpResultsFiles(self):
        self.SetUpGNUPlotFiles()
        self.SetUpHDF5Files()

    def SetUpGNUPlotFiles(self):
        if os.path.exists("temperature_alpha.txt"):
            os.remove("temperature_alpha.txt")
        self.temperature_alpha_file = open("temperature_alpha.txt", "a")
        if os.path.exists("time_alpha.txt"):
            os.remove("time_alpha.txt")
        self.time_alpha_file = open("time_alpha.txt", "a")
        if os.path.exists("decomposed_volume_evolution.txt"):
            os.remove("decomposed_volume_evolution.txt")
        self.decomposed_volume_file = open("decomposed_volume_evolution.txt", "a")

    def SetUpHDF5Files(self):
        radius_2 = lambda node: node.X**2 + node.Y**2 + node.Z**2
        self.near_field_nodes = [node for node in self.main_model_part.Nodes if radius_2(node) < self.R_far**2]
        self.radii = np.sqrt(np.array([radius_2(node) for node in self.near_field_nodes]))
        self.results_filename = 'results.h5'
        self.CreateResultsFile(self.results_filename)
        self.temperature_increments = np.array([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes])        

    def Initialize(self):
        super(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver, self).Initialize()
        self.starting_time = timer.time()
        self.SetParameters()

        self.AllocateKratosMemory()

        self.list_of_decomposed_nodes_coords = []
        self.list_of_decomposed_nodes_coords_X = []
        self.list_of_decomposed_nodes_coords_Y = []

        self.SetUpResultsFiles()

        # TODO: Initial condition, ambient temperature. Do this using GUI!
        for node in self.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, self.T0)
            #node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, self.T0)
        self.IdentifyInitialSurfaceNodes()

        initial_system_energy = self.MonitorEnergy()
        
        self.RemoveElementsByAblation()

        computed_energy_after_ablation = self.MonitorEnergy()

        residual_heat = self.residual_heat_fraction * self.Q

        self.ResidualHeatStage()

        system_energy_after_residual_heat = self.MonitorEnergy()
        print("\nPulse_number:", self.pulse_number)
        print("Initial system energy:", initial_system_energy)
        print("\nComputed energy after laser:", computed_energy_after_ablation)
        print("\nResidual_heat:", residual_heat)
        print("System energy after residual heat:", system_energy_after_residual_heat)
        difference_after_and_before_residual_heat = system_energy_after_residual_heat - computed_energy_after_ablation
        print("Difference after and before residual heat:", difference_after_and_before_residual_heat, "\n")

        self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

    def IdentifyInitialSurfaceNodes(self):
        self.list_of_decomposed_nodes_ids = []
        self.list_of_decomposed_elements_ids = []            
        self.list_of_decomposed_nodes_coords_X = []
        self.list_of_decomposed_nodes_coords_Y = []

        for elem in self.main_model_part.Elements:
            for node in elem.GetNodes():
                # TODO: extremely ad-hoc!
                if node.X < 0.0000001 and node.Y <= self.R_far:
                    self.list_of_decomposed_nodes_ids.append(node.Id)
                    self.list_of_decomposed_elements_ids.append(elem.Id)    
                    self.list_of_decomposed_nodes_coords_X.append(node.X)
                    self.list_of_decomposed_nodes_coords_Y.append(node.Y)                

        self.list_of_decomposed_nodes_ids = list(set(self.list_of_decomposed_nodes_ids))
        self.list_of_decomposed_elements_ids = list(set(self.list_of_decomposed_elements_ids))

        if not self.main_model_part.HasSubModelPart("BoundaryPart"):
            self.main_model_part.CreateSubModelPart('BoundaryPart')
            self.boundary_part = self.main_model_part.GetSubModelPart('BoundaryPart')

        self.boundary_part.AddElements(self.list_of_decomposed_elements_ids)
        self.boundary_part.AddNodes(self.list_of_decomposed_nodes_ids)

    def ResidualHeatStage(self):
        if self.residual_heat_fraction:
            self.projector = SurfaceFEMProjector(self.n_surface_elements, self.R_far, self.evaporation_enthalpy, self.rho) #, delta_coefficients)
            self.q_interp = self.projector.InterpolateFunctionAndNormalize(self.EnergyPerUnitArea1D) #, 1.0)
            if self.compute_vaporisation:
                #self.original_Q_value = self.Q
                self.total_removed_energy = 0.0
                self.first_evaporation_stage_done = False
                max_iterations = 20
                iterations = 0
                self.some_elements_are_above_the_evap_temp = True
                print("Removing elements by evaporation...")
                print("Pulse number:", self.pulse_number)

                while self.some_elements_are_above_the_evap_temp:
                    #self.ResetDomainTemperature()
                    self.ImposeTemperatureIncreaseDueTo1DConduction()
                    self.RemoveElementsByEvaporation()
                    iterations += 1
                    print("Evaporated layer:", iterations)
                    if iterations > max_iterations:
                        print("******************************MAXIMUM ITERATIONS EXCEEDED!!!")
                        break
                print("Done!")
                #self.Q = self.original_Q_value
                self.number_of_decomposed_elements = 1
            else:
                self.ImposeTemperatureIncreaseDueTo1DConduction()

    def MarkElementsToEvaporate(self):
        pass

    def RemoveElementsByEvaporation(self):
        self.number_of_decomposed_elements = 0
        number_of_boundary_elements = 0
        evap_elements_centers_Y = []
        evap_elements_volumes  = []
        self.some_elements_are_above_the_evap_temp = False
        for elem in self.boundary_part.Elements:
            if elem.Is(KratosMultiphysics.ACTIVE):
                temp = elem.CalculateOnIntegrationPoints(KratosMultiphysics.TEMPERATURE, self.main_model_part.ProcessInfo)
                element_temperature = temp[0]
                if element_temperature > self.T_e:
                    self.some_elements_are_above_the_evap_temp = True
                    elem.Set(KratosMultiphysics.ACTIVE, False)
                    vol = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, self.main_model_part.ProcessInfo)
                    elem.SetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, vol[0])                    
                    element_volume = elem.GetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME)
                    self.total_removed_energy += self.rho * element_volume * self.evaporation_enthalpy
                    self.number_of_decomposed_elements += 1
                    Y_centroid = elem.GetGeometry().Center().Y
                    evap_elements_centers_Y.append(Y_centroid)
                    evap_elements_volumes.append(element_volume)
                    for node in elem.GetNodes():
                        node.SetValue(LaserDrillingApplication.DECOMPOSED_NODE, 1.0)

        number_of_problematic_elements = 0
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.ACTIVE):
                number_of_decomposed_nodes = 0
                for node in elem.GetNodes():
                    if node.GetValue(LaserDrillingApplication.DECOMPOSED_NODE):
                        number_of_decomposed_nodes +=1
                if number_of_decomposed_nodes == 3:
                    elem.Set(KratosMultiphysics.ACTIVE, False)
                    Y_centroid = elem.GetGeometry().Center().Y
                    evap_elements_centers_Y.append(Y_centroid)
                    evap_elements_volumes.append(element_volume)                    
                    elem.CalculateOnIntegrationPoints(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, self.main_model_part.ProcessInfo)
                    element_volume = elem.GetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME)
                    self.total_removed_energy += self.rho * element_volume * self.evaporation_enthalpy
                    # self.number_of_decomposed_elements += 1
                    number_of_problematic_elements += 1
        evap_elements_centers_Y = np.array(evap_elements_centers_Y)
        evap_elements_volumes  = np.array(evap_elements_volumes)
        self.evap_elements_centers_Y = evap_elements_centers_Y[evap_elements_centers_Y.argsort()]
        self.evap_elements_volumes  = evap_elements_volumes[evap_elements_centers_Y.argsort()]

        self.projector.FillUpMassMatrix()
        self.projector.FillUpDeltasRHS(self.evap_elements_centers_Y, self.evap_elements_volumes)
        self.projector.AssignDeltasToTestFunctionSupports(self.evap_elements_volumes)
        self.u = self.projector.Project()
        self.q_interp -= self.u

        for i, q in enumerate(self.q_interp):
            if q < 0.0:
                self.q_interp[i] = 0.0

        import matplotlib.pyplot as plt

        plt.grid()
        plt.plot(self.projector.X, self.u, color='blue', marker='o')
        #plt.plot(projector.X, q_cont, color='black', marker='x')
        plt.plot(self.projector.X, self.q_interp, color='red', marker='+')
        #plt.plot(self.projector.X, remaining_energy, color='green', marker='*')
        #plt.plot(self.evap_elements_centers_Y, q_eval, color='black', marker='x')
        plt.legend(["E/A lost", "q FEM proj", "remaining energy FEM", "remaining energy FEM eval"], loc="center left")
        #plt.show()

        self.AddDecomposedNodesToSurfaceList()
        self.first_evaporation_stage_done = True

    def ResetDomainTemperature(self):
        if self.first_evaporation_stage_done:
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, self.T0)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, self.T0)
            for elem in self.main_model_part.Elements:
                elem.SetValue(KratosMultiphysics.TEMPERATURE, self.T0)

    def ImposeTemperatureIncreaseDueTo1DConduction(self):
        X = self.list_of_decomposed_nodes_coords_X
        Y = self.list_of_decomposed_nodes_coords_Y
        self.minimum_characteristic_Z =  1e6
        self.maximum_characteristic_Z = -1e6
        for node in self.main_model_part.Nodes:
            radius = node.Y
            """ if not self.ablation_energy_fraction:
                distance_to_surface = 0.0
            else: """
            F = interp1d(Y, X, bounds_error=False, fill_value=0.0)
            distance_to_surface = F(radius)
            z = node.X - distance_to_surface
            if radius <= self.R_far:
                delta_temp = self.TemperatureVariationInZDueToLaser1D(radius, z)
                old_temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                new_temp = old_temp + delta_temp
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temp)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, new_temp)
        #print("\nResidual heat fraction:", self.residual_heat_fraction)
        #print("\nMaximum characteristic depth:", self.maximum_characteristic_Z)
        problem_characteristic_time_minimum_depth = self.minimum_characteristic_Z**2 / self.kappa
        problem_characteristic_time_maximum_depth = self.maximum_characteristic_Z**2 / self.kappa
        minimum_time_step_for_minimum_depth = 0.1 * problem_characteristic_time_minimum_depth
        minimum_time_step_for_maximum_depth = 0.1 * problem_characteristic_time_maximum_depth
        #print("\nThermal problem characteristic time for maximum depth:", problem_characteristic_time_maximum_depth)
        #print("\nNecessary time step for maximum depth:", minimum_time_step_for_maximum_depth, '\n')

    def EnergyPerUnitArea1D(self, radius):
        C = (1 - self.ablation_energy_fraction) * self.Q * self.K / (np.pi * (1 - np.exp(-self.K * self.R_far**2)))
        q =  C * np.exp(-self.K * radius**2)
        return q

    def InitialThermalConductionTime(self, radius):
        # This function returns the characteristic time required for the initial heat distribution from the surface to the interior
        # 4.0 (2**2) in numerator due to equation (4c) in Weber, 2014. 'Heat accumulation during pulsed laser materials processing'
        # C = 4.0 / (4.0 * np.pi * self.rho**2 * self.kappa * self.cp**2 * (self.T_e - self.T0)**2)
        # t = C * self.EnergyPerUnitArea1D(radius)**2
        return self.thermal_penetration_time

    def TemperatureVariationInZDueToLaser1D(self, radius, z):
        #q = self.EnergyPerUnitArea1D(radius)
        q = self.projector.EvaluateFEMFunction(self.q_interp, radius)
        t_penetration = self.InitialThermalConductionTime(radius)
        # 2.0 in numerator due to equation (4c) in Weber, 2014. 'Heat accumulation during pulsed laser materials processing'
        C = 2.0 * q / (self.rho * self.cp * 2.0 * np.sqrt(np.pi * self.kappa * t_penetration))
        delta_temp = C * np.exp(- z**2 / (4.0 * self.kappa * t_penetration))
        characteristic_Z = np.sqrt(4.0 * self.kappa * t_penetration)
        if characteristic_Z <= self.minimum_characteristic_Z:
            self.minimum_characteristic_Z = characteristic_Z
        if characteristic_Z >= self.maximum_characteristic_Z:
            self.maximum_characteristic_Z = characteristic_Z
        return delta_temp
    
    def FEMInterpolatedTemperatureVariationInZDueToLaser1D(self, radius, z):
        pass

    def SolveSolutionStep(self):
        super(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver, self).SolveSolutionStep()
        self.temperature_increments = np.array([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes])
        self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        decomp_vol = self.MonitorDecomposedVolume()
        decomp_vol *= 1e9 # To convert mm3 into um3
        current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.decomposed_volume_file.write(str(current_time) + " " + str(decomp_vol) + "\n")
        for elem in self.main_model_part.Elements:
            if elem.Id == self.element_id_to_study:
                temperature = elem.GetValue(KratosMultiphysics.TEMPERATURE)
                thermal_decomposition = elem.GetValue(LaserDrillingApplication.THERMAL_DECOMPOSITION)
                self.temperature_alpha_file.write(str(temperature) + " " + str(thermal_decomposition) + "\n")
                self.time_alpha_file.write(str(current_time) + " " + str(thermal_decomposition) + "\n")
                break

    def Finalize(self):
        super().Finalize()
        self.decomposed_volume_file.close()
        self.temperature_alpha_file.close()
        self.time_alpha_file.close()
        if self.plot_decomposed_volume_graph:
            self.PrintDecomposedVolumeEvolution()
        elapsed_time = timer.time() - self.starting_time
        print("\nElapsed_time:", elapsed_time, '\n')

    def PrintDecomposedVolumeEvolution(self):
        import matplotlib.pyplot as plt
        file = open(os.path.expanduser("decomposed_volume_evolution.txt"))
        lines = file.readlines()
        x, y = [], []
        for line in lines:
            x.append(line.split()[0])
            y.append(line.split()[1])
        file.close()
        plt.plot(x,y)
        plt.show()

    def MonitorDecomposedVolume(self):
        decomposed_volume = 0.0
        for elem in self.main_model_part.Elements:
            out = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, self.main_model_part.ProcessInfo)
            decomposed_volume += out[0]
        return decomposed_volume

    def MonitorEnergy(self):
        energy = 0.0
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.ACTIVE):
                # NOTE: Here out is an std::vector with all components containing the same elemental thermal energy
                out = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.THERMAL_ENERGY, self.main_model_part.ProcessInfo)
                energy += out[0]
        return energy

    def RemoveElementsByAblation(self):
        if self.ablation_energy_fraction:
            X = self.list_of_decomposed_nodes_coords_X
            Y = self.list_of_decomposed_nodes_coords_Y

            for elem in self.main_model_part.Elements:
                X_centroid = elem.GetGeometry().Center().X
                Y_centroid = elem.GetGeometry().Center().Y
                # if self.pulse_number == 1:
                #     X_interp = 0
                # else:
                F = interp1d(Y, X, bounds_error=False)
                X_interp = F(Y_centroid)
                DeltaX = X_centroid - X_interp
                d_ev = self.EvaporationDepth(Y_centroid)
                if DeltaX <= d_ev: # and Y_centroid <= self.radius_th:
                    elem.Set(KratosMultiphysics.ACTIVE, False)
                    for node in elem.GetNodes():
                        node.SetValue(LaserDrillingApplication.DECOMPOSED_NODE, 1.0)
            for elem in self.main_model_part.Elements:
                if elem.Is(KratosMultiphysics.ACTIVE):
                    number_of_decomposed_nodes = 0
                    for node in elem.GetNodes():
                        if node.GetValue(LaserDrillingApplication.DECOMPOSED_NODE):
                            number_of_decomposed_nodes += 1
                    if number_of_decomposed_nodes == 3:
                        elem.Set(KratosMultiphysics.ACTIVE, False)
            self.AddDecomposedNodesToSurfaceList()
            print('\nR_far:', self.R_far)
            print('\nRadius_th:', self.radius_th)
            print("\nDecomposed volume:", self.MonitorDecomposedVolume())

    def sortSecond(self, val):
        return val[1]

    def AddDecomposedNodesToSurfaceList(self):
        self.list_of_decomposed_nodes_ids = []
        self.list_of_decomposed_elements_ids = []
        self.list_of_decomposed_nodes_coords = []
        self.list_of_decomposed_nodes_coords_X = []
        self.list_of_decomposed_nodes_coords_Y = []
        number_of_boundary_elements = 0
        for elem in self.main_model_part.Elements:
            first_decomposed_node_found = False
            if elem.Is(KratosMultiphysics.ACTIVE):
                for node in elem.GetNodes():
                    if node.GetValue(LaserDrillingApplication.DECOMPOSED_NODE):
                        self.list_of_decomposed_nodes_ids.append(node.Id)
                        if not first_decomposed_node_found:
                            number_of_boundary_elements += 1
                            self.list_of_decomposed_elements_ids.append(elem.Id)
                            first_decomposed_node_found = True
        self.list_of_decomposed_nodes_ids = list(set(self.list_of_decomposed_nodes_ids))

        if not self.main_model_part.HasSubModelPart("BoundaryPart"):
            self.main_model_part.CreateSubModelPart('BoundaryPart')
            self.boundary_part = self.main_model_part.GetSubModelPart('BoundaryPart')
        else:
            self.main_model_part.RemoveSubModelPart(self.boundary_part)
            self.main_model_part.CreateSubModelPart('BoundaryPart')
            self.boundary_part = self.main_model_part.GetSubModelPart('BoundaryPart')
        self.boundary_part.AddElements(self.list_of_decomposed_elements_ids)
        self.boundary_part.AddNodes(self.list_of_decomposed_nodes_ids)

        for node in self.main_model_part.Nodes:
            if node.Id in self.list_of_decomposed_nodes_ids:
                X = node.X
                Y = node.Y
                coords = [X, Y]
                self.list_of_decomposed_nodes_coords.append(coords)
        self.list_of_decomposed_nodes_coords.sort(key=self.sortSecond)
        self.list_of_decomposed_nodes_coords_X = [coord[0] for coord in self.list_of_decomposed_nodes_coords]
        self.list_of_decomposed_nodes_coords_Y = [coord[1] for coord in self.list_of_decomposed_nodes_coords]
        if os.path.exists(self.decomposed_nodes_coords_filename):
            os.remove(self.decomposed_nodes_coords_filename)
        self.decomposed_nodes_coords_file = open(self.decomposed_nodes_coords_filename, "a")
        for coord in self.list_of_decomposed_nodes_coords:
            self.decomposed_nodes_coords_file.write(str(coord[0]) + " " + str(coord[1]) + "\n")
        self.decomposed_nodes_coords_file.close()

    def PenetrationDepth(self):
        F_th = self.F_th
        V = self.V
        R_th = self.radius_th
        l_s = V / (np.pi * (0.5 * R_th**2 * (np.log(self.C / F_th)) - 0.25 * self.K * R_th**4))
        return l_s

    def EvaporationDepth(self, r):
        if r >= self.radius_th:
            return 0
        else:
            q = self.C * np.exp(-self.K * r**2)
            d_ev = 0.5 * self.l_s * (np.log(q) - np.log(self.F_th))
            return d_ev

    def CreateResultsFile(self, filename):
        if os.path.exists(self.results_filename):
            os.remove(self.results_filename)
        with h5py.File(filename, 'a') as f:
            f.attrs['ambient_temperature'] = self.T0
            f.attrs['pulse_energy'] = self.Q
            f.attrs['specific_heat_capacity'] = self.cp
            f.attrs['density'] = self.rho
            f.attrs['conductivity'] = self.conductivity
            # Create a dataset to store the radii
            dataset = f.create_dataset('radii', (self.radii.shape), dtype=self.radii.dtype)
            dataset[:] = self.radii[:]
            f.create_group('temperature_increments')

    def WriteResults(self, filename, process_info):
        step = process_info[KratosMultiphysics.STEP]
        time = step = process_info[KratosMultiphysics.TIME]
        # Open the HDF5 file.
        with h5py.File(filename, 'a') as f:
            assert self.radii.shape  == self.temperature_increments.shape
            # Create a dataset to store the radii and temperatures data.
            dataset = f['/temperature_increments'].create_dataset(str(step), self.temperature_increments.shape, dtype=self.temperature_increments.dtype)
            # Write the radii and temperatures data to the dataset.
            dataset[:] = self.temperature_increments
            # Add a time label to the dataset.
            dataset.attrs["time"] = time

    def _get_element_condition_replace_settings(self):
        ## Get and check domain size
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2,3]:
            raise Exception("DOMAIN_SIZE is not set in ProcessInfo container.")
        ## Validate the replace settings
        default_replace_settings = self.GetDefaultParameters()["element_replace_settings"]
        self.settings["element_replace_settings"].ValidateAndAssignDefaults(default_replace_settings)
        ## Elements
        ## Note that we check for the elements that require substitution to allow for custom elements
        element_name = self.settings["element_replace_settings"]["element_name"].GetString()
        element_list = ["EulerianConvDiff","LaplacianElement","MixedLaplacianElement","AdjointHeatDiffusionElement","QSConvectionDiffusionExplicit","DConvectionDiffusionExplicit","LaserAxisymmetricEulerianConvectionDiffusion"]
        if element_name in element_list:
            num_nodes_elements = 0
            if (len(self.main_model_part.Elements) > 0):
                for elem in self.main_model_part.Elements:
                    num_nodes_elements = len(elem.GetNodes())
                    break
            num_nodes_elements = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_elements)
            if not num_nodes_elements:
                num_nodes_elements = domain_size + 1
            name_string = f"{element_name}{domain_size}D{num_nodes_elements}N"
            self.settings["element_replace_settings"]["element_name"].SetString(name_string)
        ## Conditions
        condition_name = self.settings["element_replace_settings"]["condition_name"].GetString()
        condition_list = ["FluxCondition","ThermalFace","AxisymmetricThermalFace","LineCondition","SurfaceCondition"]
        if condition_name in condition_list:
            num_nodes_conditions = 0
            if (len(self.main_model_part.Conditions) > 0):
                for cond in self.main_model_part.Conditions:
                    num_nodes_conditions = len(cond.GetNodes())
                    break
            num_nodes_conditions = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_conditions)
            if not num_nodes_conditions:
                num_nodes_conditions = domain_size
            name_string = f"{condition_name}{domain_size}D{num_nodes_conditions}N"
            self.settings["element_replace_settings"]["condition_name"].SetString(name_string)
        return self.settings["element_replace_settings"]
