# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

from statistics import mean, stdev
import math
import numpy as np

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")

    return GenerateVelocityFluctuationProcess(model, settings["Parameters"])


class GenerateVelocityFluctuationProcess(KratosMultiphysics.Process):
    """
    A class responsible for the velocity fluctuation output, which is an element value in Kratos.
    """

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"      : "",
                "total_time_steps" : 100,
                "ABL_friction_velocity" : 0.375,
                "seed_for_random_samples_generation": 2020,
                "lamda_unsteadiness_parameter" : 1.0
            }
            """)
        #####################################################################################
        # ABL friction velocity: constant used for generation of logarithmic inlet velocity #
        #####################################################################################
        params.ValidateAndAssignDefaults(default_settings)

        # getting the ModelPart from the Model
        self.model_part_name = params["model_part_name"].GetString()
        if self.model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        else:
            self.model_part = model[self.model_part_name]

        self.total_time_steps = params["total_time_steps"].GetInt()
        self.ABL_friction_velocity = params["ABL_friction_velocity"].GetDouble()
        self.rand_seed = params["seed_for_random_samples_generation"].GetInt()
        self.lambda_unsteadiness = param["lamda_unsteadiness_parameter"].GetDouble()

    def InitializeModelConstants(self):
        # reading constants
        constants = self.model_settings["constants"]
        self.fluid_model_part.ProcessInfo[KratosRANS.TOTAL_WAVE_NUMBER] = constants["total_wavenumber_discretization"].GetInt()

        Kratos.Logger.PrintInfo(
            self.__class__.__name__,
            "The kinematic simulation strategy is created.")

    def Check(self):
        # check whether used variables exist for debug
        self.fluid_model_part.HasSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_U)
        self.fluid_model_part.HasSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_V)
        self.fluid_model_part.HasSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_W)
        self.fluid_model_part.HasSolutionStepVariable(KratosRANS.EFFECTIVE_WAVE_NUMBER)
        self.fluid_model_part.HasSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.fluid_model_part.HasSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        self.fluid_model_part.HasSolutionStepVariable(KratosRANS.KINEMATIC_VISCOSITY)

    def _CalculateFourierVariables(self):
        k_n = list()
        energy_spectrum = list()
        a_n = list()
        b_n = list()
        omega_n = list()
        end_time = 0.0
        for node in self.model_part.Nodes:
            k_n.append(self._DiscretiseWaveNumber(node)) # k_n: list[node][total_wavenumber_discretization]
            energy_spectrum.append(self._CalculateEnergySpectrum(node, k_n)) # energy_spectrum: list[node][total_wavenumber_discretization][3]
            a_n.append(self._GenerateFouerierCoefficient(energy_spectrum, self.rand_seed)) # a_n: list[node][total_wavenumber_discretization][3]
            b_n.append(self._GenerateFouerierCoefficient(energy_spectrum, self.rand_seed+1)) # b_n: list[node][total_wavenumber_discretization][3]
            # a_n and b_n have same PDF but chage the seed so they have different values
            omega_n.append(self._CalculateAngularFrequency(energy_spectrum, k_n)) # omega_n: list[node][total_wavenumber]
            end_time_pre = 2*math.pi/min(omega_n[i])
            if end_time_pre > end_time:
                end_time = end_time_pre
        
        return [k_n, energy_spectrum, a_n, b_n, omega_n, end_time]

    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        #current_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

        theta_n = self._GenerateTheta() # theta_n: list[total_wavenumber_discretization]
        phi_n = self._GeneratePhi() # phi_n: list[total_wavenumber_discretization]
        k_n_unitvector = self._CalculateWaveNumberUnitvector(theta_n, phi_n) # k_n_unitvector: list[total_wavenumber_discretization][3]

        [k_n, energy_spectrum, a_n, b_n, omega_n, end_time] = self._CalculateFourierVariables()
        
        while current_time < end_time:
            i = 0
            for node in self.model_part.Nodes:
                u_fluc = self._FouerierSummation(node, k_n[i], k_n_unitvector, a_n[i], b_n[i], omega_n[i], current_time)
                # U should be the last step if results of URANS are used!
                u_x = node.GetSolutionStepValue(VELOCITY_X, 0) + u_fluc[0]
                u_y = node.GetSolutionStepValue(VELOCITY_Y, 0) + u_fluc[1]
                u_z = node.GetSolutionStepValue(VELOCITY_Z, 0) + u_fluc[2]
                node.SetSolutionStepValue(LAGRANGE_DISPLACEMENT_X, 0, u_x)
                node.SetSolutionStepValue(LAGRANGE_DISPLACEMENT_Y, 0, u_y)
                node.SetSolutionStepValue(LAGRANGE_DISPLACEMENT_Z, 0, u_z)
                i += 1
        
    def ExecuteFinalize(self):
        if (self.model_part.GetCommunicator().MyPID() == 0):
            self.output_file.close()

    def _GenerateTheta(self):

        np.random.seed(seed=self.rand_seed)
        samples = np.random.uniform(0,1,self.total_wavenumber_discretization)
        samples = np.arccos(1-2*samples)

        return samples.tolist()

    def _GeneratePhi(self):

        np.random.seed(seed=self.rand_seed)
        samples = np.random.uniform(0,1,self.total_wavenumber_discretization)
        samples = 2*math.pi()*samples

        return samples.tolist()

    def _CalculateWaveNumberUnitvector(self, theta_n, phi_n):

        kn = list()
        for i in range(self.total_wavenumber_discretization):
            k_x = math.sin(theta_n[i]) * math.cos(phi_n[i])
            k_y = math.sin(theta_n[i]) * math.sin(phi_n[i])
            k_z = math.cos(theta_n[i])
            kn.append([kx, ky, kz])

        return kn

    def _DiscretiseWaveNumber(self, node):

        kinematic_viscosity = self.model_part.GetValue(Kratos.KINEMATIC_VISCOSITY)

        K_node_i = node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY)
        epsiron_node_i = node.GetSolutionStepValue(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        k_1_i = 2 * math.pi * epsiron_node_i / pow(K_node_i, 1.5)
        k_N_i = pow(epsiron_node_i, 0.25) / pow(kinematic_viscosity, 0.75)
        dk = (math.log(k_N_i) - math.log(k_1_i)) / (self.total_wavenumber_discretization-1)
        k_node_i = []
        k_i = k_1_i
        for i in range(total_wavenumber_discretization):
            k_node_i.append(k_i)
            k_i = math.exp(math.log(k_1_i) + (i+1)*dk)

        return k_node_i

    def _CalculateEnergySpectrum(self, node, k_n):

        Au = node.GetSolutionStepValue(KratosRANS.SPECTRAL_CONSTANT_U)
        Av = node.GetSolutionStepValue(KratosRANS.SPECTRAL_CONSTANT_V)
        Aw = node.GetSolutionStepValue(KratosRANS.SPECTRAL_CONSTANT_W)
        ke = node.GetSolutionStepValue(KratosRANS.EFFECTIVE_WAVE_NUMBER)
        Ku = node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_U)
        Kv = node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_V)
        Kw = node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_W)
        k_kol = k_n[-1]
        E = list()

        for i in range(self.total_wavenumber_discretization):
            k_i = k_n[i]
            Eu_i = Au * 2 * Ku * math.exp(-2*((k_i/k_kol)**2)) / (ke*pow(1+(k_i/ke)**2, 0.833333333333333))
            Ev_i = Av * 2 * Kv * ((k_i/ke)**2) * math.exp(-2*((k_i/k_kol)**2)) / (ke*pow(1+(k_i/ke)**2, 1.833333333333333))
            Ew_i = Aw * 2 * Kw * ((k_i/ke)**2) * math.exp(-2*((k_i/k_kol)**2)) / (ke*pow(1+(k_i/ke)**2, 1.833333333333333))
            E.append([Eu_i, Ev_i, Ew_i])

        return E

    def _GenerateFouerierCoefficient(self, energy_spectrum, seed):

        np.random.seed(seed=seed) # where should the seed be...?
        samples = list()
        for i in range(self.total_wavenumber_discretization):
            mean_vec = [0.0, 0.0, 0.0]
            cov_matrix = [[sqrt(2*energy_spectrum[i][0]), 0.0, 0.0],
                         [0.0, sqrt(2*energy_spectrum[i][1]), 0.0],
                         [0.0, 0.0, sqrt(2*energy_spectrum[i][2])]]

            samples_i = np.random.multivariate_normal(mean_vec, cov_matrix, 1) #samples_i[1, 3]
            samples.append(samples_i.tolist())

        return samples

    def _CalculateAngularFrequency(self, energy_spectrum, k_n):

        omega_n = self.lamda_unsteadiness * np.sqrt(np.pow(k_n,3)*(energy_spectrum[0]+energy_spectrum[1]+energy_spectrum[2]))

        return omega_n.tolist()

    def _FouerierSummation(self, node, k_n, k_n_unitvector, a_n, b_n, omega_n, current_time):

        u_x = 0.0
        u_y = 0.0
        u_z = 0.0
        for i in range(self.total_wavenumber_discretization):
            A = self.__cross(a_n, k_n_unitvector)
            B = self.__cross(b_n, k_n_unitvector)
            C = (k_n_unitvector[0]*node.X + k_n_unitvector[1]*node.Y + k_n_unitvector[2]*node.Z)*k_n[i]
            u_x += A[0]*math.cos(C+omega_n[i]*current_time) + B[0]*math.sin(C+omega_n[i]*current_time)
            u_y += A[1]*math.cos(C+omega_n[i]*current_time) + B[1]*math.sin(C+omega_n[i]*current_time)
            u_z += A[2]*math.cos(C+omega_n[i]*current_time) + B[2]*math.sin(C+omega_n[i]*current_time)

        return [u_x, u_y, u_z]

    # cross product for list
    def __cross(self, a, b):

        c = list()
        c.append(a[1]*b[2]-a[2]*b[1])
        c.append(a[2]*b[0]-a[0]*b[2])
        c.append(a[0]*b[1]-a[1]*b[0])

        return c

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo(
            "CFLOutputProcess", "CFL VALUE RESULTS:")
        KratosMultiphysics.Logger.PrintInfo(
            "CFLOutputProcess", "Current time: " + result_msg)


