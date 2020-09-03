# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.RANSApplication import RansVariableUtilities

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
                "ABL_friction_velocity" : 0.375,
                "seed_for_random_samples_generation": 2020,
                "lambda_unsteadiness_parameter" : 1.0,
                "constants":{
                    "total_wave_number": 10
                    }
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

        #self.total_wave_number = self.model_part.ProcessInfo[KratosRANS.TOTAL_WAVE_NUMBER]
        self.ABL_friction_velocity = params["ABL_friction_velocity"].GetDouble()
        self.rand_seed = params["seed_for_random_samples_generation"].GetInt()
        self.lambda_unsteadiness = params["lambda_unsteadiness_parameter"].GetDouble()

    #def InitializeModelConstants(self):
        # reading constants
        #constants = self.model_settings["constants"]
        constants = params["constants"]
        self.model_part.ProcessInfo[KratosRANS.TOTAL_WAVE_NUMBER] = constants["total_wave_number"].GetInt()
        self.total_wave_number = self.model_part.ProcessInfo[KratosRANS.TOTAL_WAVE_NUMBER]

        boundary_layer_height = self.ABL_friction_velocity*10000/6
        for node in self.model_part.Nodes:
            beta_v = 1 - 0.22 * pow(math.cos(math.pi*node.Z/(2*boundary_layer_height)),4)
            beta_w = 1 - 0.45 * pow(math.cos(math.pi*node.Z/(2*boundary_layer_height)),4)
            denom = 1+(beta_v**2)+(beta_w**2)
            K = node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY, 0)
            node.SetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_U, 0, K/denom)
            node.SetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_V, 0, K*(beta_v**2)/denom)
            node.SetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_W, 0, K*(beta_w**2)/denom)

        KratosMultiphysics.Logger.PrintInfo(
            self.__class__.__name__,
            "The kinematic simulation strategy is created.")

    def Check(self):
        # check whether used variables exist for debug
        self.model_part.HasNodalSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_U)
        self.model_part.HasNodalSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_V)
        self.model_part.HasNodalSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_W)
        self.model_part.HasNodalSolutionStepVariable(KratosRANS.EFFECTIVE_WAVE_NUMBER)
        self.model_part.HasNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.model_part.HasNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        self.model_part.HasNodalSolutionStepVariable(KratosMultiphysics.KINEMATIC_VISCOSITY)
        self.model_part.HasNodalSolutionStepVariable(KratosMultiphysics.LAGRANGE_DISPLACEMENT_X)

    def _CalculateFourierVariables(self):
        k_n = {}
        energy_spectrum = {}
        a_n = {}
        b_n = {}
        omega_n = {}
        end_time = 0.0
        i = 0
        for node in self.model_part.Nodes:
            k_n.update({node.Id:self._DiscretiseWaveNumber(node)}) # k_n: dict[node]*list[total_wavenumber_discretization]
            energy_spectrum.update({node.Id:self._CalculateEnergySpectrum(node, k_n[node.Id])}) # energy_spectrum: dict[node]*list[total_wavenumber_discretization][3]
            a_n.update({node.Id:self._GenerateFouerierCoefficient(energy_spectrum[node.Id], self.rand_seed)}) # a_n: dict[node]*list[total_wavenumber_discretization][3]
            b_n.update({node.Id:self._GenerateFouerierCoefficient(energy_spectrum[node.Id], self.rand_seed+1)}) # b_n: dict[node]*list[total_wavenumber_discretization][3]
            # a_n and b_n have same PDF but chage the seed so they have different values
            omega_n.update({node.Id:self._CalculateAngularFrequency(energy_spectrum[node.Id], k_n[node.Id])}) # omega_n: dict[node]*list[total_wavenumber]
            end_time_pre = 2*math.pi/min(omega_n[node.Id])
            if end_time_pre > end_time:
                end_time = end_time_pre
            i += 1

        return [k_n, energy_spectrum, a_n, b_n, omega_n, end_time]

    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        #current_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        if (not RansVariableUtilities.IsAnalysisStepCompleted(self.model_part, "FOURIER_SERIES_VARIABLES_CALCULATION")):
            theta_n = self._GenerateTheta() # theta_n: list[total_wavenumber_discretization]
            phi_n = self._GeneratePhi() # phi_n: list[total_wavenumber_discretization]
            self.k_n_unitvector = self._CalculateWaveNumberUnitvector(theta_n, phi_n) # k_n_unitvector: list[total_wavenumber_discretization][3]

            [self.k_n, self.energy_spectrum, self.a_n, self.b_n, self.omega_n, self.end_time] = self._CalculateFourierVariables()
            RansVariableUtilities.AddAnalysisStep(self.model_part, "FOURIER_SERIES_VARIABLES_CALCULATION")
        
        print('++++++++++ calculated end_time +++++++++++++++')
        print(self.end_time)
        if current_time >= self.end_time:
            print('aaaaaaaaa')
            RansVariableUtilities.AddAnalysisStep(self.model_part, "PERIOD_EXCEEDS_ENDTIME")

        for node in self.model_part.Nodes:
            u_fluc = self._FouerierSummation(node, self.k_n[node.Id], self.k_n_unitvector, self.a_n[node.Id], self.b_n[node.Id], self.omega_n[node.Id], current_time)
            # U should be the last step if results of URANS are used!
            u_x = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0) + u_fluc[0]
            u_y = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0) + u_fluc[1]
            u_z = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, 0) + u_fluc[2]
            node.SetSolutionStepValue(KratosMultiphysics.LAGRANGE_DISPLACEMENT_X, 0, u_x)
            node.SetSolutionStepValue(KratosMultiphysics.LAGRANGE_DISPLACEMENT_Y, 0, u_y)
            node.SetSolutionStepValue(KratosMultiphysics.LAGRANGE_DISPLACEMENT_Z, 0, u_z)

        print('~generate_velocity_fluctuation.py ExecuteFinalizeSolutionStep')
    
    '''
    def ExecuteFinalize(self):
        if (self.model_part.GetCommunicator().MyPID() == 0):
            self.output_file.close()
    '''

    def _GenerateTheta(self):

        np.random.seed(seed=self.rand_seed)
        samples = np.random.uniform(0,1,self.total_wave_number)
        samples = np.arccos(1-2*samples)

        return samples.tolist()

    def _GeneratePhi(self):

        np.random.seed(seed=self.rand_seed)
        samples = np.random.uniform(0,1,self.total_wave_number)
        samples = 2*math.pi*samples

        return samples.tolist()

    def _CalculateWaveNumberUnitvector(self, theta_n, phi_n):

        kn = list()
        for i in range(self.total_wave_number):
            k_x = math.sin(theta_n[i]) * math.cos(phi_n[i])
            k_y = math.sin(theta_n[i]) * math.sin(phi_n[i])
            k_z = math.cos(theta_n[i])
            kn.append([k_x, k_y, k_z])

        return kn

    def _DiscretiseWaveNumber(self, node):

        #kinematic_viscosity = self.model_part.GetValue(KratosMultiphysics.KINEMATIC_VISCOSITY)
        kinematic_viscosity = node.GetSolutionStepValue(KratosMultiphysics.KINEMATIC_VISCOSITY)

        K_node_i = node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY)
        epsiron_node_i = node.GetSolutionStepValue(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        k_1_i = 2 * math.pi * epsiron_node_i / pow(K_node_i, 1.5)
        k_N_i = pow(epsiron_node_i, 0.25) / pow(kinematic_viscosity, 0.75)
        dk = (math.log(k_N_i) - math.log(k_1_i)) / (self.total_wave_number-1)
        k_node_i = []
        k_i = k_1_i
        for i in range(self.total_wave_number):
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
        for i in range(self.total_wave_number):
            k_i = k_n[i]
            Eu_i = Au * 2 * Ku * math.exp(-2*((k_i/k_kol)**2)) / (ke*pow(1+(k_i/ke)**2, 0.833333333333333))
            Ev_i = Av * 2 * Kv * ((k_i/ke)**2) * math.exp(-2*((k_i/k_kol)**2)) / (ke*pow(1+(k_i/ke)**2, 1.833333333333333))
            Ew_i = Aw * 2 * Kw * ((k_i/ke)**2) * math.exp(-2*((k_i/k_kol)**2)) / (ke*pow(1+(k_i/ke)**2, 1.833333333333333))
            E.append([Eu_i, Ev_i, Ew_i])

        return E

    def _GenerateFouerierCoefficient(self, energy_spectrum, seed):

        np.random.seed(seed=seed) # where should the seed be...?
        samples = list()
        for i in range(self.total_wave_number):
            mean_vec = [0.0, 0.0, 0.0]
            cov_matrix = [[math.sqrt(2*energy_spectrum[i][0]), 0.0, 0.0],
                         [0.0, math.sqrt(2*energy_spectrum[i][1]), 0.0],
                         [0.0, 0.0, math.sqrt(2*energy_spectrum[i][2])]]

            samples_i = np.random.multivariate_normal(mean_vec, cov_matrix) #samples_i[1, 3]
            samples.append(samples_i.tolist())

        return samples

    def _CalculateAngularFrequency(self, energy_spectrum, k_n):

        omega_n = list()
        for i in range(self.total_wave_number):
            omega_n.append(self.lambda_unsteadiness * math.sqrt(pow(k_n[i],3)*(energy_spectrum[i][0]+energy_spectrum[i][1]+energy_spectrum[i][2])))

        return omega_n

    def _FouerierSummation(self, node, k_n, k_n_unitvector, a_n, b_n, omega_n, current_time):

        u_x = 0.0
        u_y = 0.0
        u_z = 0.0
        for i in range(self.total_wave_number):
            A = self.__cross(a_n[i], k_n_unitvector[i])
            B = self.__cross(b_n[i], k_n_unitvector[i])
            C = (k_n_unitvector[i][0]*node.X + k_n_unitvector[i][1]*node.Y + k_n_unitvector[i][2]*node.Z)*k_n[i]
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


