# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

from statistics import mean, stdev
import math 

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")

    return CFLOutputProcess(model, settings["Parameters"])


class GenerateVelocityFluctuation(KratosMultiphysics.Process):
    """
    A class responsible for the velocity fluctuation output, which is an element value in Kratos.
    """

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"      : "",
                "constants":
                {
                "total_wavenumber_discretization" : 100
                }
                "time_step" : 0.1,
                "ABL_friction_velocity" : 0.656
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

        self.total_wavenumber_discretization = params["total_wavenumber_discretization"].GetInt()
        self.time_step = params["time_step"].GetDouble()
        self.ABL_friction_velocity = params["ABL_friction_velocity"].GetDouble()
    
    def InitializeModelConstants(self):
        # reading constants
        constants = self.model_settings["constants"]
        self.fluid_model_part.ProcessInfo[KratosRANS.TOTAL_WAVE_NUMBER] = constants["total_wavenumber_discretization"].GetInt()

        Kratos.Logger.PrintInfo(
            self.__class__.__name__,
            "The kinematic simulation strategy is created.")

    def AddVariables(self):
        # adding kinematic simulation specific variables
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_U)
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_V)
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_W)
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosRANS.EFFECTIVE_WAVE_NUMBER)
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_U)
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_V)
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_W)

        super(GenerateVelocityFluctuation, self).AddVariables()

    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        current_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

        if((current_time >= self.interval[0]) and (current_time < self.interval[1])) and (current_step % self.output_step == 0):
            cfl_value = self._EvaluateCFL()

            if (self.model_part.GetCommunicator().MyPID() == 0):
                output = self._SummarizeCFL(cfl_value)
                output_vals = [format(val, self.format) for val in output]

                # not formatting time in order to not lead to problems with time recognition
                # in the file writer when restarting
                output_vals.insert(0, str(current_time))

                res_labels = ["time: ", "mean: ", "std: ", "max: ", "cfl" +
                              "{:.1f}".format(self.cfl_output_limit) + ": ", "cfl1.0: "]

                if (self.print_to_screen):

                    result_msg = 'CFL evaluation for model part ' + \
                        self.model_part_name + '\n'
                    result_msg += ', '.join([a+b for a,
                                             b in zip(res_labels, output_vals)])
                    self._PrintToScreen(result_msg)

    
    def ExecuteFinalize(self):
        if (self.model_part.GetCommunicator().MyPID() == 0):
            self.output_file.close()

    def _DiscretiseWaveNumber(self):

        kinematic_viscosity = self.model_part.GetValue(Kratos.KINEMATIC_VISCOSITY)
        
        K = []
        Epsiron = []
        k = []
        for node in self.model_part.Nodes:
            K_node_i = node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY)
            epsiron_node_i = node.GetSolutionStepValue(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
            K.append(K_node_i)
            epsiron.append(epsiron_node_i)
            k_1_i = 2 * math.pi * epsiron_node_i / pow(K_node_i, 1.5)
            k_N_i = pow(epsiron_node_i, 0.25) / pow(kinematic_viscosity, 0.25)
            dk = (math.log10(k_N_i) - math.log10(k_1_i)) / total_wavenumber_discretization
            k_node_i = []
            k_i = k_1_i
            for i in range(total_wavenumber_discretization):
                k_node_i.append(k_i)
                k_i = math.exp(math.log10(k_1_i) + (i+1)*dk)
            k.append(k_node_i) 
        
        return k

    def _EstimateConstantsOfEnergySpectrum(self): 
    # assumption: z is height from ground, ground z = 0
    # !! need to be solved using c++ elements !! 
    # here Turbulent kinetic energy u,v,w is calculated to each node and will be used in elements!

        h = ABL_friction_velocity * 10000 / 6
        beta_v = []
        beta_w = []
        for node in self.model_part.Nodes:
            z = node.Z
            beta_v_node_i = 1-0.22*pow(math.cos(math.pi*z/(2*h)), 4)
            beta_w_node_i = 1-0.45*pow(math.cos(math.pi*z/(2*h)), 4) 
            beta_v.append(beta_v_node_i)
            beta_w.append(beta_w_node_i)
            K_node_i = node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY)
            Ku_node_i = (1/(1+beta_v_node_i**2+beta_w_node_i**2)) * K_node_i
            Kv_node_i = (beta_v_node_i**2/(1+beta_v_node_i**2+beta_w_node_i**2)) * K_node_i
            Kw_node_i = (beta_w_node_i**2/(1+beta_v_node_i**2+beta_w_node_i**2)) * K_node_i
            node.SetSolutionStepValue(KratosRNAS.TURBULENT_KINETIC_ENERGY_U, Ku_node_i)
            node.SetSolutionStepValue(KratosRNAS.TURBULENT_KINETIC_ENERGY_V, Kv_node_i)
            node.SetSolutionStepValue(KratosRNAS.TURBULENT_KINETIC_ENERGY_W, Kw_node_i)


    def _GetFileHeader(self):
        header = '# CFL for model part ' + self.model_part_name + \
            '| CFL_threshold: ' + str(self.cfl_output_limit) + '\n'
        header += '# Time Mean Std Max HowMany>' + \
            "{:.1f}".format(self.cfl_output_limit) + ' [%] HowMany>1.0 [%]\n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo(
            "CFLOutputProcess", "CFL VALUE RESULTS:")
        KratosMultiphysics.Logger.PrintInfo(
            "CFLOutputProcess", "Current time: " + result_msg)

    def _CalculateWithRespectToThreshold(self, x):

        y = [val for val in x if val < self.cfl_output_limit]
        y1 = [val for val in x if val < 1.0]
        # % of element with cfl above threshold
        how_many = ((len(x)-len(y))/len(x))*100
        # % of element with cfl above 1
        how_many1 = ((len(x)-len(y1))/len(x))*100

        # quantifying the mean and std for values below the threshold
        y_mean = mean(y)
        y_std = stdev(y)

        # qunatifying the global max
        x_max = max(x)
        return [y_mean, y_std, x_max, how_many, how_many1]

    def _EvaluateCFL(self):

        if (self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            KratosCFD.EstimateDtUtility2D.CalculateLocalCFL(self.model_part)
        else:
            KratosCFD.EstimateDtUtility3D.CalculateLocalCFL(self.model_part)

        local_cfl = []
        for elem in self.model_part.Elements:
            local_cfl.append(elem.GetValue(KratosMultiphysics.CFL_NUMBER))

        local_cfl = self.model_part.GetCommunicator().GetDataCommunicator().GathervDoubles(local_cfl, 0)

        return local_cfl

    def _SummarizeCFL(self, local_cfl):

        global_cfl = []
        for k in local_cfl:
            global_cfl.extend(k)

        cfl_mean, cfl_std, cfl_max, cfl_how_many, cfl_how_many1 = self._CalculateWithRespectToThreshold(
            global_cfl)

        return [cfl_mean, cfl_std, cfl_max, cfl_how_many, cfl_how_many1]
