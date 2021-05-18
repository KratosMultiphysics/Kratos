# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_criteria import CoSimulationConvergenceCriteria

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
from numpy import linalg as la

def Create(settings, solvers):
    cs_tools.SettingsTypeCheck(settings)
    return AbsoluteNormResidualConvergenceCriteria(settings, solvers)

class AbsoluteNormResidualConvergenceCriteria(CoSimulationConvergenceCriteria):
    def __init__(self, settings, solvers):
        super().__init__(settings)

        if self.settings["criteria_composition"].GetString() != "energy_conjugate":
            self.__RaiseException('Energy conjugate criteria composition requires energy conjugate variables to be specified in "data_name" and "conjugate_data_name".')

        self.iteration = 1
        self.abs_tolerance = self.settings["abs_tolerance"].GetDouble()
        self.ignore_first_convergence = self.settings["ignore_first_convergence"].GetBool()

        # Determine if we are looking at the energy difference between two domains (solvers), or just one
        self.solver_vec = [solvers[settings["solver"].GetString()]]
        is_dual_domain = False
        for criteria_option in settings["criteria_options"]:
            if criteria_option.GetString() == "domain_difference":
                is_dual_domain = True
                break
        if is_dual_domain:
            is_error = False
            if settings.Has("solver_domain_two"):
                solver_domain_two = settings["solver_domain_two"].GetString()
                if solver_domain_two == "UNSPECIFIED":
                    is_error = True
                else:
                    self.solver_vec.append(solvers[solver_domain_two])
            else:
                is_error = True
            if is_error:
                self.__RaiseException('Domain difference requires "solver_domain_two" to be set to the second domain.')

        # Setup interface data matrix (general form)
        self.interface_data = [None]*len(solvers)
        for solver_index in range(0,len(self.interface_data)):
            self.interface_data[solver_index] = [self.solver_vec[solver_index].GetInterfaceData(settings["data_name"].GetString())]
            self.interface_data[solver_index].append(self.solver_vec[solver_index].GetInterfaceData(settings["conjugate_data_name"].GetString()))

        self.second_domain_data_sign = 1.0
        if "swap_second_domain_data_sign" in settings["criteria_options"].GetStringArray():
            self.second_domain_data_sign = -1.0

        settings.RemoveValue("data_name")
        settings.RemoveValue("solver")

        if not settings.Has("label"):
            settings.AddEmptyValue("label").SetString(colors.bold('{}.{}'.format(self.interface_data[0][0].solver_name, self.interface_data[0][0].name)))
        self.label = self.settings["label"].GetString()

    def IsConverged(self):
        # Compute energy scalar on interface
        current_data = 0.0

        for solver_index in range(0,len(self.interface_data)):
            #check length of data vectors are the same
            interface_energy = 0.0;
            data_1 = self.interface_data[solver_index][0].GetData()
            data_2 = self.interface_data[solver_index][1].GetData()
            if len(data_1) != len(data_2):
                self.__RaiseException('Data vector lengths for conjugate criteria composition must be identical, but they are different!')
            else:
                for i in range(0,len(data_1)):
                    interface_energy += data_1[i]*data_2[i]
            if solver_index == 0:
                current_data = interface_energy
            else:
                current_data -= self.second_domain_data_sign*interface_energy #assumes domain_difference

        abs_norm = la.norm(current_data)

        if self.ignore_first_convergence and self.iteration == 1:
            is_converged = False
        else:
            is_converged = abs_norm < self.abs_tolerance

        self.iteration += 1

        info_msg = ""

        if self.echo_level > 1:
            info_msg  = 'Convergence '

            if self.label != "":
                info_msg += 'for "{}": '.format(self.label)

            if is_converged:
                info_msg += colors.green("ACHIEVED")
            else:
                info_msg += colors.red("NOT ACHIEVED")

        if self.echo_level > 2:
            info_msg += '\n\t abs-norm = {:.2e} | abs-tol = {}'.format(abs_norm, self.abs_tolerance)

        if info_msg != "":
            cs_tools.cs_print_info(self._ClassName(), info_msg)

        return is_converged

    def Check(self):
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "abs_tolerance"            : 1e-5,
            "criteria_options"         : [],
            "conjugate_data_name"      : "UNSPECIFIED",
            "solver_domain_two"        : "UNSPECIFIED",
            "criteria_composition"     : "UNSPECIFIED",
            "use_wrapper"              : false,
            "ignore_first_convergence" : false,
            "label"                    : ""
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

