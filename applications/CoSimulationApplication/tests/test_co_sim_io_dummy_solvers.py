from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

import os, subprocess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestCoSimIODummySolvers(KratosUnittest.TestCase):

    def setUp(self):
        with open("co_sim_io_dummy_solver_test/cosim_parameters.json", 'r') as parameter_file:
            self.cosim_parameters = KM.Parameters(parameter_file.read())
        self.start_external_solver = True
        self.echo_level_ext_solver = 0 # can be added to each test if desired

    def test_weak_coupling_Cpp_solver(self):
        self.name_executable = "dummy_solver_cpp"
        self.execute_test_weak_coupling()

    def test_weak_coupling_Cpp_solver_run_external(self):
        self.name_executable = "dummy_solver_cpp"
        self.start_external_solver = False
        self.execute_test_weak_coupling()

    def test_strong_coupling_Cpp_solver(self):
        self.name_executable = "dummy_solver_cpp"
        self.execute_test_strong_coupling()

    def test_co_sim_orchestrated_weak_coupling_Cpp_solver(self):
        self.name_executable = "dummy_solver_cpp"
        self.execute_test_co_sim_orchestrated_weak_coupling()

    def test_co_sim_orchestrated_strong_coupling_Cpp_solver(self):
        self.name_executable = "dummy_solver_cpp"
        self.execute_test_co_sim_orchestrated_strong_coupling()


    # def test_weak_coupling_C_solver(self):
    #     self.skipTest("This test is not finished")
    #     self.name_executable = "dummy_solver_c"
    #     self.execute_test_weak_coupling()

    # def test_strong_coupling_C_solver(self):
    #     self.skipTest("This test is not finished")
    #     self.name_executable = "dummy_solver_c"
    #     self.execute_test_strong_coupling()

    # def test_co_sim_orchestrated_weak_coupling_C_solver(self):
    #     self.skipTest("This test is not finished")
    #     self.name_executable = "dummy_solver_c"
    #     self.execute_test_co_sim_orchestrated_weak_coupling()

    # def test_co_sim_orchestrated_strong_coupling_C_solver(self):
    #     self.skipTest("This test is not finished")
    #     self.name_executable = "dummy_solver_c"
    #     self.execute_test_co_sim_orchestrated_strong_coupling()


    def test_weak_coupling_Fortran_solver(self):
        self.name_executable = "dummy_solver_fortran"
        self.execute_test_weak_coupling()

    def test_strong_coupling_Fortran_solver(self):
        self.name_executable = "dummy_solver_fortran"
        self.execute_test_strong_coupling()

    # def test_co_sim_orchestrated_weak_coupling_Fortran_solver(self):
    #     self.skipTest("This test is not finished")
    #     self.name_executable = "dummy_solver_fortran"
    #     self.execute_test_co_sim_orchestrated_weak_coupling()

    # def test_co_sim_orchestrated_strong_coupling_Fortran_solver(self):
    #     self.skipTest("This test is not finished")
    #     self.name_executable = "dummy_solver_fortran"
    #     self.execute_test_co_sim_orchestrated_strong_coupling()


    def execute_test_weak_coupling(self):
        self.coupling_level = 1
        self.__ModifyExternalSolverCoSimSettings()
        self.execute_test()

    def execute_test_strong_coupling(self):
        self.coupling_level = 2
        self.__ModifyCoSimParametersForStrongCoupling()
        self.__ModifyExternalSolverCoSimSettings()
        self.execute_test()

    def execute_test_co_sim_orchestrated_weak_coupling(self):
        self.coupling_level = 3
        self.__ModifyExternalSolverCoSimSettings()
        self.execute_test()

    def execute_test_co_sim_orchestrated_strong_coupling(self):
        self.coupling_level = 3
        self.__ModifyCoSimParametersForStrongCoupling()
        self.__ModifyExternalSolverCoSimSettings()
        self.execute_test()


    def __ModifyExternalSolverCoSimSettings(self):
        self.external_solver_start_command = "../../../install/libs/" + self.name_executable # TODO find better solution, this can be different btw different installations

        ext_solver_params = self.cosim_parameters["solver_settings"]["solvers"]["external_dummy_solver"]["solver_wrapper_settings"]

        ext_solver_params.AddEmptyValue("start_external_solver").SetBool(self.start_external_solver)
        ext_solver_params.AddEmptyValue("controlling_external_solver").SetBool(self.coupling_level==3) # this means co-sim-orchestrated coupling
        if self.start_external_solver:
            ext_solver_params.AddEmptyValue("external_solver_start_command").SetString(self.external_solver_start_command)
            ext_solver_params.AddEmptyList("external_solver_arguments")
            ext_solver_params["external_solver_arguments"].Append(str(self.coupling_level))
            ext_solver_params["external_solver_arguments"].Append(str(self.echo_level_ext_solver))

    def __ModifyCoSimParametersForStrongCoupling(self):
        self.cosim_parameters["solver_settings"]["type"].SetString("coupled_solvers.gauss_seidel_strong")
        self.cosim_parameters["solver_settings"].AddEmptyValue("num_coupling_iterations").SetInt(3)

        self.cosim_parameters["solver_settings"].AddValue("convergence_criteria", KM.Parameters("""[
            {
                "type"          : "relative_norm_previous_residual",
                "solver"        : "helpers_dummy_kratos_solver",
                "data_name"     : "interface_pressure"
            }
        ]"""))

    def execute_test(self):
        if self.start_external_solver:
            CoSimulationAnalysis(self.cosim_parameters).Run()
        else:
            subprocess_cmd = [self.external_solver_start_command, str(self.coupling_level), str(self.echo_level_ext_solver)]
            p = subprocess.Popen(subprocess_cmd)
            CoSimulationAnalysis(self.cosim_parameters).Run()
            p.communicate()
            self.assertEqual(p.returncode, 0)

if __name__ == '__main__':
    KratosUnittest.main()
