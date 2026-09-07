import KratosMultiphysics as Kratos
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.StructuralMechanicsApplication.step_controller import StepController, DefaultStepController, Factory

class StructuralMechanicsLoadSteppingAnalysis(StructuralMechanicsAnalysis):
    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """

        computing_mp: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()

        def __get_serializer():
            return Kratos.FileSerializer("serialization", Kratos.SerializerTraceType.SERIALIZER_NO_TRACE, True)

        self.is_converged = False
        while self.KeepAdvancingSolutionLoop():
            time_begin = self.time # current step lower bound
            # current step upper bound.
            #   - This also clones a new time step
            #   - This sets the DELTA_TIME (time_end - time_begin)
            time_end = self._AdvanceTime()

            # create the appropriate step controller
            if self.project_parameters.Has("step_controller_settings"):
                step_controller = Factory(time_begin, time_end, self.project_parameters["step_controller_settings"])
            else:
                step_controller = DefaultStepController(time_begin, time_end, Kratos.Parameters("""{}"""))

            self.time = time_end

            __get_serializer().Save(computing_mp.FullName(), computing_mp)

            # first try to solve for the final time.
            self.InitializeSolutionStep()
            self.is_converged = self._GetSolver().SolveSolutionStep()

            if not self.is_converged:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Did not converge for time = {self.time}.")

            current_step_controller_time = time_begin
            while not step_controller.IsCompleted(current_step_controller_time, self.is_converged):
                if self.is_converged:
                    Kratos.Logger.PrintInfo(self.__class__.__name__, f"Step at time = [{time_begin}, {self.time}] converged.")
                    # so the sub-step converged.

                    # first finalize the success full step
                    self.FinalizeSolutionStep()
                    self.OutputSolutionStep()

                    # make the next step ready for solving
                    computing_mp.ProcessInfo[Kratos.STEP] += 1

                    # here we will correctly set TIME and DELTA_TIME
                    # and put current values in the previous time step
                    time_begin = self.time
                    self.time = step_controller.GetNextStep(self.time, self.is_converged)
                    computing_mp.CloneTimeStep(self.time)

                    # now collect the nodes positions, which may have moved
                    # can be used at a later time to reset the nodal positions.
                    __get_serializer().Save(computing_mp.FullName(), computing_mp)
                else:
                    Kratos.Logger.PrintInfo(self.__class__.__name__, f"Step at time = {self.time} did not converge.")
                    # here we need to reset the coordinates of the mesh, if someone has used the move_mesh_flag = true
                    # reset the mesh coordinates
                    __get_serializer().Load(computing_mp.FullName(), computing_mp)

                    # sub_step did not converge
                    # do not advance in step. get a new sub-step
                    self.time = step_controller.GetNextStep(time_begin, self.is_converged)
                    computing_mp.ProcessInfo[Kratos.TIME] = self.time
                    computing_mp.ProcessInfo[Kratos.DELTA_TIME] = self.time - time_begin

                self.InitializeSolutionStep()
                self.is_converged = self._GetSolver().SolveSolutionStep()

                current_step_controller_time = self.time

            # we need the last finalize solution step, because once the Solver solves, and converges, and
            # [t_begin, t_end] is reached, it will no longer go in to the stepping while loop.
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = StructuralMechanicsLoadSteppingAnalysis(model, parameters)
    simulation.Run()
