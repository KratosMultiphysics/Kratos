import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics import KratosUnittest

# python imports
import os
import sys
import subprocess

from pathlib import Path


def GetPython3Command():
    """Return the name of the python command, can be used with subprocess."""
    sys_executable = sys.executable
    if sys_executable: # was found
        return sys_executable
    raise Exception("The python command could not be determined!")


def ReadModelPart(mdpa_file_name, model_part, settings=None):
    """This method is designed to read a ModelPart.

    The ModelPart is filled with respective nodes, conditions, elements
    from mdpa file using either MPI or Serial reading depending on run type.

    Args:
        mdpa_file_name (str): Name of the mdpa file (without ".mdpa" extension)
        model_part (Kratos.ModelPart): ModelPart to be filled
    """
    if model_part.NumberOfNodes() > 0:
        raise Exception("ModelPart must not contain Nodes!")

    communicator = Kratos.Testing.GetDefaultDataCommunicator()
    if communicator.IsDistributed():
        ReadDistributedModelPart(mdpa_file_name, model_part, settings)
    else:
        ReadSerialModelPart(mdpa_file_name, model_part)


def ReadSerialModelPart(mdpa_file_name, model_part):
    """Reads mdpa file

    This method reads mdpa file and fills given model_part accordingly without MPI

    Args:
        mdpa_file_name (str): Name of the mdpa file (without ".mdpa" extension)
        model_part (Kratos.ModelPart): ModelPart to be filled
    """
    import_flags = Kratos.ModelPartIO.READ | Kratos.ModelPartIO.SKIP_TIMER
    Kratos.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)


@KratosUnittest.skipIfApplicationsNotAvailable("MetisApplication")
def ReadDistributedModelPart(mdpa_file_name, model_part, importer_settings=None):
    """Reads mdpa file

    This method reads mdpa file and fills given model_part accordingly using MPI

    Args:
        mdpa_file_name (str): Name of the mdpa file (without ".mdpa" extension)
        model_part (Kratos.ModelPart): ModelPart to be filled
    """

    from KratosMultiphysics.mpi import distributed_import_model_part_utility
    model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)

    if importer_settings is None:
        importer_settings = Kratos.Parameters("""{
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": \"""" + mdpa_file_name + """\",
                "partition_in_memory" : true
            },
            "echo_level" : 0
        }""")

    model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(
        model_part, importer_settings)
    model_part_import_util.ImportModelPart()
    model_part_import_util.CreateCommunicators()


def PrintTestHeader(application):
    Kratos.Logger.Flush()
    print(f"\nRunning {application} tests", file=sys.stderr, flush=True)


def PrintTestFooter(application, exit_code):
    Kratos.Logger.Flush()
    appendix = f" with exit code {exit_code}!" if exit_code != 0 else "."
    print(f"Completed {application} tests{appendix}\n", file=sys.stderr, flush=True)


def PrintTestSummary(exit_codes):
    Kratos.Logger.Flush()
    print("Test results summary:", file=sys.stderr, flush=True)
    max_test_name_length = len(max(exit_codes.keys(), key=len))+1
    for test, exit_code in exit_codes.items():
        result_string = "OK" if exit_code == 0 else "FAILED"
        pretty_name = test.ljust(max_test_name_length)
        print(f"  {pretty_name}: {result_string}", file=sys.stderr, flush=True)
    sys.stderr.flush()


class Commander(object):
    def __init__(self):
        self.process = None
        self.exitCodes = {}

    def PrintOutput(self, output, channel):
        ''' Prints the output of the process line by line.
            This prevents the github actions log from buffering too much and casuing output to be printed out of order.

            Ideallt this should detect if we are in the ci and only use the line/print in that case.
        '''

        for line in output.decode('utf8').split('\n'):
            print(line, file=channel)

    def TestToAppName(self, application):
        ''' Converts the name of a test suit into an application
        '''
        return application[6:-8] + "Application"

    def MPITestToAppName(self, application):
        ''' Converts the name of a test suit into an application
        '''
        return application.replace("MPI","")[6:-8] + "Application"
    
    def _RunTest(self, test_suit_name, command, timer, working_dir=os.getcwd()):
        # Print test header
        PrintTestHeader(test_suit_name)

        process_stdout = None
        process_stderr = None

        try:
            self.process = subprocess.Popen(command, cwd=working_dir, stdout=subprocess.PIPE)
        except OSError:
            print(f'[Error]: Unable to execute "{command}"', file=sys.stderr)
            self.exitCodes[test_suit_name] = 1
        except ValueError:
            # Command does exist, but the arguments are invalid (It sohuld never enter here. Just to be safe)
            print(f'[Error]: Invalid arguments when calling "{command}"', file=sys.stderr)
            self.exitCodes[test_suit_name] = 1
        else:
            # Used instead of wait to "soft-block" the process and prevent deadlocks
            # and capture the first exit code different from OK
            try:
                # Pipe the output of the process until the timer is reached
                process_stdout, process_stderr = self.process.communicate(timeout=timer)
                
                # Capture the result of the process
                self.exitCodes[test_suit_name] = int(self.process.returncode != 0)
            except subprocess.TimeoutExpired:
                # Timeout reached
                self.process.kill()
                print(f'[Error]: Tests for {test_suit_name} took too long. Process Killed.', file=sys.stderr)
                self.exitCodes[test_suit_name] = 1
            except Exception as e:
                # Unknown error
                print(f"[Error]: Unhandled exception while running {test_suit_name} tests: {e}", file=sys.stderr)
                self.exitCodes[test_suit_name] = 1
            finally:
                if process_stdout:
                    self.PrintOutput(process_stdout, sys.stderr)
                if process_stderr:
                    self.PrintOutput(process_stderr, sys.stderr)

        # Exit message
        PrintTestFooter(test_suit_name, self.process.returncode)
        
    def RunPythonTests(self, applications, level, verbose, command, timer):
        ''' Calls the script that will run the tests.

        Input
        -----
        application: list of strings
            Name of the ap.

        level: string
            minimum level of the test that will be run if possible.

        verbose: int
            detail of the ouptut. The grater the verbosity level, the greate the
            detail will be.

        command: string
            command to be used to call the tests. Ex: Python, Python3

        timer: integer
            limit time considered to execute the tests

        '''

        fullApplicationList = ["KratosCore"] + kratos_utils.GetListOfAvailableApplications()

        # If no applications are selected by the user, run all applications
        if not applications:
            applications = fullApplicationList

        # Iterate over the list of all applications an execute the ones selected by the user
        for application in fullApplicationList:
            if application in applications:
                    
                    if application == "KratosCore":
                        test_script = Path(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath())) / Path("kratos") / Path("tests") / Path("test_{}.py".format(application))
                    else:
                        test_script = Path(Kratos.KratosPaths.kratos_applications) / application / Path("tests") / Path("test_{}.py".format(application))

                    if os.path.isfile(test_script):
                        # Run all the tests in the executable
                        self._RunTest(
                            test_suit_name=application, 
                            command=filter(None, [
                                command, 
                                test_script, 
                                f"-v{verbose}", 
                                f"-l{level}"
                            ]),
                            timer=timer,
                            working_dir=os.path.dirname(os.path.abspath(test_script))
                        )
                    else:
                        if verbose > 0:
                            print(
                                '[Warning]: No test script found for {}'.format(
                                    application),
                                file=sys.stderr)
                            sys.stderr.flush()
                        if verbose > 1:
                            print(
                                '  expected file: "{}"'.format(
                                    test_script),
                                file=sys.stderr)
                            sys.stderr.flush()

    def RunCppTests(self, applications, timer, config):
        ''' Calls the cpp tests directly
        '''

        # Iterate over all executables that are not mpi dependant and execute them.
        print(kratos_utils.GetKratosMultiphysicsPath())
        for test_suite in os.listdir(os.path.join(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()), "test")):
            filename = str(Path(os.fsdecode(test_suite)).with_suffix('')) # Name of the file without extension 
            binfname = os.fsdecode(test_suite)                            # Name of the file with extension

            working_dir = os.getcwd()

            if filename in config and "working_dir" in config[filename]:
                working_dir = config[filename]["working_dir"]

            # Skip mpi tests
            if ("MPI" not in filename and self.TestToAppName(filename) in applications) or filename == "KratosCoreTest":
                
                # Run all the tests in the executable
                self._RunTest(
                    test_suit_name=filename, 
                    command=[
                        os.path.join(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()),"test",binfname)
                    ], 
                    timer=timer,
                    working_dir=working_dir
                )

    def RunMPIPythonTests(self, applications, mpi_command, mpi_flags, num_processes_flag, num_processes, level, verbose, command, timer):

        fullApplicationList = ["KratosMPICore"] + kratos_utils.GetListOfAvailableApplications()

        # mpi_flags may need to be passed using quotes by some executors. This removes the quotes if they are present.
        mpi_flags = mpi_flags.split(" ")

        # If no applications are selected by the user, run all applications
        if not applications:
            applications = fullApplicationList

        # Iterate over the list of all applications an execute the ones selected by the user
        for application in fullApplicationList:
            if application in applications:

                if application == "KratosMPICore":
                    test_script = Path(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath())) / Path("kratos") / Path("mpi") / Path("tests") / Path("test_{}.py".format(application))
                else:
                    test_script = Path(Kratos.KratosPaths.kratos_applications) / application / Path("tests") / Path("test_{}.py".format(application + "_mpi"))

                if Path.is_file(test_script):
                    # Run all the tests in the executable
                    self._RunTest(
                        test_suit_name=application, 
                        command=filter(None, [
                            mpi_command, 
                            *mpi_flags, 
                            num_processes_flag, 
                            str(num_processes), 
                            command, 
                            test_script,
                            "--using-mpi",
                            f"-v{verbose}", 
                            f"-l{level}"
                        ]),
                        timer=timer,
                        working_dir=os.path.dirname(os.path.abspath(str(test_script)))
                    )
                    
                else:
                    if verbose > 0:
                        print('[Warning]: No test script found for {}'.format(application), file=sys.stderr, flush=True)
                    if verbose > 1:
                        print('  expected file: "{}"'.format(test_script), file=sys.stderr, flush=True)

    def RunMPICppTests(self, applications, verbosity, mpi_command, mpi_flags, num_processes_flag, num_processes):
        ''' Calls the mpi cpp tests directly
        '''

        # mpi_flagss may need to be passed using quotes by some executors. This removes the quotes if they are present.
        mpi_flags = mpi_flags.split(" ")

        # Iterate over all executables that are mpi dependant and execute them.
        for test_suite in os.listdir(os.path.join(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()), "test")):
            filename = str(Path(os.fsdecode(test_suite)).with_suffix('')) # Name of the file without extension 
            binfname = os.fsdecode(test_suite)                            # Name of the file with extension
            
            # Skip non-mpi tests
            if ("MPI" in filename and self.MPITestToAppName(filename) in applications) or filename == "KratosMPICoreTest":

                # Run all the tests in the executable
                self._RunTest(
                    test_suit_name=filename, 
                    command=filter(None, [
                        mpi_command, 
                        *mpi_flags, 
                        num_processes_flag, 
                        str(num_processes), 
                        os.path.join(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()),"test",binfname)
                    ]),
                    timer=300)
