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

    def TestToAppName(self, application):
        ''' Converts the name of a test suit into an application
        '''
        return application[6:-8] + "Application"

    def MPITestToAppName(self, application):
        ''' Converts the name of a test suit into an application
        '''
        return application.replace("MPI","")[6:-8] + "Application"
    
    def _RunTest(self, test_suit_name, level, verbose, command, timer, workin_dir=os.getcwd()):
        # Print test header
        PrintTestHeader(test_suit_name)

        self.process = subprocess.Popen(
            command,
            cwd=workin_dir, stdout=subprocess.PIPE
        )

        # Used instead of wait to "soft-block" the process and prevent deadlocks
        # and capture the first exit code different from OK
        try:
            # timeout should not be a problem for cpp, but we leave it just in case
            timer = int(180)
            process_stdout, process_stderr = self.process.communicate(timeout=timer)
        except subprocess.TimeoutExpired:
            # Timeout reached
            self.process.kill()
            print('[Error]: Tests for {} took too long. Process Killed.'.format(test_suit_name), file=sys.stderr)
            self.exitCodes[test_suit_name] = 1
        else:
            if process_stdout:
                print(process_stdout.decode('utf8'), file=sys.stdout, flush=True)
            if process_stderr:
                print(process_stderr.decode('utf8'), file=sys.stderr, flush=True)

        # Exit message
        PrintTestFooter(test_suit_name, self.process.returncode)

        # Running out of time in the tests will send the error code -15. We may want to skip
        # that one in a future. Right now will throw everything different from 0.
        self.exitCodes[test_suit_name] = int(self.process.returncode != 0)
        
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

        fullApplicationList = kratos_utils.GetListOfAvailableApplications()

        # If no applications are selected by the user, run all applications
        if not applications:
            applications = fullApplicationList

        # Iterate over the list of all applications an execute the ones selected by the user
        for application in fullApplicationList:
            if application in applications:
                appNormalizedPath = application.lower().replace('_', '')

                possiblePaths = [
                    {'Found': p, 'FoundNormalized': p.split('/')[-1].lower().replace('_', ''), 'Expected': application, 'ExpectedNormalized': appNormalizedPath} for p in os.listdir(Kratos.KratosPaths.kratos_applications) if p.split('/')[-1].lower().replace('_', '') == appNormalizedPath
                ]

                if len(possiblePaths) < 1:
                    if verbose > 0:
                        print(
                            '[Warning]: No directory found for {}'.format(
                                application),
                            file=sys.stderr)
                        sys.stderr.flush()
                elif len(possiblePaths) > 1:
                    if verbose > 0:
                        print('Unable to determine correct path for {}'.format(application), file=sys.stderr)
                        print(
                            'Please try to follow the standard naming convention \'FooApplication\' Snake-Capital string  without symbols.',
                            file=sys.stderr)
                    if verbose > 1:
                        print('Several possible options were found:', file=sys.stderr)
                        for p in possiblePaths:
                            print('\t', p, file=sys.stderr)
                else:
                    script = Kratos.KratosPaths.kratos_applications+'/'+possiblePaths[0]['Found']+'/tests/'+'test_'+application+'.py'
                    print(script, file=sys.stderr)

                    if possiblePaths[0]['Found'] != possiblePaths[0]['Expected']:
                        print(
                            '[Warning]: Application has been found in "{}" directory but it was expected in "{}". Please check the naming convention.'.format(
                                possiblePaths[0]['Found'],
                                possiblePaths[0]['Expected']),
                            file=sys.stderr)

                    if os.path.isfile(script):
                        try:
                            self.process = subprocess.Popen([
                                command,
                                script,
                                '-l'+level,
                                '-v'+str(verbose)
                            ], stdout=subprocess.PIPE, cwd=os.path.dirname(os.path.abspath(script)))
                        except OSError:
                            # Command does not exist
                            print('[Error]: Unable to execute {}'.format(command), file=sys.stderr)
                            self.exitCodes[application] = 1
                        except ValueError:
                            # Command does exist, but the arguments are invalid (It sohuld never enter here. Just to be safe)
                            print('[Error]: Invalid arguments when calling {} {} {} {}'.format(command, script, '-l'+level, '-v'+str(verbose)), file=sys.stderr)
                            self.exitCodes[application] = 1
                        else:
                            # Used instead of wait to "soft-block" the process and prevent deadlocks
                            # and capture the first exit code different from OK
                            try:
                                process_stdout, process_stderr = self.process.communicate(timeout=timer)
                            except subprocess.TimeoutExpired:
                                # Timeout reached
                                self.process.kill()
                                print('[Error]: Tests for {} took too long. Process Killed.'.format(application), file=sys.stderr)
                                self.exitCodes[application] = 1
                            else:
                                if process_stdout:
                                    print(process_stdout.decode('utf8'), file=sys.stdout)
                                if process_stderr:
                                    print(process_stderr.decode('utf8'), file=sys.stderr)

                            # Running out of time in the tests will send the error code -15. We may want to skip
                            # that one in a future. Right now will throw everything different from 0.
                            self.exitCodes[application] = int(self.process.returncode != 0)
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
                                    script),
                                file=sys.stderr)
                            sys.stderr.flush()

    def RunCppTests(self, applications, verbosity = 1):
        ''' Calls the cpp tests directly
        '''

        # Iterate over all executables that are not mpi dependant and execute them.
        for test_suite in os.listdir(os.path.join(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()), "test")):
            filename = os.fsdecode(test_suite)

            # Skip mpi tests
            if ("MPI" not in filename and self.TestToAppName(filename) in applications) or filename == "KratosCoreTest":
                
                # Run all the tests in the executable
                self._RunTest(
                    test_suit_name=filename, 
                    level="nightly", 
                    verbose=1, 
                    command=[
                        os.path.join(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()),"test",filename)
                    ], 
                    timer=180)


    def RunMPIPythonTests(self, applications, mpi_command, mpi_flags, num_processes_flag, num_processes, level, verbose, command, timer):

        fullApplicationList = kratos_utils.GetListOfAvailableApplications()

        # If no applications are selected by the user, run all applications
        if not applications:
            applications = fullApplicationList

        # Iterate over the list of all applications an execute the ones selected by the user
        for application in fullApplicationList:
            if application in applications:

                test_script = Path(Kratos.KratosPaths.kratos_applications) / application / Path("tests") / Path("test_{}.py".format(application + "_mpi"))

                print(application, test_script, Path.is_file(test_script))

                if Path.is_file(test_script):
                    full_command = "{} {} {} {} {} {} --using-mpi -v{} -l{}".format(mpi_command, mpi_flags, num_processes_flag, num_processes, command, test_script, verbose, level)
                    try:
                        self.process = subprocess.Popen([
                            full_command
                        ], shell=True,
                        stdout=subprocess.PIPE,
                        cwd=os.path.dirname(os.path.abspath(str(test_script))))
                    except:
                        print('[Error]: Unable to execute "{}"'.format(full_command), file=sys.stderr)
                        self.exitCodes[application] = 1
                    else:
                        # Used instead of wait to "soft-block" the process and prevent deadlocks
                        # and capture the first exit code different from OK
                        try:
                            process_stdout, process_stderr = self.process.communicate(timeout=timer)
                        except subprocess.TimeoutExpired:
                            # Timeout reached
                            self.process.kill()
                            print('[Error]: Tests for {} took too long. Process Killed.'.format(application), file=sys.stderr)
                            self.exitCodes[application] = 1
                        else:
                            if process_stdout:
                                print(process_stdout.decode('utf8'), file=sys.stdout)
                            if process_stderr:
                                print(process_stderr.decode('utf8'), file=sys.stderr)

                        # Running out of time in the tests will send the error code -15. We may want to skip
                        # that one in a future. Right now will throw everything different from 0.
                        self.exitCodes[application] = int(self.process.returncode != 0)
                else:
                    if verbose > 0:
                        print('[Warning]: No test script found for {}'.format(application), file=sys.stderr, flush=True)
                    if verbose > 1:
                        print('  expected file: "{}"'.format(test_script), file=sys.stderr, flush=True)

    def RunMPICppTests(self, applications, verbosity, mpi_command, mpi_flags, num_processes_flag, num_processes):
        ''' Calls the mpi cpp tests directly
        '''

        # Iterate over all executables that are mpi dependant and execute them.
        for test_suite in os.listdir(os.path.join(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()), "test")):
            filename = os.fsdecode(test_suite)
            
            # Skip non-mpi tests
            if ("MPI" in filename and self.MPITestToAppName(filename) in applications) or filename == "KratosMPICoreTest":

                # Run all the tests in the executable
                self._RunTest(
                    test_suit_name=filename, 
                    level="nightly", 
                    verbose=verbosity, 
                    command=filter(None, [
                        mpi_command, 
                        mpi_flags, 
                        num_processes_flag, 
                        str(num_processes), 
                        os.path.join(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()),"test",filename)]
                    ),
                    timer=180)
