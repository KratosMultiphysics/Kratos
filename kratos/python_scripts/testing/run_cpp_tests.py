import os
import sys
import subprocess

import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils

for test_suite in os.listdir(os.path.join(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()), "test")):
    filename = os.fsdecode(test_suite)
    print(f"Running tests for {filename} ...")

    # Skip mpi tests
    if "MPI" in filename:
        pass
    else:
        # Run all the tests in the executable
        process = subprocess.Popen([
            os.path.join(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()),"test",filename)
        ], stdout=subprocess.PIPE)

        # Used instead of wait to "soft-block" the process and prevent deadlocks
        # and capture the first exit code different from OK
        try:
            # timeout should not be a problem for cpp, but we leave it just in case
            timer = int(90)
            process_stdout, process_stderr = process.communicate(timeout=timer)
        except subprocess.TimeoutExpired:
            # Timeout reached
            process.kill()
            print('[Error]: Tests for {} took too long. Process Killed.'.format(application), file=sys.stderr)
            exitCode = 1
        else:
            if process_stdout:
                print(process_stdout.decode('utf8'), file=sys.stdout)
            if process_stderr:
                print(process_stderr.decode('utf8'), file=sys.stderr)

        # Running out of time in the tests will send the error code -15. We may want to skip
        # that one in a future. Right now will throw everything different from 0.
        exitCode = int(process.returncode != 0)

