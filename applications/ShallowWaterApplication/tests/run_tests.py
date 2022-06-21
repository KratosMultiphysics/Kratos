import argparse
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication
import test_ShallowWaterApplication as python_tests

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--level', default='all', choices=['all', 'cpp', 'nightly', 'small', 'validation'])
parser.add_argument('-v', '--verbosity', default=2, type=int, choices=[0, 1, 2])
parser.add_argument('--timing', action='store_true')
parser.add_argument('--using-mpi', action='store_true')
args = parser.parse_args() # Note: the options are redefined to avoid conflict with KratosUnittest

def run_cpp_tests():
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.FAILED_TESTS_OUTPUTS)
    KratosMultiphysics.Tester.RunTestSuite("ShallowWaterApplicationFastSuite")

try:
    run_cpp_tests()
    cpp = 'OK'
except RuntimeError:
    cpp = 'Failed'

if not args.level == 'cpp':
    try:
        python_tests.run()
    except (RuntimeError, SystemExit) as e:
        if e.code:
            py = 'Failed'
        else:
            py = 'OK'

print()
print('ShallowWaterApplication:')
print('cpp tests ..........', cpp)
if not args.level == 'cpp':
    print('python tests .......', py)
