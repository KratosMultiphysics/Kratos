import argparse
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication
import test_ShallowWaterApplication as python_tests

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--level', default='all', choices=['all', 'cpp', 'small', 'nightly', 'validation'])
parser.add_argument('-v', '--verbosity', default=2, type=int, choices=[0, 1, 2])
parser.add_argument('--timing', action='store_true')
parser.add_argument('--using-mpi', action='store_true')
parser.add_argument('-t', '--cpp_test', type=str) # This option is incompatible with KratosUnittest
args = parser.parse_args() # Note: the options are redefined to avoid conflict with KratosUnittest

if args.cpp_test and args.level != 'cpp':
    args.level = 'cpp'
    KratosMultiphysics.Logger.PrintInfo('[run_tests] parse arguments','Setting "-level" to "cpp"')

cpp_verbosity = {
    0 : KratosMultiphysics.Tester.Verbosity.QUITE,
    1 : KratosMultiphysics.Tester.Verbosity.PROGRESS,
    2 : KratosMultiphysics.Tester.Verbosity.FAILED_TESTS_OUTPUTS,
}

def run_cpp_tests():
    KratosMultiphysics.Tester.SetVerbosity(cpp_verbosity[args.verbosity])
    if args.cpp_test:
        KratosMultiphysics.Tester.RunTestCases(args.cpp_test)
    else:
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

if args.verbosity > 0 and args.level != 'cpp':
    print()
    print('ShallowWaterApplication:')
    print('cpp tests ..........', cpp)
    if not args.level == 'cpp':
        print('python tests .......', py)
