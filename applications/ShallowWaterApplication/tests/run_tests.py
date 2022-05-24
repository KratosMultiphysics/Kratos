import argparse
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication
import test_ShallowWaterApplication as python_tests

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--level', default='all', choices=['all', 'cpp', 'nightly', 'small', 'validation'])
args = parser.parse_args()

try:
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.FAILED_TESTS_OUTPUTS)
    KratosMultiphysics.Tester.RunTestSuite("ShallowWaterApplicationFastSuite")
    cpp = 'OK'
except:
    cpp = 'Failed'

if not args.level == 'cpp':
    try:
        python_tests.run()
        py = 'OK'
    except:
        py = 'Failed'

print()
print('ShallowWaterApplication tests:')
print('cpp tests ..........', cpp)
if not args.level == 'cpp':
    print('python tests .......', py)
