import argparse
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication
import test_ShallowWaterApplication as python_tests

parser = argparse.ArgumentParser(description='run the ShallowWaterApplication tests. By default, both cpp and python tests are executed.')
options = parser.add_mutually_exclusive_group()
options.add_argument('-c', '--cpp', dest='py', action='store_false', help='run the cpp tests.')
options.add_argument('-p', '--py', dest='cpp', action='store_false', help='run the python tests.')
args = parser.parse_args()

print(args) # para ver que pasa

if args.cpp:
    try:
        KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.FAILED_TESTS_OUTPUTS)
        KratosMultiphysics.Tester.RunTestSuite("ShallowWaterApplicationFastSuite")
        cpp = 'OK'
    except:
        cpp = 'Failed'

if args.py:
    try:
        python_tests.run()
        py = 'OK'
    except:
        py = 'Failed'

print()
print('ShallowWaterApplication tests:')
if args.cpp:
    print('cpp tests ..........', cpp)
if args.py:
    print('python tests .......', py)
