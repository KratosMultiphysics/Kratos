import sys
import os
import math

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsSettlementTests(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    @KratosUnittest.skip("unit test skipped in this branch, "
                         "as for some reason, the velocity vector is required in this quasi static test")
    def test_Abc_1_1_0_True_Deformations(self):
        test_name = 'test_Abc_1_1_0_True_Deformations'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        test_helper.run_stages(file_path, 2)
        node = 1
        cwd = os.getcwd()
        os.chdir(file_path)
        kratos_res = getTimesDisplacement(node)
        self.assertEqual(128, len(kratos_res))
        analytical_res = getComparisonResult("analytical.csv")

        # convert to objective stress
        analytical_error = compareResult(kratos_res, convertToObjectiveStress(analytical_res, 100))

        self.assertLess(analytical_error, 7.5,
                        "Analytical Comparison Failed (Av Diff >7.5%): {0:.2f}%".format(analytical_error))

        os.chdir(cwd)

def getComparisonResult(file):
    with open(os.path.join('.', file), 'r') as fo:
        fo.readline() # header
        timeDisplacement = []
        for line in fo:
            cells = line.split(";")
            timeDisplacement.append([float(cells[1]), float(cells[2])])
    return timeDisplacement

def interpolate(baseData, xValue):
    # This assumes a sorted increasing order of x_values
    for value in baseData:
        if value[0] <= xValue:
            x0 = value[0]
            y0 = value[1]
        if value[0] >= xValue:
            x1 = value[0]
            y1 = value[1]
            break

    if (x1 == x0):
        return y1

    return y0 * (1 - (xValue - x0) / (x1 - x0)) + y1 * (1 - (x1 - xValue) / (x1 - x0))

def getTimesDisplacement(node):
    timeDisplacement =[]
    with open(os.path.join('.', "mesh1.post.res"), 'r') as fo:
        for line in fo:
            if line.startswith('Result "DISPLACEMENT" "Kratos"'):
                time = float(line.split()[3])
                if time < 0.0:
                    continue
                fo.readline() # header
                for _ in range(node):
                    subline = fo.readline()
                    yDisp = float(subline.split()[2])
                    timeDisplacement.append([time, yDisp])
    return timeDisplacement

def compareResult(kratos_res, comp_res):
    ErrorRelativePercent = []
    for point in kratos_res:
        if point[0] < 0 + 2.5 or point[0] > 10000 - 2.5:
            continue
        y = interpolate(comp_res, point[0])
        if not math.isnan(y) and y != 0:
            ErrorRelativePercent.append(abs(100 * y/point[1] - 100))
            print(y, point[1], point[0], ErrorRelativePercent[-1])
    return sum(ErrorRelativePercent)/len(ErrorRelativePercent)

def convertToObjectiveStress(result, depth):
    # ABC
    # S = H0 * (1 - exp(-eps(t))) => S = H0 * (eps(t) / (1 + eps(t))

    objectiveStress = [];
    for point in result:
        y = -point[1]
        et = -math.log(1 - (y / depth))
        s = -depth * (et / (1 + et))
        objectiveStress.append([point[0], s])

    return objectiveStress


if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    small_suite = suites['small'] # These tests are executed by the continuous integration tool
    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsSettlementTests]))
    all_suite = suites['all']
    all_suite.addTests(small_suite)
    KratosUnittest.runTests(suites)
