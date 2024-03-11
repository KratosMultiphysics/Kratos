import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsConditionTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the regression on a previously obtained value.
    """
    etalon_value1 = -76430.423509075
    etalon_value2 = -0.0005
    etalon_value3 = 1356.95764908908

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_PwNormalFluxCondition2D2N(self):
        test_name = 'test_PwNormalFluxCondition2D2N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        pres = pressure[2]
        self.assertAlmostEqual(self.etalon_value1, pres)

    def test_PwNormalFluxCondition2D3N(self):
        test_name = 'test_PwNormalFluxCondition2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        pres = pressure[2]
        self.assertAlmostEqual(self.etalon_value1, pres)

    def test_PwNormalFluxCondition2D4N(self):
        test_name = 'test_PwNormalFluxCondition2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        pres = pressure[2]
        self.assertAlmostEqual(self.etalon_value1, pres)
        
    def test_PwNormalFluxCondition2D5N(self):
        test_name = 'test_PwNormalFluxCondition2D5N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        pres = pressure[2]
        self.assertAlmostEqual(self.etalon_value1, pres)
        
        
    def test_UPwNormalFluxCondition2D2N(self):
        test_name = 'test_UPwNormalFluxCondition2D2N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        pres = pressure[2]
        self.assertAlmostEqual(self.etalon_value3, pres)

    def test_UPwNormalFluxCondition2D3N(self):
        test_name = 'test_UPwNormalFluxCondition2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        pres = pressure[2]
        self.assertAlmostEqual(self.etalon_value3, pres)

    def test_UPwNormalFluxCondition2D4N(self):
        test_name = 'test_UPwNormalFluxCondition2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        pres = pressure[2]
        self.assertAlmostEqual(self.etalon_value3, pres)
        
    def test_UPwNormalFluxCondition2D5N(self):
        test_name = 'test_UPwNormalFluxCondition2D5N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        pres = pressure[2]
        self.assertAlmostEqual(self.etalon_value3, pres)
        
        
    def test_LineNormalFluidFluxDiffOrderCondition2D3N(self):
        test_name = 'test_LineNormalFluidFluxDiffOrderCondition2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        pres = pressure[2]
        self.assertAlmostEqual(self.etalon_value3, pres)

        
    def test_UPwFaceLoadCondition2D2N(self):
        test_name = 'test_UPwFaceLoadCondition2D2N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_UPwFaceLoadCondition2D3N(self):
        test_name = 'test_UPwFaceLoadCondition2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_UPwFaceLoadCondition2D4N(self):
        test_name = 'test_UPwFaceLoadCondition2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_UPwFaceLoadCondition2D5N(self):
        test_name = 'test_UPwFaceLoadCondition2D5N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
        
    def test_LineLoadDiffOrderCondition2D3N(self):
        test_name = 'test_LineLoadDiffOrderCondition2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_LineLoadDiffOrderCondition2D4N(self):
        test_name = 'test_LineLoadDiffOrderCondition2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_LineLoadDiffOrderCondition2D5N(self):
        test_name = 'test_LineLoadDiffOrderCondition2D5N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
        
    def test_UPwNormalFaceLoadCondition2D2N(self):
        test_name = 'test_UPwNormalFaceLoadCondition2D2N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_UPwNormalFaceLoadCondition2D3N(self):
        test_name = 'test_UPwNormalFaceLoadCondition2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_UPwNormalFaceLoadCondition2D4N(self):
        test_name = 'test_UPwNormalFaceLoadCondition2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_UPwNormalFaceLoadCondition2D5N(self):
        test_name = 'test_UPwNormalFaceLoadCondition2D5N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
        
    def test_LineNormalLoadDiffOrderCondition2D3N(self):
        test_name = 'test_LineNormalLoadDiffOrderCondition2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_LineNormalLoadDiffOrderCondition2D4N(self):
        test_name = 'test_LineNormalLoadDiffOrderCondition2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_LineNormalLoadDiffOrderCondition2D5N(self):
        test_name = 'test_LineNormalLoadDiffOrderCondition2D5N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
        
    def test_AxisymmetricUPwNormalFaceLoadCondition2D2N(self):
        test_name = 'test_AxisymmetricUPwNormalFaceLoadCondition2D2N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_AxisymmetricUPwNormalFaceLoadCondition2D3N(self):
        test_name = 'test_AxisymmetricUPwNormalFaceLoadCondition2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_AxisymmetricUPwNormalFaceLoadCondition2D4N(self):
        test_name = 'test_AxisymmetricUPwNormalFaceLoadCondition2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])
        
    def test_AxisymmetricUPwNormalFaceLoadCondition2D5N(self):
        test_name = 'test_AxisymmetricUPwNormalFaceLoadCondition2D5N'
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        disp = displacement[2]
        self.assertAlmostEqual(self.etalon_value2, disp[1])

if __name__ == '__main__':
    KratosUnittest.main()
