import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import os
try:
    from mappers.test_nearest import Case1D, Case2D, Case3DSphere, Case3DSinc
except:
    from test_nearest import Case1D, Case2D, Case3DSphere, Case3DSinc


class TestMapperRadialBasis(KratosUnittest.TestCase):

    def test_mapper_radial_basis(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_radial_basis.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())
        par_mapper = parameters['mapper']

        gui = 1  # *** gui gives problems when running all tests?

        # 1D case: square-root grid + linear function
        """
        n_from = 14, n_to = 5 
            => max error = 5.8e-5
        """
        n_from, n_to = 14, 5
        case = Case1D(cs_data_structure, n_from, n_to)

        par_mapper['settings'].SetArray('directions', ['Z'])
        mapper = cs_tools.CreateInstance(par_mapper)
        mapper.Initialize(case.model_part_from, case.model_part_to)
        mapper((case.model_part_from, case.var_from),
               (case.model_part_to, case.var_to))

        self.assertTrue(case.check(tolerance=1e-4))

        if gui:
            case.plot()

        # 2D case: circle + linear function
        """
        n_from = 33, n_to = 22 
            => max error = 0.003
        """
        n_from, n_to = 33, 22
        case = Case2D(cs_data_structure, n_from, n_to)

        par_mapper['settings'].SetArray('directions', ['X', 'Y'])
        mapper = cs_tools.CreateInstance(par_mapper)
        mapper.Initialize(case.model_part_from, case.model_part_to)
        mapper((case.model_part_from, case.var_from),
               (case.model_part_to, case.var_to))

        self.assertTrue(case.check(tolerance=0.005))

        if gui:
            case.plot()

        # 3D case: sphere + sine function
        """
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
            => max error = 3.2e-4
        """
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
        case = Case3DSphere(cs_data_structure, n_theta_from, n_phi_from, n_theta_to, n_phi_to)

        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])
        mapper = cs_tools.CreateInstance(par_mapper)
        mapper.Initialize(case.model_part_from, case.model_part_to)
        mapper((case.model_part_from, case.var_from),
               (case.model_part_to, case.var_to))

        self.assertTrue(case.check(tolerance=5e-4))

        if gui:
            case.plot()

        # 3D case: sinc + linear vector function
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
            => max error = 0.05
            
        n_x_from, n_y_from = 50, 50
        n_x_to, n_y_to = 60, 60
            => max error = 0.0024
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
        case = Case3DSinc(cs_data_structure, n_x_from, n_y_from, n_x_to, n_y_to)

        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])
        mapper = cs_tools.CreateInstance(par_mapper)
        mapper.Initialize(case.model_part_from, case.model_part_to)
        mapper((case.model_part_from, case.var_from),
               (case.model_part_to, case.var_to))

        for tmp in case.check(tolerance=0.1):
            self.assertTrue(tmp)

        if gui:
            case.plot()


if __name__ == '__main__':
    KratosUnittest.main()
