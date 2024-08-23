import os
import json
import numpy as np
import math

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper
from strip_load_semi_analytical_solution import StripLoad
from point_load_semi_analytical_solution import PointLoad


class KratosGeoMechanicsDynamicsLongTests(KratosUnittest.TestCase):
    """
    This class contains long-running tests to check dynamic calculations
    """

    def test_constant_strip_load_2d(self):
        """
        Tests a constant strip load on a 10m X 10m  block. The block is loaded with a constant strip load of 1kN/m and
        a length of 1m.
    
        The solution should roughly follow the semi-analytical solution presented in the following publication:
        An Introduction to Soil Dynamics , Verruijt A., 2009, Delft University of Technology, Chapter 12.2
    
    
        """
        # set to true if the results should be checked with the analytical solution
        CHECK_RESULTS = False
    
        test_name = 'test_constant_strip_load_2d'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
    
        simulation = test_helper.run_kratos(file_path)
    
        # get calculated results
        with open(os.path.join(file_path, "json_output.json")) as fp:
            calculated_result = json.load(fp)
    
        # get expected results
        with open(os.path.join(file_path, "expected_json_output.json")) as fp:
            expected_result = json.load(fp)
    
        # check if results are as expected
        self.assertTrue(test_helper.are_dictionaries_almost_equal(expected_result, calculated_result))
    
        # If this test is altered, results between the analytical solution and the calculated solution can be compared
        # below
        if CHECK_RESULTS:
            # get properties from material file as not all variables are registered in python
            with open(os.path.join(file_path, "MaterialParameters_stage_1.json")) as fp:
                material_parameters = json.load(fp)
    
            material_variables = material_parameters["properties"][0]["Material"]["Variables"]
            process_info = simulation.model.GetModelPart('PorousDomain.porous_computational_model_part').ProcessInfo
    
            top_surface = 10.0  # follows from the mesh
            line_load_length = 1.0  # follows from the mesh
            load_value = -1000.0  # follows from table in the mdpa file
            end_time = process_info.GetValue(Kratos.TIME)
    
            # initialise semi-analytical solution
            analytical_solution = StripLoad(material_variables["YOUNG_MODULUS"],
                                            material_variables["POISSON_RATIO"],
                                            (1 - material_variables["POROSITY"]) * material_variables["DENSITY_SOLID"],
                                            load_value)
    
            # get calculated json output
            json_model_part = simulation.model.GetModelPart('PorousDomain.json_output')
    
            elements = json_model_part.Elements
            x_coords = []
            calculated_vert_stresses = []
            analytical_vert_stresses = []
            for element in elements:
                # get centroid
                centroid = element.GetGeometry().Center()
                x_coord = centroid.X
                y_coord = centroid.Y
                depth_check = top_surface - y_coord
                x_coords.append(centroid.X)
    
                vert_stress_analytic = analytical_solution.calculate_vertical_stress(x_coord,
                                                                                     depth_check,
                                                                                     end_time,
                                                                                     line_load_length,
                                                                                     load_value)
    
                analytical_vert_stresses.append(-vert_stress_analytic)
    
                # get element mean vertical cauchy stress
                element_stress = element.CalculateOnIntegrationPoints(Kratos.CAUCHY_STRESS_VECTOR,
                                                                      json_model_part.ProcessInfo)
                vert_stress = test_helper.compute_mean_list([gauss_stress[1] for gauss_stress in element_stress])
    
                calculated_vert_stresses.append(vert_stress)
    
            # visualise results if matplotlib is installed
            try:
                import matplotlib.pyplot as plt
                plt.plot(x_coords, analytical_vert_stresses, 'o', color="r")
                plt.plot(x_coords, calculated_vert_stresses, 'o', color="b")
    
                plt.legend(["Analytical", "Calculated"])
    
                plt.show()
    
            except ImportError:
                print("Matplotlib not installed. Cannot visualise results")


    def constant_point_load(self, test_name):
        """
        Tests a constant point load on a 10m X 10m block in axisymmetric coordinates. The block is loaded with a constant point load of 1kN/m.

        The solution should roughly follow the semi-analytical solution presented in the following publication:
        An Introduction to Soil Dynamics , Verruijt A., 2009, Delft University of Technology, Chapter 13.2


        """
        # set to true if the results should be checked with the analytical solution
        CHECK_RESULTS = False

        file_path = test_helper.get_file_path(os.path.join('.', 'test_dynamic', test_name))

        simulation = test_helper.run_kratos(file_path)

        # get calculated results
        with open(os.path.join(file_path, "json_output.json")) as fp:
            calculated_result = json.load(fp)

        # get expected results
        with open(os.path.join(file_path, "expected_json_output.json")) as fp:
            expected_result = json.load(fp)

        # check if results are as expected
        self.assertTrue(test_helper.are_dictionaries_almost_equal(expected_result, calculated_result))

        # If this test is altered, results between the analytical solution and the calculated solution can be compared
        # below
        if CHECK_RESULTS:
            # get properties from material file as not all variables are registered in python
            with open(os.path.join(file_path, "MaterialParameters.json")) as fp:
                material_parameters = json.load(fp)

            material_variables = material_parameters["properties"][0]["Material"]["Variables"]
            process_info = simulation.model.GetModelPart('PorousDomain.porous_computational_model_part').ProcessInfo

            load_value = -1000.0  # follows from table in the mdpa file

            # initialise semi-analytical solution
            analytical_solution = PointLoad(material_variables["YOUNG_MODULUS"],
                                            material_variables["POISSON_RATIO"],
                                            (1 - material_variables["POROSITY"]) * material_variables["DENSITY_SOLID"],
                                            load_value)

            # get calculated json output
            json_model_part = simulation.model.GetModelPart('PorousDomain.json_output')

            nodes = json_model_part.Nodes
            calculated_vert_displacement = []
            analytical_vert_displacement = []
            time = calculated_result["TIME"]
            
            shear_modulus = material_variables["YOUNG_MODULUS"] / (2 * (1 + material_variables["POISSON_RATIO"]))
            density = (1 - material_variables["POROSITY"]) * material_variables["DENSITY_SOLID"]
            cs = math.sqrt(shear_modulus / density)
            
            for node in nodes:
                x_coord = node.X
                z_coord = node.Z
                radius = math.sqrt(x_coord * x_coord + z_coord * z_coord)
                time_normalized = np.asarray(time) * cs / radius
                node_id = node.Id
                vert_displ = calculated_result[f"NODE_{node_id}"]["DISPLACEMENT_Y"]
                vert_displ = - np.asarray(vert_displ) * shear_modulus * radius / load_value
                calculated_vert_displacement.append(vert_displ)
                
                for time_step in time:
                    vert_displ_analytic = analytical_solution.calculate_vertical_displacement(radius, time_step, load_value)
                    analytical_vert_displacement.append(-vert_displ_analytic)

            # visualise results if matplotlib is installed
            try:
                import matplotlib.pyplot as plt
                plt.plot(time_normalized, analytical_vert_displacement, '-', color="r")
                plt.plot(time_normalized, calculated_vert_displacement[0], '-', color="b")

                plt.legend(["Analytical", "Calculated"])
                plt.show()

            except ImportError:
                print("Matplotlib not installed. Cannot visualise results")


    def test_constant_point_load_2d_axi(self):
        self.constant_point_load('test_constant_point_load_2d_axi')
