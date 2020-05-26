from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestSearchMPMParticle(KratosUnittest.TestCase):

    def _generate_particle_element(self, current_model, dimension, geometry_element, is_structured, is_fine=False):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # Initialize model part
        ## Material model part definition
        material_point_model_part = current_model.CreateModelPart("dummy_name")
        material_point_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)
        self.process_info = material_point_model_part.ProcessInfo

        ## Initial material model part definition
        initial_mesh_model_part = current_model.CreateModelPart("Initial_dummy_name")
        initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        ## Grid model part definition
        grid_model_part = current_model.CreateModelPart("Background_Grid")
        grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        # Create Background Grid
        sub_background = grid_model_part.CreateSubModelPart("test")
        if is_structured:
            self._create_background_nodes_structured(sub_background, dimension, geometry_element)
        else:
            self._create_background_nodes_unstructured(sub_background, dimension, geometry_element, is_fine)

        self._create_background_elements(sub_background,dimension, geometry_element, is_structured)

        # Create element and nodes
        sub_mp = initial_mesh_model_part.CreateSubModelPart("test")
        sub_mp.GetProperties()[1].SetValue(KratosParticle.PARTICLES_PER_ELEMENT, 1)
        if is_structured:
            self._create_nodes_structured(sub_mp, dimension, geometry_element)
        else:
            self._create_nodes_unstructured(sub_mp, dimension, geometry_element, is_fine)

        self._create_elements(sub_mp,dimension, geometry_element)

        # Set active
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mesh_model_part.Elements)

        # Generate MP Elements
        KratosParticle.GenerateMaterialPointElement(grid_model_part, initial_mesh_model_part, material_point_model_part, False)


    def _create_nodes_structured(self, model_part, dimension, geometry_element):
        if geometry_element == "Triangle":
            model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
            model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
            if (dimension == 3):
                model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
        elif geometry_element == "Quadrilateral":
            model_part.CreateNewNode(1, -0.5, -0.5, 0.0)
            model_part.CreateNewNode(2,  0.5, -0.5, 0.0)
            model_part.CreateNewNode(3,  0.5,  0.5, 0.0)
            model_part.CreateNewNode(4, -0.5,  0.5, 0.0)
            if (dimension == 3):
                model_part.CreateNewNode(5, -0.5, -0.5, 1.0)
                model_part.CreateNewNode(6,  0.5, -0.5, 1.0)
                model_part.CreateNewNode(7,  0.5,  0.5, 1.0)
                model_part.CreateNewNode(8, -0.5,  0.5, 1.0)


    def _create_background_nodes_structured(self, model_part, dimension, geometry_element):
        self._create_nodes_structured(model_part, dimension, geometry_element)
        if geometry_element == "Triangle":
            model_part.CreateNewNode(5, 1.0, 1.0, 0.0)
            if (dimension == 3):
                model_part.CreateNewNode(6, 1.0, 0.0, 1.0)
                model_part.CreateNewNode(7, 0.0, 1.0, 1.0)
                model_part.CreateNewNode(8, 1.0, 1.0, 1.0)
        elif geometry_element == "Quadrilateral":
            model_part.CreateNewNode(9 , 1.5, -0.5, 0.0)
            model_part.CreateNewNode(10, 1.5,  0.5, 0.0)
            if (dimension == 3):
                model_part.CreateNewNode(11, 1.5, -0.5, 1.0)
                model_part.CreateNewNode(12, 1.5,  0.5, 1.0)


    def _create_nodes_unstructured(self, model_part, dimension, geometry_element, is_fine):
        if is_fine:
            modulus=1.e-7
        else:
            modulus=1

        if geometry_element == "Triangle":
            model_part.CreateNewNode(1, 0.9*modulus, 2.9*modulus, 0.0*modulus)
            model_part.CreateNewNode(2, 0.0*modulus, 0.0*modulus, 0.0*modulus)
            model_part.CreateNewNode(3, 2.9*modulus, 2.4*modulus, 0.0*modulus)
            if (dimension == 3):
                model_part.CreateNewNode(4, 0.5*modulus, 0.5*modulus, 1.5*modulus)
        elif geometry_element == "Quadrilateral":
            model_part.CreateNewNode(1, -0.734*modulus, -0.621*modulus, 0.0*modulus)
            model_part.CreateNewNode(2,  0.497*modulus, -0.432*modulus, 0.0*modulus)
            model_part.CreateNewNode(3,  0.587*modulus,  0.402*modulus, 0.0*modulus)
            model_part.CreateNewNode(4, -0.809*modulus,  0.522*modulus, 0.0*modulus)
            if (dimension == 3):
                model_part.CreateNewNode(5, -0.621*modulus, -0.734*modulus, 0.93*modulus)
                model_part.CreateNewNode(6,  0.432*modulus, -0.497*modulus, 1.11*modulus)
                model_part.CreateNewNode(7,  0.402*modulus,  0.587*modulus, 0.89*modulus)
                model_part.CreateNewNode(8, -0.522*modulus,  0.809*modulus, 1.21*modulus)


    def _create_background_nodes_unstructured(self, model_part, dimension, geometry_element, is_fine):
        self._create_nodes_unstructured(model_part, dimension, geometry_element, is_fine)

        if is_fine:
            modulus=1.e-7
        else:
            modulus=1

        if geometry_element == "Triangle":
            model_part.CreateNewNode(5, 2.1*modulus, 0.1*modulus, 0.0*modulus)
        elif geometry_element == "Quadrilateral":
            model_part.CreateNewNode(9,  1.343*modulus, -0.451*modulus, 0.0*modulus)
            model_part.CreateNewNode(10, 1.512*modulus,  0.392*modulus, 0.0*modulus)
            if (dimension == 3):
                model_part.CreateNewNode(11, 1.742*modulus, -0.620*modulus, 0.999*modulus)
                model_part.CreateNewNode(12, 1.520*modulus,  0.671*modulus, 1.120*modulus)


    def _create_elements(self, model_part, dimension, geometry_element):
        if geometry_element == "Triangle":
            if (dimension == 2):
                model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])
            if (dimension == 3):
                model_part.CreateNewElement("Element3D4N", 1, [1,2,3,4], model_part.GetProperties()[1])
        elif geometry_element == "Quadrilateral":
            if (dimension == 2):
                model_part.CreateNewElement("Element2D4N", 1, [1,2,3,4], model_part.GetProperties()[1])
            if (dimension == 3):
                model_part.CreateNewElement("Element3D8N", 1, [1,2,3,4,5,6,7,8], model_part.GetProperties()[1])


    def _create_background_elements(self, model_part, dimension, geometry_element, is_structured):
        self._create_elements(model_part, dimension, geometry_element)
        if geometry_element == "Triangle":
            if (dimension == 2):
                model_part.CreateNewElement("Element2D3N", 2, [2,3,5], model_part.GetProperties()[1])
            if (dimension == 3):
                if (is_structured):
                    model_part.CreateNewElement("Element3D4N", 2, [2,8,4,6], model_part.GetProperties()[1])
                    model_part.CreateNewElement("Element3D4N", 3, [4,8,3,7], model_part.GetProperties()[1])
                    model_part.CreateNewElement("Element3D4N", 4, [2,5,3,8], model_part.GetProperties()[1])
                    model_part.CreateNewElement("Element3D4N", 5, [8,3,2,4], model_part.GetProperties()[1])
                else:
                    model_part.CreateNewElement("Element3D4N", 2, [2,3,5,4], model_part.GetProperties()[1])
        elif geometry_element == "Quadrilateral":
            if (dimension == 2):
                model_part.CreateNewElement("Element2D4N", 2, [2,9,10,3], model_part.GetProperties()[1])
            if (dimension == 3):
                model_part.CreateNewElement("Element3D8N", 2, [2,9,10,3,6,11,12,7], model_part.GetProperties()[1])


    def _move_and_search_element(self, current_model, new_coordinate, max_num_results = 1000, specific_tolerance = 1.e-5):
        # Get model part
        material_point_model_part = current_model.GetModelPart("dummy_name")
        grid_model_part           = current_model.GetModelPart("Background_Grid")

        # Apply before search
        for mpm in material_point_model_part.Elements:
            mpm.SetValuesOnIntegrationPoints(KratosParticle.MP_COORD, [new_coordinate], self.process_info)

        # Search element
        KratosParticle.SearchElement(grid_model_part, material_point_model_part, max_num_results, specific_tolerance)

    def _check_connectivity(self, current_model, expected_connectivity_node=[]):
        # Get model part
        material_point_model_part = current_model.GetModelPart("dummy_name")
        grid_model_part           = current_model.GetModelPart("Background_Grid")

        # Check the searched node as expected connectivity
        if not expected_connectivity_node:
            for mpm in material_point_model_part.Elements:
                self.assertEqual(mpm.GetNodes(), [])
        else:
            for mpm in material_point_model_part.Elements:
                for i in range (len(expected_connectivity_node)):
                    self.assertEqual(mpm.GetNode(i).Id, grid_model_part.GetNode(expected_connectivity_node[i]).Id)
                    self.assertEqual(mpm.GetNode(i).X, grid_model_part.GetNode(expected_connectivity_node[i]).X)
                    self.assertEqual(mpm.GetNode(i).Y, grid_model_part.GetNode(expected_connectivity_node[i]).Y)
                    self.assertEqual(mpm.GetNode(i).Z, grid_model_part.GetNode(expected_connectivity_node[i]).Z)

    def test_SearchMPMParticleTriangle2DStructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Triangle", is_structured=True)

        new_coordinate = [0.5, 0.5, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3])

        new_coordinate = [0.50001, 0.50001, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,3,5])

        new_coordinate = [1.00001, 1.00001, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMPMParticleTriangle3DStructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Triangle", is_structured=True)

        new_coordinate = [0.5, 0.25, 0.20]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4])

        new_coordinate = [0.90, 0.55, 0.90]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,8,4,6])

        new_coordinate = [0.10, 0.90, 0.55]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [4,8,3,7])

        new_coordinate = [0.90, 0.90, 0.55]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,5,3,8])

        new_coordinate = [0.50, 0.50, 0.50]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [8,3,2,4])

        new_coordinate = [1.0001, 1.0001, 1.0001]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMPMParticleQuadrilateral2DStructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Quadrilateral", is_structured=True)

        new_coordinate = [-0.11111, 0.12345, 1.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4])

        new_coordinate = [0.6, 0.12345, 1.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,9,10,3])

        new_coordinate = [1.00001, 1.00001, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMPMParticleQuadrilateral3DStructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Quadrilateral", is_structured=True)

        new_coordinate = [0.5, 0.25, 0.20]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4,5,6,7,8])

        new_coordinate = [0.7, 0.35, 0.3]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,9,10,3,6,11,12,7])

        new_coordinate = [0.50001, 0.50001, 0.50001]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMPMParticleTriangle2DUnstructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Triangle", is_structured=False)

        new_coordinate = [1.31967, 1.85246, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3])

        new_coordinate = [1.72951, 0.491803, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,3,5])

        new_coordinate = [3.00001, 3.00001, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)


    def test_SearchMPMParticleTriangle3DUnstructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Triangle", is_structured=False)

        new_coordinate = [1.31967, 1.85246, 0.1]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4])

        new_coordinate = [1.72951, 0.491803, 0.1]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,3,5,4])

        new_coordinate = [3.00001, 3.00001, 1.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMPMParticleQuadrilateral2DUnstructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Quadrilateral", is_structured=False)

        new_coordinate = [-0.11111, 0.12345, 1.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4])

        new_coordinate = [0.6, 0.12345, 1.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,9,10,3])

        new_coordinate = [1.00001, 1.00001, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMPMParticleQuadrilateral3DUnstructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Quadrilateral", is_structured=False)

        new_coordinate = [0.5, 0.25, 0.20]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4,5,6,7,8])

        new_coordinate = [0.7, 0.35, 0.3]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,9,10,3,6,11,12,7])

        new_coordinate = [0.70001, 0.20001, 1.20001]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMPMParticleTriangle2DUnstructuredFine(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Triangle", is_structured=False, is_fine=True)

        new_coordinate = [1.31967e-7, 1.85246e-7, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3])

        new_coordinate = [1.72951e-7, 0.491803e-7, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,3,5])

        new_coordinate = [3.00001e-7, 3.00001e-7, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMPMParticleTriangle3DUnstructuredFine(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Triangle", is_structured=False, is_fine=True)

        new_coordinate = [1.31967e-7, 1.85246e-7, 1.0e-8]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4])

        new_coordinate = [1.72951e-7, 0.491803e-7, 0.1e-7]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,3,5,4])

        new_coordinate = [3.00001e-7, 3.00001e-7, 1.0e-7]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMPMParticleQuadrilateral2DUnstructuredFine(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Quadrilateral", is_structured=False, is_fine=True)

        new_coordinate = [-0.11111e-7, 0.12345e-7, 1.0e-7]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4])

        new_coordinate = [0.6e-7, 0.12345e-7, 1.0e-7]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,9,10,3])

        new_coordinate = [1.00001e-7, 1.00001e-7, 0.0]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMPMParticleQuadrilateral3DUnstructuredFine(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Quadrilateral", is_structured=False, is_fine=True)

        new_coordinate = [0.5e-7, 0.25e-7, 0.20e-7]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4,5,6,7,8])

        new_coordinate = [0.7e-7, 0.35e-7, 0.3e-7]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,9,10,3,6,11,12,7])

        new_coordinate = [0.70001e-7, 0.20001e-7, 1.20001e-7]
        self._move_and_search_element(current_model, new_coordinate)
        self._check_connectivity(current_model)

if __name__ == '__main__':
    KratosUnittest.main()
