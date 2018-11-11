from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestSearchMPMParticle(KratosUnittest.TestCase):

    def _generate_particle_element(self, current_model, dimension, geometry_element, is_structured, is_fine=False):
        # KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # Initialize model part
        ## Material model part definition
        material_model_part = current_model.CreateModelPart("dummy_name")
        material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        ## Initial material model part definition
        initial_material_model_part = current_model.CreateModelPart("Initial_dummy_name")
        initial_material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        ## Grid model part definition
        grid_model_part = current_model.CreateModelPart("Background_Grid")
        grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        # Create Background Grid
        sub_background = grid_model_part.CreateSubModelPart("test")
        if is_structured:
            self._create_nodes_structured(sub_background, dimension, geometry_element)
        else:
            self._create_nodes_unstructured(sub_background, dimension, geometry_element, is_fine)

        self._create_elements(sub_background,dimension, geometry_element)

        # Create element and nodes
        sub_mp = initial_material_model_part.CreateSubModelPart("test")
        if is_structured:
            self._create_nodes_structured(sub_mp, dimension, geometry_element)
        else:
            self._create_nodes_unstructured(sub_mp, dimension, geometry_element, is_fine)

        self._create_elements(sub_mp,dimension, geometry_element)

        # Initialize linear_solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()

        # Initialize element
        if geometry_element == "Triangle":
            if (dimension == 2):
                new_element = KratosParticle.CreateUpdatedLagragian2D3N()
            else:
                new_element = KratosParticle.CreateUpdatedLagragian3D4N()
        elif geometry_element == "Quadrilateral":
            if (dimension == 2):
                new_element = KratosParticle.CreateUpdatedLagragian2D4N()
            else:
                new_element = KratosParticle.CreateUpdatedLagragian3D8N()

        # Initialize solver
        if(dimension==2):
            self.solver = KratosParticle.MPM2D(grid_model_part, initial_material_model_part, material_model_part, linear_solver, new_element, False, "static", geometry_element, 1, False, False)
        else:
            self.solver = KratosParticle.MPM3D(grid_model_part, initial_material_model_part, material_model_part, linear_solver, new_element, False, "static", geometry_element, 1, False, False)


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

    def _create_nodes_unstructured(self, model_part, dimension, geometry_element, is_fine):
        if is_fine:
            modulus=1.e-7
        else:
            modulus=1

        if geometry_element == "Triangle":
            model_part.CreateNewNode(1, -2.3750791794090045*modulus, 1.31342           *modulus, 0.0*modulus)
            model_part.CreateNewNode(2, -1.2908699999999997*modulus, 2.3976291794090043*modulus, 0.0*modulus)
            model_part.CreateNewNode(3, -1.2206021540087075*modulus, 1.2431521540087074*modulus, 0.0*modulus)
            if (dimension == 3):
                model_part.CreateNewNode(4, -1.5782651431694477*modulus, 1.595285586264502*modulus, 0.6883104121306771*modulus)
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

    def _create_elements(self, model_part, dimension, geometry_element):
        if geometry_element == "Triangle":
            if (dimension == 2):
                model_part.CreateNewElement("UpdatedLagrangian2D3N", 1, [1,2,3], model_part.GetProperties()[1])
            if (dimension == 3):
                model_part.CreateNewElement("UpdatedLagrangian3D4N", 1, [1,2,3,4], model_part.GetProperties()[1])
        elif geometry_element == "Quadrilateral":
            if (dimension == 2):
                model_part.CreateNewElement("UpdatedLagrangian2D4N", 1, [1,2,3,4], model_part.GetProperties()[1])
            if (dimension == 3):
                model_part.CreateNewElement("UpdatedLagrangian3D8N", 1, [1,2,3,4,5,6,7,8], model_part.GetProperties()[1])

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, model_part.Elements)

    def _search_element_and_check(self,current_model):
        # Get model part
        material_model_part = current_model.GetModelPart("dummy_name")
        grid_model_part     = current_model.GetModelPart("Background_Grid")

        # Check before search
        # for mpm in material_model_part.Elements:
        #     KratosMultiphysics.Logger.PrintWarning(mpm.Id)

        # Search element
        self.solver.SearchElement(grid_model_part, material_model_part)

        # Check after search
        # for mpm in material_model_part.Elements:
        #     KratosMultiphysics.Logger.PrintWarning(mpm.Id)

    def test_SearchMPMParticleTriangle2DStructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Triangle", is_structured=True)
        self._search_element_and_check(current_model)

    def test_SearchMPMParticleTriangle3DStructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Triangle", is_structured=True)

    def test_SearchMPMParticleQuadrilateral2DStructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Quadrilateral", is_structured=True)

    def test_SearchMPMParticleQuadrilateral3DStructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Quadrilateral", is_structured=True)

    def test_SearchMPMParticleTriangle2DUnstructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Triangle", is_structured=False)

    def test_SearchMPMParticleTriangle3DUnstructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Triangle", is_structured=False)

    def test_SearchMPMParticleQuadrilateral2DUnstructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Quadrilateral", is_structured=False)

    def test_SearchMPMParticleQuadrilateral3DUnstructured(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Quadrilateral", is_structured=False)

    def test_SearchMPMParticleTriangle2DUnstructuredFine(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Triangle", is_structured=False, is_fine=True)

    def test_SearchMPMParticleTriangle3DUnstructuredFine(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Triangle", is_structured=False, is_fine=True)

    def test_SearchMPMParticleQuadrilateral2DUnstructuredFine(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Quadrilateral", is_structured=False, is_fine=True)

    def test_SearchMPMParticleQuadrilateral3DUnstructuredFine(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Quadrilateral", is_structured=False, is_fine=True)


if __name__ == '__main__':
    KratosUnittest.main()
