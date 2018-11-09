from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestGenerateMPMParticle(KratosUnittest.TestCase):

    def _generate_particle_element(self, current_model, dimension, geometry_element, num_particle):
        # Initialize element
        if geometry_element == "Triangle":
            if (dimension == 2):
                self.new_element = KratosParticle.CreateUpdatedLagragian2D3N()
            else:
                self.new_element = KratosParticle.CreateUpdatedLagragian3D4N()
        elif geometry_element == "Quadrilateral":
            if (dimension == 2):
                self.new_element = KratosParticle.CreateUpdatedLagragian2D4N()
            else:
                self.new_element = KratosParticle.CreateUpdatedLagragian3D8N()

        # Initialize solver
        # if(dimension==2):
        #     self.solver = KratosParticle.MPM2D(self.grid_model_part, self.initial_material_model_part, self.material_model_part, self.linear_solver, self.new_element, self.move_mesh_flag, "static", geometry_element, num_particle, False, False)
        # else:
        #     self.solver = KratosParticle.MPM3D(self.grid_model_part, self.initial_material_model_part, self.material_model_part, self.linear_solver, self.new_element, self.move_mesh_flag, "static", geometry_element, num_particle, False, False)

        pass

    def test_GenerateMPMParticleTriangle2D_1(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Triangle", num_particle=1)

    def test_GenerateMPMParticleTriangle2D_3(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Triangle", num_particle=3)

    def test_GenerateMPMParticleTriangle2D_6(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Triangle", num_particle=6)

    def test_GenerateMPMParticleTriangle2D_16(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Triangle", num_particle=16)

    def test_GenerateMPMParticleTriangle3D_1(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Triangle", num_particle=1)

    def test_GenerateMPMParticleTriangle3D_3(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Triangle", num_particle=3)

    def test_GenerateMPMParticleTriangle3D_6(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Triangle", num_particle=6)

    def test_GenerateMPMParticleTriangle3D_16(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Triangle", num_particle=16)

    def test_GenerateMPMParticleQuadrilateral2D_1(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=1)

    def test_GenerateMPMParticleQuadrilateral2D_4(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=4)

    def test_GenerateMPMParticleQuadrilateral2D_9(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=9)

    def test_GenerateMPMParticleQuadrilateral2D_16(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=16)

    def test_GenerateMPMParticleQuadrilateral3D_1(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=1)

    def test_GenerateMPMParticleQuadrilateral3D_4(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=4)

    def test_GenerateMPMParticleQuadrilateral3D_9(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=9)

    def test_GenerateMPMParticleQuadrilateral3D_16(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=16)

if __name__ == '__main__':
    KratosUnittest.main()
