import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestParticleCreatorDestructor(KratosUnittest.TestCase):

    def setUp(self):
        self.current_model = Kratos.Model()
        self.spheres_model_part = self.current_model.CreateModelPart("SpheresPart")
        self.spheres_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.spheres_model_part.AddNodalSolutionStepVariable(Kratos.ANGULAR_VELOCITY)
        self.spheres_model_part.AddNodalSolutionStepVariable(Kratos.RADIUS)
        self.spheres_model_part.AddNodalSolutionStepVariable(Kratos.PARTICLE_MATERIAL)
        self.spheres_model_part.AddNodalSolutionStepVariable(Kratos.NODAL_MASS)

        properties = Kratos.Properties(0)
        properties[Kratos.YOUNG_MODULUS] = 3.331
        properties[DEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEM_D_Hertz_viscous_Coulomb"
        properties[DEM.FRICTION] = 0.0
        properties[Kratos.POISSON_RATIO] = 0.0
        properties[DEM.DAMPING_GAMMA] = 0.0
        properties[DEM.COEFFICIENT_OF_RESTITUTION] = 0.0

        self.ModifyProperties(properties)

        self.spheres_model_part.AddProperties(properties)

        DEM.PropertiesProxiesManager().CreatePropertiesProxies(self.spheres_model_part)

        self.creator_destructor = DEM.ParticleCreatorDestructor()

    def ModifyProperties(self, properties, param = 0):

        if not param:
            DiscontinuumConstitutiveLaw = getattr(DEM, properties[DEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME])()
            DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties, False)

        scheme = DEM.SymplecticEulerScheme()
        scheme.SetTranslationalIntegrationSchemeInProperties(properties, False)
        scheme.SetRotationalIntegrationSchemeInProperties(properties, False)

    def test_CreateSphericParticle1(self):

        coordinates = Kratos.Array3()
        coordinates[0] = 1.0
        coordinates[1] = 2.0
        coordinates[2] = 3.0

        properties = self.spheres_model_part.GetProperties()[0]
        radius = 0.1
        element_name = "SphericParticle3D"

        created_element = self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)
        created_node = created_element.GetNodes()[0]

        #Direct verification with the returned pointer
        self.assertEqual(created_element.Properties[Kratos.YOUNG_MODULUS], 3.331)
        self.assertEqual(created_node.X, 1.0)
        self.assertEqual(created_node.Y, 2.0)
        self.assertEqual(created_node.Z, 3.0)

        #Indirect verification looping over the ModelPart (also verifies that the element has been added to the ModelPart
        counter = 0
        for node in self.spheres_model_part.Nodes:
            counter += 1
            self.assertEqual(node.X, 1.0)
            self.assertEqual(node.Y, 2.0)
            self.assertEqual(node.Z, 3.0)

        self.assertEqual(counter, 1)

        counter = 0
        for element in self.spheres_model_part.Elements:
            counter += 1
            self.assertEqual(element.Properties[Kratos.YOUNG_MODULUS], 3.331)

        self.assertEqual(counter, 1)


    def test_CreateSphericParticle2(self):
        pass


if __name__ == '__main__':
    KratosUnittest.main()
