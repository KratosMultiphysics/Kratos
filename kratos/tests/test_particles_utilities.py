
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import numpy as np

class TestParticlesUtilities(KratosUnittest.TestCase):

    def test_count_particles_in_nodes(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("Main")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        particles_mp = current_model.CreateModelPart("Particles")

        # 3 - 4
        # | / |
        # 1 - 2
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,0.0,0.0)
        mp.CreateNewNode(3,0.0,1.0,0.0)
        mp.CreateNewNode(4,1.0,1.0,0.0)

        mp.CreateNewElement("Element2D3N", 1, [1,4,3], mp.GetProperties()[1])
        mp.CreateNewElement("Element2D3N", 2, [1,2,4], mp.GetProperties()[1])

        #generate particles
        particles_mp.CreateNewNode(1,0.0,0.1,0.0)
        particles_mp.CreateNewNode(2,0.2,0.3,0.0)
        particles_mp.CreateNewNode(3,0.7,0.8,0.0)

        locator = KratosMultiphysics.BinBasedFastPointLocator2D(mp)
        locator.UpdateSearchDatabase()

        #non historical version
        KratosMultiphysics.ParticlesUtilities.CountParticlesInNodesNonHistorical(locator,mp, particles_mp, KratosMultiphysics.AUX_INDEX, 1e-5)

        self.assertEqual(mp.Nodes[1].GetValue(KratosMultiphysics.AUX_INDEX), 3.0)
        self.assertEqual(mp.Nodes[2].GetValue(KratosMultiphysics.AUX_INDEX), 0.0)
        self.assertEqual(mp.Nodes[3].GetValue(KratosMultiphysics.AUX_INDEX), 3.0)
        self.assertEqual(mp.Nodes[4].GetValue(KratosMultiphysics.AUX_INDEX), 3.0)

        #test historical version
        particles_mp.CreateNewNode(4,0.5,0.1,0.0)

        KratosMultiphysics.ParticlesUtilities.CountParticlesInNodesHistorical(locator,mp, particles_mp, KratosMultiphysics.TEMPERATURE, 1e-5)

        self.assertEqual(mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 4.0)
        self.assertEqual(mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 1.0)
        self.assertEqual(mp.Nodes[3].GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 3.0)
        self.assertEqual(mp.Nodes[4].GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 4.0)

    def test_interpolate_particles(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("Main")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        # 3 - 4
        # | / |
        # 1 - 2
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,0.0,0.0)
        mp.CreateNewNode(3,0.0,1.0,0.0)
        mp.CreateNewNode(4,1.0,1.0,0.0)

        mp.CreateNewElement("Element2D3N", 1, [1,4,3], mp.GetProperties()[1])
        mp.CreateNewElement("Element2D3N", 2, [1,2,4], mp.GetProperties()[1])

        #generate particles
        coords = np.array([[-0.01,0.1,0.0],[0.2,0.3,0.0],[0.7,0.8,0.0]])

        locator = KratosMultiphysics.BinBasedFastPointLocator2D(mp)
        locator.UpdateSearchDatabase()

        for node in mp.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0,node.X)
            node.SetValue(KratosMultiphysics.AUX_INDEX,node.Y)

        #historical version
        is_inside, values = KratosMultiphysics.ParticlesUtilities.InterpolateValuesAtCoordinatesHistorical(locator,coords, KratosMultiphysics.TEMPERATURE, 1e-5)

        self.assertEqual(is_inside[0], False)
        self.assertEqual(is_inside[1], True)
        self.assertEqual(is_inside[2], True)
        #note that we cannot check the value for the first coordinate as is_inside[0]=False
        self.assertEqual(values[1], 0.2)
        self.assertEqual(values[2], 0.7)

        #non historical version
        is_inside, values = KratosMultiphysics.ParticlesUtilities.InterpolateValuesAtCoordinatesNonHistorical(locator,coords, KratosMultiphysics.AUX_INDEX, 1e-5)

        self.assertEqual(is_inside[0], False)
        self.assertEqual(is_inside[1], True)
        self.assertEqual(is_inside[2], True)
        #note that we cannot check the value for the first coordinate as is_inside[0]=False
        self.assertEqual(values[1], 0.3)
        self.assertEqual(values[2], 0.8)

    def test_mark_outsider_particles(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("Main")

        particles_mp = current_model.CreateModelPart("Particles")
        particles_mp.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        # 3 - 4
        # | / |
        # 1 - 2
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,0.0,0.0)
        mp.CreateNewNode(3,0.0,1.0,0.0)
        mp.CreateNewNode(4,1.0,1.0,0.0)

        mp.CreateNewElement("Element2D3N", 1, [1,4,3], mp.GetProperties()[1])
        mp.CreateNewElement("Element2D3N", 2, [1,2,4], mp.GetProperties()[1])

        #generate particles
        #generate particles
        particles_mp.CreateNewNode(1,-1.0,0.1,0.0)
        particles_mp.CreateNewNode(2,0.2,0.3,0.0)
        particles_mp.CreateNewNode(3,0.7,0.8,0.0)

        locator = KratosMultiphysics.BinBasedFastPointLocator2D(mp)
        locator.UpdateSearchDatabase()

        for node in particles_mp.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0,1.0)
            node.SetValue(KratosMultiphysics.AUX_INDEX,2.0)

        #historical version
        KratosMultiphysics.ParticlesUtilities.MarkOutsiderParticlesHistorical(locator,particles_mp, KratosMultiphysics.TEMPERATURE, -1.0,1e-5)
        self.assertEqual(particles_mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), -1.0)
        self.assertEqual(particles_mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 1.0)
        self.assertEqual(particles_mp.Nodes[3].GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 1.0)


        #non historical version
        KratosMultiphysics.ParticlesUtilities.MarkOutsiderParticlesNonHistorical(locator,particles_mp, KratosMultiphysics.AUX_INDEX, -2.0,1e-5)
        self.assertEqual(particles_mp.Nodes[1].GetValue(KratosMultiphysics.AUX_INDEX), -2.0)
        self.assertEqual(particles_mp.Nodes[2].GetValue(KratosMultiphysics.AUX_INDEX), 2.0)
        self.assertEqual(particles_mp.Nodes[3].GetValue(KratosMultiphysics.AUX_INDEX), 2.0)

    def test_classify_particles_in_elements(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("Main")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        particles_mp = current_model.CreateModelPart("Particles")

        # 3 - 4
        # | / |
        # 1 - 2
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,0.0,0.0)
        mp.CreateNewNode(3,0.0,1.0,0.0)
        mp.CreateNewNode(4,1.0,1.0,0.0)

        mp.CreateNewElement("Element2D3N", 1, [1,4,3], mp.GetProperties()[1])
        mp.CreateNewElement("Element2D3N", 2, [1,2,4], mp.GetProperties()[1])

        #generate particles
        particles_mp.CreateNewNode(1,0.0,0.1,0.0)
        particles_mp.CreateNewNode(2,0.2,0.3,0.0)
        particles_mp.CreateNewNode(3,0.7,0.6,0.0)
        particles_mp.CreateNewNode(4,0.1,0.15,0.0)
        particles_mp.Nodes[1].SetValue(KratosMultiphysics.AUX_INDEX,0.0)
        particles_mp.Nodes[2].SetValue(KratosMultiphysics.AUX_INDEX,1.0)
        particles_mp.Nodes[3].SetValue(KratosMultiphysics.AUX_INDEX,1.0)
        particles_mp.Nodes[4].SetValue(KratosMultiphysics.AUX_INDEX,2.0)

        locator = KratosMultiphysics.BinBasedFastPointLocator2D(mp)
        locator.UpdateSearchDatabase()

        #non historical version
        KratosMultiphysics.ParticlesUtilities.ClassifyParticlesInElementsNonHistorical(locator,
            mp,
            particles_mp,
            3, #to be changed
            KratosMultiphysics.AUX_INDEX,
            KratosMultiphysics.MARKER_LABELS,
            1e-9
            )

        self.assertEqual(mp.Elements[1].GetValue(KratosMultiphysics.MARKER_LABELS)[0], 1)
        self.assertEqual(mp.Elements[1].GetValue(KratosMultiphysics.MARKER_LABELS)[1], 1)
        self.assertEqual(mp.Elements[1].GetValue(KratosMultiphysics.MARKER_LABELS)[2], 1)
        self.assertEqual(mp.Elements[2].GetValue(KratosMultiphysics.MARKER_LABELS)[0], 0)
        self.assertEqual(mp.Elements[2].GetValue(KratosMultiphysics.MARKER_LABELS)[1], 1)
        self.assertEqual(mp.Elements[2].GetValue(KratosMultiphysics.MARKER_LABELS)[0], 0)

if __name__ == '__main__':
    KratosUnittest.main()
