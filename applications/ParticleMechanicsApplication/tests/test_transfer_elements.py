
import KratosMultiphysics
from KratosMultiphysics import KratosUnittest
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
data_comm = KratosMultiphysics.DataCommunicator.GetDefault()

class TestTransferObjects(KratosUnittest.TestCase):

    def _set_up_model_parts(self, current_model, dimension):
        # Initialize model part
        ## Material model part definition
        self.material_point_model_part = current_model.CreateModelPart("MPMModelPart")
        self.material_point_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        ## Initial material model part definition
        self.initial_mesh_model_part = current_model.CreateModelPart("Initial")
        self.initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        ## Grid model part definition
        self.grid_model_part = current_model.CreateModelPart("Background_Grid")
        self.grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)


    def _generate_elements(self, current_model, dimension, geometry_element, num_particle):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # Create element and nodes for background grids
        if self.grid_model_part.HasSubModelPart("Sub_Background_Grid"):
            sub_background = self.grid_model_part.GetSubModelPart("Sub_Background_Grid")
        else:
            sub_background = self.grid_model_part.CreateSubModelPart("Sub_Background_Grid")

        self._create_nodes(sub_background, dimension, geometry_element)
        self._create_elements(sub_background,dimension, geometry_element)

        # Create element and nodes for initial meshes
        if self.initial_mesh_model_part.HasSubModelPart("Elements"):
            sub_mp = self.initial_mesh_model_part.GetSubModelPart("Elements")
        else:
            sub_mp = self.initial_mesh_model_part.CreateSubModelPart("Elements")
            sub_mp.GetProperties()[1].SetValue(KratosParticle.PARTICLES_PER_ELEMENT, num_particle)
            sub_mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY, 1000.0)

        self._create_nodes(sub_mp, dimension, geometry_element)
        self._create_elements(sub_mp,dimension, geometry_element)

        # Generate MP Elements
        KratosParticle.GenerateMaterialPointElement(self.grid_model_part, self.initial_mesh_model_part, self.material_point_model_part, False)

    def _create_nodes(self, initial_mp, dimension, geometry_element):
        if geometry_element == "Triangle":
            initial_mp.CreateNewNode(1, 0.0, 0.0, 0.0)
            initial_mp.CreateNewNode(2, 1.0, 0.0, 0.0)
            initial_mp.CreateNewNode(3, 0.0, 1.0, 0.0)
            if (dimension == 3):
                initial_mp.CreateNewNode(4, 0.0, 0.0, 1.0)
        elif geometry_element == "Quadrilateral":
            initial_mp.CreateNewNode(1, -0.5, -0.5, 0.0)
            initial_mp.CreateNewNode(2,  0.5, -0.5, 0.0)
            initial_mp.CreateNewNode(3,  0.5,  0.5, 0.0)
            initial_mp.CreateNewNode(4, -0.5,  0.5, 0.0)
            if (dimension == 3):
                initial_mp.CreateNewNode(5, -0.5, -0.5, 1.0)
                initial_mp.CreateNewNode(6,  0.5, -0.5, 1.0)
                initial_mp.CreateNewNode(7,  0.5,  0.5, 1.0)
                initial_mp.CreateNewNode(8, -0.5,  0.5, 1.0)

    def _create_elements(self, initial_mp, dimension, geometry_element):
        if geometry_element == "Triangle":
            if (dimension == 2):
                initial_mp.CreateNewElement("Element2D3N", 1, [1,2,3], initial_mp.GetProperties()[1])
            if (dimension == 3):
                initial_mp.CreateNewElement("Element3D4N", 1, [1,2,3,4], initial_mp.GetProperties()[1])
        elif geometry_element == "Quadrilateral":
            if (dimension == 2):
                initial_mp.CreateNewElement("Element2D4N", 1, [1,2,3,4], initial_mp.GetProperties()[1])
            if (dimension == 3):
                initial_mp.CreateNewElement("Element3D8N", 1, [1,2,3,4,5,6,7,8], initial_mp.GetProperties()[1])

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mp.Elements)

    def _check_elements(self, mp, dimension):
        SomeProcessInfo = KratosMultiphysics.ProcessInfo()
        for el in mp.Elements:
            # Material Id
            material_id = el.CalculateOnIntegrationPoints(KratosParticle.MP_MATERIAL_ID, SomeProcessInfo)
            self.assertEqual(material_id[0], 1)
            # Check geometry
            self.assertEqual( el.GetGeometry().WorkingSpaceDimension(), dimension)
            self.assertEqual( el.GetGeometry().LocalSpaceDimension(), dimension)
            self.assertEqual( el.GetGeometry().IntegrationPointsNumber(), 1)
            jacobian = el.GetGeometry().Jacobian(0)
            shape_functions_values = el.GetGeometry().ShapeFunctionsValues()
            shape_functions_derivatives = el.GetGeometry().ShapeFunctionDerivatives(1,0)
            if dimension is 2:
                jacobian_ref = KratosMultiphysics.Matrix(2,2,0.0)
                jacobian_ref[0,0] = 1.0
                jacobian_ref[1,1] = 1.0
                shape_functions_values_ref = KratosMultiphysics.Matrix(1,3,1.0/3.0)
                shape_functions_derivatives_ref = KratosMultiphysics.Matrix(3,2,0.0)
                shape_functions_derivatives_ref[0,0] = -1.0
                shape_functions_derivatives_ref[0,1] = -1.0
                shape_functions_derivatives_ref[1,0] = 1.0
                shape_functions_derivatives_ref[2,1] = 1.0

                self.assertMatrixAlmostEqual( jacobian,jacobian_ref,7)
                self.assertMatrixAlmostEqual( shape_functions_values,shape_functions_values_ref,7)
                self.assertMatrixAlmostEqual( shape_functions_derivatives, shape_functions_derivatives_ref, 7)
                self.assertVectorAlmostEqual( el.GetGeometry().Center(), [1.0/3.0, 1.0/3.0, 0.0])
            else:
                jacobian_ref = KratosMultiphysics.Matrix(3,3)
                jacobian_ref[0,0] = 1.0
                jacobian_ref[1,1] = 1.0
                jacobian_ref[2,2] = 1.0
                shape_functions_values_ref = KratosMultiphysics.Matrix(1,4,0.25)
                shape_functions_derivatives_ref = KratosMultiphysics.Matrix(4,3,0.0)
                shape_functions_derivatives_ref[0,0] = -1.0
                shape_functions_derivatives_ref[0,1] = -1.0
                shape_functions_derivatives_ref[0,2] = -1.0
                shape_functions_derivatives_ref[1,0] = 1.0
                shape_functions_derivatives_ref[2,1] = 1.0
                shape_functions_derivatives_ref[3,2] = 1.0
                self.assertMatrixAlmostEqual( jacobian,jacobian_ref,7)
                self.assertMatrixAlmostEqual( shape_functions_values,shape_functions_values_ref,7)
                self.assertMatrixAlmostEqual( shape_functions_derivatives, shape_functions_derivatives_ref, 7)
                self.assertVectorAlmostEqual( el.GetGeometry().Center(), [0.25, 0.25, 0.25])
            #add the rest fro add_geometreis_to_python.cpp
            # Material Point variables
            xg = el.CalculateOnIntegrationPoints(KratosParticle.MPC_COORD, SomeProcessInfo)
            mass = el.CalculateOnIntegrationPoints(KratosParticle.MP_MASS, SomeProcessInfo)
            volume = el.CalculateOnIntegrationPoints(KratosParticle.MP_VOLUME, SomeProcessInfo)
            if dimension is 2:
                self.assertVectorAlmostEqual( xg[0], [1.0/3.0, 1.0/3.0, 0.0], 7)
                self.assertAlmostEqual( mass[0], 500.0, 7)
                self.assertAlmostEqual( volume[0], 0.5, 7)
            else:
                self.assertVectorAlmostEqual( xg[0], [0.25, 0.25, 0.25], 7)
                self.assertAlmostEqual( mass[0], 166.6666666666666666, 7)
                self.assertAlmostEqual( volume[0], 0.1666666666666666, 7)

            density = el.CalculateOnIntegrationPoints(KratosParticle.MP_DENSITY, SomeProcessInfo)
            self.assertAlmostEqual( density[0], 1000.0, 7)
            displacement = el.CalculateOnIntegrationPoints(KratosParticle.MP_DISPLACEMENT, SomeProcessInfo)
            self.assertVectorAlmostEqual( displacement[0], [0.1, 2.21, 3.0], 7)
            velocity = el.CalculateOnIntegrationPoints(KratosParticle.MP_VELOCITY, SomeProcessInfo)
            self.assertVectorAlmostEqual( velocity[0], [0.5, 2.25, 3.5], 7)
            acceleration = el.CalculateOnIntegrationPoints(KratosParticle.MP_ACCELERATION, SomeProcessInfo)
            self.assertVectorAlmostEqual( acceleration[0], [0.2, 2.22, 3.2], 7)
            volume_acceleration = el.CalculateOnIntegrationPoints(KratosParticle.MP_VOLUME_ACCELERATION, SomeProcessInfo)
            self.assertVectorAlmostEqual( volume_acceleration[0], [0.3, 2.32, 1.2], 7)
            chauchy_stress_vector = el.CalculateOnIntegrationPoints(KratosParticle.MP_CAUCHY_STRESS_VECTOR, SomeProcessInfo)
            self.assertVectorAlmostEqual( chauchy_stress_vector[0], [1.2,2.0,3.45], 7)
            almansi_strain_vector = el.CalculateOnIntegrationPoints(KratosParticle.MP_ALMANSI_STRAIN_VECTOR, SomeProcessInfo)
            self.assertVectorAlmostEqual( almansi_strain_vector[0], [1.6,2.0,1.45], 7)



    def _transfer_elements(self, dimension, geometry_element):
        ''' One element is created in rank=0 and send to all other processes.
            The test is passed if the recieving processes hold the element with the correct properties.'''
        current_model = KratosMultiphysics.Model()
        rank = data_comm.Rank()
        size = data_comm.Size()
        send_elements = []

        for i in range(size):
            send_elements.append( KratosMultiphysics.ElementsArray() )

        self._set_up_model_parts(current_model,dimension)

        if( rank == 0): # Sender
            self._generate_elements(current_model, dimension, geometry_element, num_particle=1)
            mp = current_model.GetModelPart("MPMModelPart")
            SomeProcessInfo = KratosMultiphysics.ProcessInfo()
            for i in range(size):
                if i is not rank:
                    for element in mp.Elements:
                        send_elements[i].append(element)
                        for el in send_elements[i]:
                            #Give elements some pseudo variables
                            el.SetValuesOnIntegrationPoints(KratosParticle.MP_DISPLACEMENT, [[0.1, 2.21, 3.0]], SomeProcessInfo)
                            el.SetValuesOnIntegrationPoints(KratosParticle.MP_VELOCITY, [[0.5, 2.25, 3.5]], SomeProcessInfo)
                            el.SetValuesOnIntegrationPoints(KratosParticle.MP_ACCELERATION, [[0.2, 2.22, 3.2]], SomeProcessInfo)
                            el.SetValuesOnIntegrationPoints(KratosParticle.MP_VOLUME_ACCELERATION, [[0.3, 2.32, 1.2]], SomeProcessInfo)
                            cauchy_stress_vector = [KratosMultiphysics.Vector([1.2,2.0,3.45])]
                            el.SetValuesOnIntegrationPoints(KratosParticle.MP_CAUCHY_STRESS_VECTOR, cauchy_stress_vector, 0, SomeProcessInfo)
                            almansi_strain_vector = [KratosMultiphysics.Vector([1.6,2.0,1.45])]
                            el.SetValuesOnIntegrationPoints(KratosParticle.MP_ALMANSI_STRAIN_VECTOR, almansi_strain_vector, 0, SomeProcessInfo)
        else: # Recievers
            #Make sure all ModelParts have same SubmodelParts
            mp = current_model.GetModelPart("MPMModelPart")
            mp.CreateSubModelPart("Elements")

        sub_mp = mp.GetSubModelPart("Elements")
        KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(sub_mp)
        #Send elements from 0 to all other
        KratosParticle.TransferElements(sub_mp, send_elements)
        #Check if model part hold the correct elements
        if rank is 0:
            self.assertEqual(mp.NumberOfElements(), 0) #Check if element was removed after sent
        else:
            self.assertEqual(mp.NumberOfElements(), 1) #Check if element was added
            self._check_elements(mp,dimension) #Check properties of element

    def test_transfer_elements_triangle2D3N(self):
        self._transfer_elements(dimension=2,geometry_element="Triangle")

    def test_transfer_elements_triangle3D4N(self):
        self._transfer_elements(dimension=3,geometry_element="Triangle")







if __name__ == '__main__':
    KratosUnittest.main()