
import KratosMultiphysics
from KratosMultiphysics import KratosUnittest
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MPMApplication as KratosMPM
data_comm = KratosMultiphysics.DataCommunicator.GetDefault()

class TestTransferElements(KratosUnittest.TestCase):
    ''' This class provides all required methods to test the MPM_MPI_Utilities::TransferElements function.
        Different tests for all available elements are performed.
        New developed elements must be added here.
    '''
    def _set_up_model_parts(self, current_model, dimension, is_pqmpm):
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

        if is_pqmpm:
            self.grid_model_part.ProcessInfo[KratosMPM.IS_PQMPM] = True

    def _generate_particle_elements(self, current_model, dimension, geometry_element, num_particle, is_mixed_formulation):
        # Create element and nodes for background grids
        if self.grid_model_part.HasSubModelPart("Sub_Background_Grid"):
            sub_background = self.grid_model_part.GetSubModelPart("Sub_Background_Grid")
        else:
            sub_background = self.grid_model_part.CreateSubModelPart("Sub_Background_Grid")

        self._create_nodes(sub_background, dimension, geometry_element)
        self._create_elements(sub_background,dimension, geometry_element)

        # Create element and nodes for initial meshes
        if self.initial_mesh_model_part.HasSubModelPart("Elements"):
            sub_initial = self.initial_mesh_model_part.GetSubModelPart("Elements")
        else:
            sub_initial = self.initial_mesh_model_part.CreateSubModelPart("Elements")
            sub_initial.GetProperties()[1].SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, num_particle)
            sub_initial.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY, 1000.0)

        self._create_nodes(sub_initial, dimension, geometry_element)
        self._create_elements(sub_initial,dimension, geometry_element)

        # Generate MP Elements
        KratosMPM.GenerateMaterialPointElement(self.grid_model_part, self.initial_mesh_model_part, self.material_point_model_part, is_mixed_formulation)

    def _create_nodes(self, mp, dimension, geometry_element):
        if geometry_element == "Triangle":
            mp.CreateNewNode(1, 0.0, 0.0, 0.0)
            mp.CreateNewNode(2, 1.0, 0.0, 0.0)
            mp.CreateNewNode(3, 0.0, 1.0, 0.0)
            if (dimension == 3):
                mp.CreateNewNode(4, 0.0, 0.0, 1.0)
        elif geometry_element == "Quadrilateral":
            mp.CreateNewNode(1, -0.5, -0.5, 0.0)
            mp.CreateNewNode(2,  0.5, -0.5, 0.0)
            mp.CreateNewNode(3,  0.5,  0.5, 0.0)
            mp.CreateNewNode(4, -0.5,  0.5, 0.0)
            if (dimension == 3):
                mp.CreateNewNode(5, -0.5, -0.5, 1.0)
                mp.CreateNewNode(6,  0.5, -0.5, 1.0)
                mp.CreateNewNode(7,  0.5,  0.5, 1.0)
                mp.CreateNewNode(8, -0.5,  0.5, 1.0)

    def _create_elements(self, mp, dimension, geometry_element):
        if geometry_element == "Triangle":
            if (dimension == 2):
                mp.CreateNewElement("Element2D3N", 1, [1,2,3], mp.GetProperties()[1])
            if (dimension == 3):
                mp.CreateNewElement("Element3D4N", 1, [1,2,3,4], mp.GetProperties()[1])
        elif geometry_element == "Quadrilateral":
            if (dimension == 2):
                mp.CreateNewElement("Element2D4N", 1, [1,2,3,4], mp.GetProperties()[1])
            if (dimension == 3):
                mp.CreateNewElement("Element3D8N", 1, [1,2,3,4,5,6,7,8], mp.GetProperties()[1])

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, mp.Elements)

    def _check_elements(self, mp, dimension, element_type, is_mixed_formulation):
        process_info = KratosMultiphysics.ProcessInfo()
        for el in mp.Elements:
            # Check material Id
            material_id = el.CalculateOnIntegrationPoints(KratosMPM.MP_MATERIAL_ID, process_info)
            self.assertEqual(material_id[0], 1)
            # Check geometry
            self.assertEqual( el.GetGeometry().WorkingSpaceDimension(), dimension)
            self.assertEqual( el.GetGeometry().LocalSpaceDimension(), dimension)
            self.assertEqual( el.GetGeometry().IntegrationPointsNumber(), 1)
            jacobian = el.GetGeometry().Jacobian(0)
            shape_functions_values = el.GetGeometry().ShapeFunctionsValues()
            shape_functions_derivatives = el.GetGeometry().ShapeFunctionDerivatives(1,0)
            center = el.GetGeometry().Center()
            if dimension == 2:
                jacobian_ref = KratosMultiphysics.Matrix(2,2,0.0)
                if element_type == "Triangle":
                    jacobian_ref[0,0] = 1.0
                    jacobian_ref[1,1] = 1.0
                    shape_functions_values_ref = KratosMultiphysics.Matrix(1,3,1.0/3.0)
                    shape_functions_derivatives_ref = KratosMultiphysics.Matrix(3,2,0.0)
                    shape_functions_derivatives_ref[0,0] = -1.0
                    shape_functions_derivatives_ref[0,1] = -1.0
                    shape_functions_derivatives_ref[1,0] = 1.0
                    shape_functions_derivatives_ref[2,1] = 1.0
                    center_ref = [1.0/3.0, 1.0/3.0, 0.0]
                else:
                    jacobian_ref[0,0] = 0.5
                    jacobian_ref[1,1] = 0.5
                    shape_functions_values_ref = KratosMultiphysics.Matrix(1,4,0.25)
                    shape_functions_derivatives_ref = KratosMultiphysics.Matrix(4,2,0.25)
                    shape_functions_derivatives_ref[0,0] = -0.25
                    shape_functions_derivatives_ref[0,1] = -0.25
                    shape_functions_derivatives_ref[1,1] = -0.25
                    shape_functions_derivatives_ref[3,0] = -0.25
                    center_ref = [0.0, 0.0, 0.0]
            else:
                jacobian_ref = KratosMultiphysics.Matrix(3,3,0.0)
                if element_type == "Triangle":
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
                    center_ref = [0.25, 0.25, 0.25]
                else:
                    jacobian_ref[0,0] = 0.5
                    jacobian_ref[1,1] = 0.5
                    jacobian_ref[2,2] = 0.5
                    shape_functions_values_ref = KratosMultiphysics.Matrix(1,8,0.125)
                    shape_functions_derivatives_ref = KratosMultiphysics.Matrix(8,3,0.125)
                    helper_vector = [-1,-1,-1,1,-1,-1,1,1,-1,-1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,1,1]
                    count = 0
                    for i in range(8):
                        for j in range(3):
                            shape_functions_derivatives_ref[i,j] *= helper_vector[count]
                            count = count + 1
                    center_ref = [0.0, 0.0, 0.5]

            self.assertMatrixAlmostEqual( jacobian,jacobian_ref,7)
            self.assertMatrixAlmostEqual( shape_functions_values,shape_functions_values_ref,7)
            self.assertMatrixAlmostEqual( shape_functions_derivatives, shape_functions_derivatives_ref, 7)
            self.assertVectorAlmostEqual( center, center_ref)
            # Material Point variables
            xg = el.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD, process_info)
            mass = el.CalculateOnIntegrationPoints(KratosMPM.MP_MASS, process_info)
            volume = el.CalculateOnIntegrationPoints(KratosMPM.MP_VOLUME, process_info)
            if dimension == 2:
                if element_type == "Triangle":
                    self.assertVectorAlmostEqual( xg[0], [1.0/3.0, 1.0/3.0, 0.0], 7)
                    self.assertAlmostEqual( mass[0], 500.0, 7)
                    self.assertAlmostEqual( volume[0], 0.5, 7)
                else:
                    self.assertVectorAlmostEqual( xg[0], [0.0, 0.0, 0.0], 7)
                    self.assertAlmostEqual( mass[0], 1000.0, 7)
                    self.assertAlmostEqual( volume[0], 1.0, 7)
            else:
                if element_type == "Triangle":
                    self.assertVectorAlmostEqual( xg[0], [0.25, 0.25, 0.25], 7)
                    self.assertAlmostEqual( mass[0], 166.6666666666666666, 7)
                    self.assertAlmostEqual( volume[0], 0.1666666666666666, 7)
                else:
                    self.assertVectorAlmostEqual( xg[0], [0.0, 0.0, 0.5], 7)
                    self.assertAlmostEqual( mass[0], 1000.0, 7)
                    self.assertAlmostEqual( volume[0], 1.0, 7)

            density = el.CalculateOnIntegrationPoints(KratosMPM.MP_DENSITY, process_info)
            self.assertAlmostEqual( density[0], 1000.0, 7)
            displacement = el.CalculateOnIntegrationPoints(KratosMPM.MP_DISPLACEMENT, process_info)
            self.assertVectorAlmostEqual( displacement[0], [0.1, 2.21, 3.0], 7)
            velocity = el.CalculateOnIntegrationPoints(KratosMPM.MP_VELOCITY, process_info)
            self.assertVectorAlmostEqual( velocity[0], [0.5, 2.25, 3.5], 7)
            acceleration = el.CalculateOnIntegrationPoints(KratosMPM.MP_ACCELERATION, process_info)
            self.assertVectorAlmostEqual( acceleration[0], [0.2, 2.22, 3.2], 7)
            volume_acceleration = el.CalculateOnIntegrationPoints(KratosMPM.MP_VOLUME_ACCELERATION, process_info)
            self.assertVectorAlmostEqual( volume_acceleration[0], [0.3, 2.32, 1.2], 7)
            chauchy_stress_vector = el.CalculateOnIntegrationPoints(KratosMPM.MP_CAUCHY_STRESS_VECTOR, process_info)
            self.assertVectorAlmostEqual( chauchy_stress_vector[0], [1.2,2.0,3.45], 7)
            almansi_strain_vector = el.CalculateOnIntegrationPoints(KratosMPM.MP_ALMANSI_STRAIN_VECTOR, process_info)
            self.assertVectorAlmostEqual( almansi_strain_vector[0], [1.6,2.0,1.45], 7)
            if is_mixed_formulation:
                # updated_langrangian_up members
                pressure = el.CalculateOnIntegrationPoints(KratosMPM.MP_PRESSURE, process_info)
                self.assertAlmostEqual(pressure[0], 3.3)

    def _transfer_elements(self, dimension, geometry_element, is_mixed_formulation, is_pqmpm):
        ''' One element is created in rank=0 and send to all other processes.
            The test is passed if the recieving processes hold the element with the correct properties.'''
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        rank = data_comm.Rank()
        size = data_comm.Size()

        current_model = KratosMultiphysics.Model()
        send_elements = []

        for i in range(size):
            send_elements.append( KratosMultiphysics.ElementsArray() )

        self._set_up_model_parts(current_model,dimension,is_pqmpm)
        num_particle = 1
        if( rank == 0): # Sender
            self._generate_particle_elements(current_model, dimension, geometry_element, num_particle, is_mixed_formulation)
            mp = current_model.GetModelPart("MPMModelPart")
            process_info = KratosMultiphysics.ProcessInfo()
            for i in range(size):
                if i != rank:
                    for element in mp.Elements:
                        send_elements[i].append(element)
                        for el in send_elements[i]:
                            #Give elements some pseudo variables
                            el.SetValuesOnIntegrationPoints(KratosMPM.MP_DISPLACEMENT, [[0.1, 2.21, 3.0]], process_info)
                            el.SetValuesOnIntegrationPoints(KratosMPM.MP_VELOCITY, [[0.5, 2.25, 3.5]], process_info)
                            el.SetValuesOnIntegrationPoints(KratosMPM.MP_ACCELERATION, [[0.2, 2.22, 3.2]], process_info)
                            el.SetValuesOnIntegrationPoints(KratosMPM.MP_VOLUME_ACCELERATION, [[0.3, 2.32, 1.2]], process_info)
                            cauchy_stress_vector = [KratosMultiphysics.Vector([1.2,2.0,3.45])]
                            el.SetValuesOnIntegrationPoints(KratosMPM.MP_CAUCHY_STRESS_VECTOR, cauchy_stress_vector, 0, process_info)
                            almansi_strain_vector = [KratosMultiphysics.Vector([1.6,2.0,1.45])]
                            el.SetValuesOnIntegrationPoints(KratosMPM.MP_ALMANSI_STRAIN_VECTOR, almansi_strain_vector, 0, process_info)
                            if is_mixed_formulation:
                                el.SetValuesOnIntegrationPoints(KratosMPM.MP_PRESSURE, [3.3], process_info)
        else: # Recievers
            #Make sure all ModelParts have same SubmodelParts
            mp = current_model.GetModelPart("MPMModelPart")
            mp.CreateSubModelPart("Elements")

        sub_mp = mp.GetSubModelPart("Elements")
        KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(sub_mp)
        #Send elements from rank=0 to all other
        KratosMPM.MPM_MPI_Utilities.TransferElements(sub_mp, send_elements)
        #Check if model_parts hold the correct elements
        if rank == 0:
            self.assertEqual(mp.NumberOfElements(), 0) #Check if element was removed after sent
        else:
            self.assertEqual(mp.NumberOfElements(), 1) #Check if element was added
            self._check_elements(mp,dimension,geometry_element,is_mixed_formulation) #Check properties of element

    ## 2-Dimensions
    # Triangle
    def test_transfer_elements_triangle2D3N(self):
        #updated_langrangian
        self._transfer_elements(dimension=2,geometry_element="Triangle", is_mixed_formulation=False, is_pqmpm = False)

    def test_transfer_elements_triangle2D3N_up(self):
        #updated_langrangian_up
        self._transfer_elements(dimension=2,geometry_element="Triangle", is_mixed_formulation=True, is_pqmpm = False)

    def test_transfer_elements_triangle2D3N_pq(self):
        #updated_langrangian_pq
        self._transfer_elements(dimension=2,geometry_element="Triangle", is_mixed_formulation=False, is_pqmpm = True)

    # Quadrilateral
    def test_transfer_elements_quadrilateral2D4N(self):
        #updated_langrangian
        self._transfer_elements(dimension=2,geometry_element="Quadrilateral", is_mixed_formulation=False, is_pqmpm = False)

    def test_transfer_elements_quadrilateral2D4N_pq(self):
        #updated_langrangian_pq
        self._transfer_elements(dimension=2,geometry_element="Quadrilateral", is_mixed_formulation=False, is_pqmpm = True)

    ## 3-Dimensions
    # Triangle
    def test_transfer_elements_triangle3D4N(self):
        #updated_langrangian
        self._transfer_elements(dimension=3,geometry_element="Triangle", is_mixed_formulation=False, is_pqmpm = False)

    def test_transfer_elements_triangle3D4N_pq(self):
        #updated_langrangian_pq
        self._transfer_elements(dimension=3,geometry_element="Triangle", is_mixed_formulation=False, is_pqmpm = True)

    # Quadrilateral
    def test_transfer_elements_quadrilateral3D8N(self):
        #updated_langrangian
        self._transfer_elements(dimension=3,geometry_element="Quadrilateral", is_mixed_formulation=False, is_pqmpm = False)

    def test_transfer_elements_quadrilateral3D8N_pq(self):
        #updated_langrangian_pq
        self._transfer_elements(dimension=3,geometry_element="Quadrilateral", is_mixed_formulation=False, is_pqmpm = True)

if __name__ == '__main__':
    KratosUnittest.main()