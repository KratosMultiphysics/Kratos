import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtils

class FluidAuxiliaryUtilitiesTest(UnitTest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        fluid_model_part = self.model.CreateModelPart("FluidModelPart")
        fluid_model_part.SetBufferSize(2)
        fluid_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 2)
        fluid_model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        fluid_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        Kratos.ModelPartIO("Cavity/square5").ReadModelPart(fluid_model_part)

    def testCalculateFluidVolume(self):
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        fluid_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidVolume(fluid_model_part)
        self.assertAlmostEqual(fluid_volume, 1.0, 12)

    def testCalculateFluidPositiveVolume(self):
        # Set fluid level set
        level_set_y = 1.0/3.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Calculate the fluid positive volume
        fluid_positive_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidPositiveVolume(fluid_model_part)
        self.assertAlmostEqual(fluid_positive_volume, 1.0 - level_set_y, 12)

    def testCalculateFluidNegativeVolume(self):
        # Set fluid level set
        level_set_y = 2.0/3.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Calculate the fluid negative volume
        fluid_negative_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(fluid_model_part)
        self.assertAlmostEqual(fluid_negative_volume, level_set_y, 12)

    def testCalculateFlowRatePositiveSkin(self):
        # Set fluid level set
        level_set_y = 2.0/3.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY, 0, [1.0,0.0,0.0])
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Call the tetrahedral mesh orientation process to calculate the normals and neighbours
        tmoc = Kratos.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        Kratos.TetrahedralMeshOrientationCheck(fluid_model_part, throw_errors, flags).Execute()

        # Calculate the left wall flow
        ref_positive_flow_rate_left = -1.0/3.0
        left_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_left_wall")
        for cond in left_skin_model_part.Conditions:
            cond.Set(Kratos.INLET, True)
        positive_flow_rate_left = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRatePositiveSkin(left_skin_model_part, Kratos.INLET)
        self.assertAlmostEqual(positive_flow_rate_left, ref_positive_flow_rate_left, 12)

        # Calculate the right wall flow
        ref_positive_flow_rate_right = 1.0/3.0
        right_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_right_wall")
        for cond in right_skin_model_part.Conditions:
            cond.Set(Kratos.OUTLET, True)
        positive_flow_rate_right = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRatePositiveSkin(right_skin_model_part, Kratos.OUTLET)
        self.assertAlmostEqual(positive_flow_rate_right, ref_positive_flow_rate_right, 12)

        # Calculate the top wall flow (without flag)
        ref_positive_flow_rate_top = 0.0
        top_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_top_wall")
        positive_flow_rate_top = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRatePositiveSkin(top_skin_model_part)
        self.assertAlmostEqual(positive_flow_rate_top, ref_positive_flow_rate_top, 12)

    def testCalculateFlowRateNegativeSkin(self):
        # Set fluid level set
        level_set_y = 2.0/3.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY, 0, [1.0,0.0,0.0])
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Call the tetrahedral mesh orientation process to calculate the normals and neighbours
        tmoc = Kratos.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        Kratos.TetrahedralMeshOrientationCheck(fluid_model_part, throw_errors, flags).Execute()

        # Calculate the left wall flow
        ref_negative_flow_rate_left = -2.0/3.0
        left_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_left_wall")
        for cond in left_skin_model_part.Conditions:
            cond.Set(Kratos.INLET, True)
        negative_flow_rate_left = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(left_skin_model_part, Kratos.INLET)
        self.assertAlmostEqual(negative_flow_rate_left, ref_negative_flow_rate_left, 12)

        # Calculate the right wall flow
        ref_negative_flow_rate_right = 2.0/3.0
        right_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_right_wall")
        for cond in right_skin_model_part.Conditions:
            cond.Set(Kratos.OUTLET, True)
        negative_flow_rate_right = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(right_skin_model_part, Kratos.OUTLET)
        self.assertAlmostEqual(negative_flow_rate_right, ref_negative_flow_rate_right, 12)

        # Calculate the top wall flow (without flag)
        ref_negative_flow_rate_top = 0.0
        top_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_top_wall")
        negative_flow_rate_top = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(top_skin_model_part)
        self.assertAlmostEqual(negative_flow_rate_top, ref_negative_flow_rate_top, 12)

    def testMapVelocityFromSkinToVolumeRBF(self):
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.GetValue(Kratos.EMBEDDED_VELOCITY)
            node.SetSolutionStepValue(Kratos.DISTANCE,0,node.X-0.5)

        skin_model_part = self.model.CreateModelPart("skin_part")
        skin_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        skin_model_part.CreateNewNode(1,0.0,0.0,0.0)
        skin_model_part.CreateNewNode(2,0.1,0.0,0.0)
        skin_model_part.CreateNewNode(3,0.1,0.1,0.0)
        ref_v = Kratos.Vector([1.0,2.0,3.0])
        for node in skin_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY,0,ref_v)

        KratosFluid.FluidAuxiliaryUtilities.MapVelocityFromSkinToVolumeRBF(fluid_model_part, skin_model_part, 1.2)

        for elem in fluid_model_part.Elements:
            npos=0
            nneg=0
            for node in elem.GetGeometry():
                if(node.GetSolutionStepValue(Kratos.DISTANCE)>=0):
                    npos+=1
                else:
                    nneg+=1

            if(npos>0 and nneg>0):
                for node in elem.GetGeometry():
                    v = node.GetValue(Kratos.EMBEDDED_VELOCITY)
                    self.assertVectorAlmostEqual(v,ref_v)

    def testFindMaximumEdgeLength(self):
        calculate_nodal_neighbours = True
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        max_edge_length = KratosFluid.FluidAuxiliaryUtilities.FindMaximumEdgeLength(fluid_model_part, calculate_nodal_neighbours)
        ref_max_edge_length = 0.2828427124746191
        self.assertAlmostEqual(max_edge_length, ref_max_edge_length, 12)
    
    def testCalculateFluidNegativeCutVolume(self):
        # Set fluid level set
        level_set_y = 1.0/2.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Calculate the fluid negative volume
        fluid_negative_cut_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidCutElementsNegativeVolume(fluid_model_part)
        self.assertAlmostEqual(fluid_negative_cut_volume,0.1, 12)
    def testCalculateFluidPositiveCutVolume(self):
        # Set fluid level set
        level_set_y = 1.0/2.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)
        # Calculate the fluid negative volume
        fluid_positive_cut_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidCutElementPositiveVolume(fluid_model_part)
        print(fluid_positive_cut_volume)
        self.assertAlmostEqual(fluid_positive_cut_volume, 0.1, 12)

    def tearDown(self):
        KratosUtils.DeleteFileIfExisting("Cavity/square5.time")


if __name__ == '__main__':
    UnitTest.main()
