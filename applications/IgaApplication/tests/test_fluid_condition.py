import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.FluidDynamicsApplication as DFA
import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_creation_utility import TestCreationUtility

class TestSupportFluidCondition(KratosUnittest.TestCase):

    def create_condition(self, model_part, degree, integration_point):
        # Define properties
        props = model_part.CreateNewProperties(0)
        props.SetValue(IGA.PENALTY_FACTOR, 100.0)

        # Assign constitutive law
        law = DFA.Newtonian2DLaw()
        props.SetValue(KM.CONSTITUTIVE_LAW, law)

        # Create quadrature point geometry on a curve
        geometry = TestCreationUtility.GetQuadraturePointGeometryOnCurve(model_part, degree, integration_point)

        # Create the condition and add it to the model part
        condition = model_part.CreateNewCondition("SupportFluidCondition", 1, geometry, props)
        model_part.AddCondition(condition)

        return condition

    def test_SupportFluidConditionP3(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        # Add required variables
        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        process_info = model_part.ProcessInfo

        # Define integration point in parametric space
        ipt = [0.333333333333333, 0.05, 0.0, 0.086963711284364]

        # Create the support condition
        condition = self.create_condition(model_part, 3, ipt)

        # Set prescribed displacement value
        u_D = KM.Vector(3)
        u_D[0] = 0.1
        u_D[1] = -0.5
        u_D[2] = 0.0
        condition.SetValue(KM.VELOCITY, u_D)

        # Add DOFs
        for node in model_part.Nodes:
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.PRESSURE)

        # Initialize and step the condition
        condition.Initialize(process_info)
        
        # Compute local system
        lhs = KM.Matrix()
        rhs = KM.Vector()
        condition.CalculateLocalSystem(lhs, rhs, process_info)

        # expected values
        expected_LHS = [0.9082884211890723, 0.1914419633305822, 0.1132391371476317, -0.7120647472279207, -0.1114426344123852, 0.0253471024528432, -0.1214269163593457, -0.0198272566535484, 0.0018912059880185, -0.0047186079794445, -0.0007794380537414, 0.0000470357521900, 0.1717725084961046, -0.0093879512287520, 0.0303423353503050, -0.2068244902336144, -0.0434444351586120, 0.0067917356327083, -0.0337320793440627, -0.0063261916428638, 0.0005067471172104, -0.0012940885407889, -0.0002340561806795, 0.0000126031918147]
        expected_RHS = [1.6580668225671398,0.0000000000000000,-0.0147966462120511,0.3711366112802327,0.0000000000000000,0.0125886221611062,0.0276913616825390,0.0000000000000000,0.0021256511581339,0.0006887055318957,0.0000000000000000,0.0000823728928110,0.4442776661037025,0.0000000000000000,-0.0039647494032081,0.0994457552741618,0.0000000000000000,0.0033731111418889,0.0074198780001544,0.0000000000000000,0.0005695665112123,0.0001845380910943,0.0000000000000000,0.0000220717501069]

        tolerance = 1e-7

        for i in range(rhs.Size()):
                self.assertAlmostEqual(lhs[0, i], expected_LHS[i], delta=tolerance)

        # Confronto RHS
        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

        
if __name__ == '__main__':
    KratosUnittest.main()