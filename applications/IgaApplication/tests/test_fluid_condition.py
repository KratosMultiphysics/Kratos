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
        geometry = TestCreationUtility.GetQuadraturePointGeometryOnCurveP2(model_part, degree, integration_point)

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
        mesh_size = [0.1, 0.1]
        condition.SetValue(IGA.KNOT_SPAN_SIZES, mesh_size)

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
        expected_LHS = [1.1193281795841600e+02, 0.0000000000000000e+00, 2.2329836010140803e-18, 1.1193281795841582e+02, 0.0000000000000000e+00, 2.2329836010140768e-18, 2.7983204489603917e+01, 0.0000000000000000e+00, 5.5824590025351834e-19, 1.1782401890359578e+01, 0.0000000000000000e+00, 2.3505090536990317e-19, 1.1782401890359560e+01, 0.0000000000000000e+00, 2.3505090536990278e-19, 2.9456004725898861e+00, 0.0000000000000000e+00, 5.8762726342475611e-20, 3.1006320764104151e-01, 0.0000000000000000e+00, 6.1855501413132414e-21, 3.1006320764104101e-01, 0.0000000000000000e+00, 6.1855501413132309e-21, 7.7515801910260154e-02, 0.0000000000000000e+00, 1.5463875353283055e-21]
        expected_RHS = [2.7905688687693718e+01,-1.3952844343846857e+02,0.0000000000000000e+00,2.7905688687693676e+01,-1.3952844343846837e+02,0.0000000000000000e+00,6.9764221719234083e+00,-3.4882110859617036e+01,0.0000000000000000e+00,2.9374409144940756e+00,-1.4687204572470376e+01,0.0000000000000000e+00,2.9374409144940712e+00,-1.4687204572470355e+01,0.0000000000000000e+00,7.3436022862351658e-01,-3.6718011431175830e+00,0.0000000000000000e+00,7.7301076697212506e-02,-3.8650538348606250e-01,0.0000000000000000e+00,7.7301076697212381e-02,-3.8650538348606189e-01,0.0000000000000000e+00,1.9325269174303068e-02,-9.6626345871515348e-02,0.0000000000000000e+00]

        tolerance = 1e-7

        for i in range(rhs.Size()):
                self.assertAlmostEqual(lhs[0, i], expected_LHS[i], delta=tolerance)

        # Confronto RHS
        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

        
if __name__ == '__main__':
    KratosUnittest.main()