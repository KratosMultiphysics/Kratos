from __future__ import print_function, absolute_import, division

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils

dependencies_are_available = KratosUtils.CheckIfApplicationsAvailable("StructuralMechanicsApplication")
if dependencies_are_available:
    import KratosMultiphysics.StructuralMechanicsApplication

class TestSpecificationsUtilities(KratosUnittest.TestCase):
    def test_specifications_utilities_elements(self):
        current_model = KratosMultiphysics.Model()

        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)
        KratosMultiphysics.SpecificationsUtilities.AddMissingVariables(model_part)
        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)
        
    def test_specifications_utilities_conditions(self):
        current_model = KratosMultiphysics.Model()

        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.CreateNewCondition("Condition3D", 1, [1,2,3], model_part.GetProperties()[1])
        
        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)
        KratosMultiphysics.SpecificationsUtilities.AddMissingVariables(model_part)
        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)

    @KratosUnittest.skipUnless(dependencies_are_available,"StructuralMechanicsApplication is not available")
    def test_specifications_utilities_elements_dependencies(self):
        current_model = KratosMultiphysics.Model()

        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.CreateNewElement("SmallDisplacementElement2D3N", 1, [1,2,3], model_part.GetProperties()[1])
        
        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)
        KratosMultiphysics.SpecificationsUtilities.AddMissingVariables(model_part)
        #self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 3)

    @KratosUnittest.skipUnless(dependencies_are_available,"StructuralMechanicsApplication is not available")
    def test_specifications_utilities_conditions_dependencies(self):
        current_model = KratosMultiphysics.Model()

        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.CreateNewCondition("SurfaceLoadCondition3D3N", 1, [1,2,3], model_part.GetProperties()[1])
        
        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)
        KratosMultiphysics.SpecificationsUtilities.AddMissingVariables(model_part)
        #self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 3)
        
if __name__ == '__main__':
    KratosUnittest.main()
