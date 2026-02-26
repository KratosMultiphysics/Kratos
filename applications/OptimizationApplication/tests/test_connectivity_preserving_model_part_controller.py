import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.model_part_controllers.connectivity_preserving_model_part_controller import ConnectivityPreservingModelPartController


class TestConnectivityPreservingModelPartController(kratos_unittest.TestCase):
    def test_Duplication1(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("source")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        model_part.CreateNewNode(1, 0, 0, 0)
        model_part.CreateNewNode(2, 1, 0, 0)
        model_part.CreateNewNode(3, 1, 1, 0)
        model_part.CreateNewNode(4, 0, 1, 0)

        properties = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("Element3D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewElement("Element3D3N", 2, [1, 3, 4], properties)

        parameters = Kratos.Parameters("""{
            "transformation_settings": [
                {
                    "source_model_part_name"     : "source",
                    "destination_model_part_name": "destination",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : ""
                }
            ]
        }""")
        ConnectivityPreservingModelPartController(model, parameters).ImportModelPart()
        self.__CheckModelParts(model_part, model["destination"])

    def test_Duplication2(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("source")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        model_part.CreateNewNode(1, 0, 0, 0)
        model_part.CreateNewNode(2, 1, 0, 0)
        model_part.CreateNewNode(3, 1, 1, 0)
        model_part.CreateNewNode(4, 0, 1, 0)

        properties = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("Element3D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewElement("Element3D3N", 2, [1, 3, 4], properties)

        sub_model_part = model_part.CreateSubModelPart("sub_source")
        sub_model_part.AddElements([2])
        sub_model_part.AddNodes([1, 3, 4])

        parameters = Kratos.Parameters("""{
            "transformation_settings": [
                {
                    "source_model_part_name"     : "source.sub_source",
                    "destination_model_part_name": "destination.sub_source",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : ""
                }
            ]
        }""")
        ConnectivityPreservingModelPartController(model, parameters).ImportModelPart()
        self.__CheckModelParts(sub_model_part, model["destination"])
        self.__CheckModelParts(sub_model_part, model["destination.sub_source"])

    def test_Duplication3(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("source")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        model_part.CreateNewNode(1, 0, 0, 0)
        model_part.CreateNewNode(2, 1, 0, 0)
        model_part.CreateNewNode(3, 1, 1, 0)
        model_part.CreateNewNode(4, 0, 1, 0)

        properties = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("Element3D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewElement("Element3D3N", 2, [1, 3, 4], properties)

        sub_model_part_1 = model_part.CreateSubModelPart("sub_source_1")
        sub_model_part_1.AddElements([1])
        sub_model_part_1.AddNodes([1, 2, 3])

        sub_model_part_2 = model_part.CreateSubModelPart("sub_source_2")
        sub_model_part_2.AddElements([2])
        sub_model_part_2.AddNodes([1, 3, 4])

        parameters = Kratos.Parameters("""{
            "transformation_settings": [
                {
                    "source_model_part_name"     : "source.sub_source_1",
                    "destination_model_part_name": "destination.sub_source_1",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : ""
                },
                {
                    "source_model_part_name"     : "source.sub_source_2",
                    "destination_model_part_name": "destination.sub_source_2",
                    "destination_element_name"   : "HelmholtzVectorSurfaceElement3D3N",
                    "destination_condition_name" : ""
                }
            ]
        }""")
        ConnectivityPreservingModelPartController(model, parameters).ImportModelPart()
        self.__CheckModelParts(model_part, model["destination"])
        self.__CheckModelParts(sub_model_part_1, model["destination.sub_source_1"])
        self.__CheckModelParts(sub_model_part_2, model["destination.sub_source_2"])

    def test_Duplication4(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("source")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        model_part.CreateNewNode(1, 0, 0, 0)
        model_part.CreateNewNode(2, 1, 0, 0)
        model_part.CreateNewNode(3, 1, 1, 0)
        model_part.CreateNewNode(4, 0, 1, 0)

        properties = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("Element3D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewElement("Element3D3N", 2, [1, 3, 4], properties)
        model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [1, 3, 4], properties)

        sub_model_part_1 = model_part.CreateSubModelPart("sub_source_1")
        sub_model_part_1.AddElements([1])
        sub_model_part_1.AddConditions([1])
        sub_model_part_1.AddNodes([1, 2, 3])

        sub_model_part_2 = model_part.CreateSubModelPart("sub_source_2")
        sub_model_part_2.AddElements([2])
        sub_model_part_1.AddConditions([2])
        sub_model_part_2.AddNodes([1, 3, 4])

        parameters = Kratos.Parameters("""{
            "transformation_settings": [
                {
                    "source_model_part_name"     : "source.sub_source_1",
                    "destination_model_part_name": "destination.sub_source_1",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : "HelmholtzSurfaceShapeCondition3D3N"
                },
                {
                    "source_model_part_name"     : "source.sub_source_2",
                    "destination_model_part_name": "destination.sub_source_2",
                    "destination_element_name"   : "HelmholtzVectorSurfaceElement3D3N",
                    "destination_condition_name" : "HelmholtzSurfaceShapeCondition3D3N"
                }
            ]
        }""")
        ConnectivityPreservingModelPartController(model, parameters).ImportModelPart()
        self.__CheckModelParts(model_part, model["destination"])
        self.__CheckModelParts(sub_model_part_1, model["destination.sub_source_1"])
        self.__CheckModelParts(sub_model_part_2, model["destination.sub_source_2"])

    def test_Duplication5(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("source")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        model_part.CreateNewNode(1, 0, 0, 0)
        model_part.CreateNewNode(2, 1, 0, 0)
        model_part.CreateNewNode(3, 1, 1, 0)
        model_part.CreateNewNode(4, 0, 1, 0)

        properties = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("Element3D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewElement("Element3D3N", 2, [1, 3, 4], properties)

        sub_model_part = model_part.CreateSubModelPart("sub_source")
        sub_model_part.AddElements([2])
        sub_model_part.AddNodes([1, 3, 4])

        dest_model_part = model.CreateModelPart("destination")
        dest_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)

        parameters = Kratos.Parameters("""{
            "transformation_settings": [
                {
                    "source_model_part_name"     : "source.sub_source",
                    "destination_model_part_name": "destination.sub_source",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : ""
                }
            ]
        }""")
        ConnectivityPreservingModelPartController(model, parameters).ImportModelPart()
        self.__CheckModelParts(sub_model_part, model["destination.sub_source"])

    def test_Duplication6(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("source")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        model_part.CreateNewNode(1, 0, 0, 0)
        model_part.CreateNewNode(2, 1, 0, 0)
        model_part.CreateNewNode(3, 1, 1, 0)
        model_part.CreateNewNode(4, 0, 1, 0)

        properties = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("Element3D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewElement("Element3D3N", 2, [1, 3, 4], properties)

        sub_model_part = model_part.CreateSubModelPart("sub_source")
        sub_model_part.AddElements([2])
        sub_model_part.AddNodes([1, 3, 4])

        dest_model_part = model.CreateModelPart("destination")
        dest_model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)

        parameters = Kratos.Parameters("""{
            "transformation_settings": [
                {
                    "source_model_part_name"     : "source.sub_source",
                    "destination_model_part_name": "destination.sub_source",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : ""
                }
            ]
        }""")

        with self.assertRaises(RuntimeError):
            ConnectivityPreservingModelPartController(model, parameters).ImportModelPart()

    def test_Duplication7(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("source")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        model_part.CreateNewNode(1, 0, 0, 0)
        model_part.CreateNewNode(2, 1, 0, 0)
        model_part.CreateNewNode(3, 1, 1, 0)
        model_part.CreateNewNode(4, 0, 1, 0)
        model_part.CreateNewNode(5, 2, 0, 0)
        model_part.CreateNewNode(6, 2, 1, 0)

        properties = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("Element3D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewElement("Element3D3N", 2, [1, 3, 4], properties)
        model_part.CreateNewElement("Element3D3N", 3, [2, 5, 6], properties)
        model_part.CreateNewElement("Element3D3N", 4, [2, 6, 3], properties)

        sub_model_part_1 = model_part.CreateSubModelPart("sub_source_1")
        sub_model_part_1.AddElements([2,3,4])
        sub_model_part_1.AddNodes([1, 2, 3, 4, 5, 6])

        sub_sub_model_part = sub_model_part_1.CreateSubModelPart("sub_sub_source")
        sub_sub_model_part.AddElements([3,4])
        sub_sub_model_part.AddNodes([2, 3, 5, 6])

        sub_model_part_2 = model_part.CreateSubModelPart("sub_source_2")
        sub_model_part_2.AddElements([1])
        sub_model_part_2.AddNodes([1, 2, 3])

        parameters = Kratos.Parameters("""{
            "transformation_settings": [
                {
                    "source_model_part_name"     : "source.sub_source_1",
                    "destination_model_part_name": "destination.sub_source_1",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : ""
                },
                {
                    "source_model_part_name"     : "source.sub_source_1.sub_sub_source",
                    "destination_model_part_name": "destination.sub_source_1.sub_sub_source",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : ""
                },
                {
                    "source_model_part_name"     : "source.sub_source_2",
                    "destination_model_part_name": "destination.sub_source_2",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : ""
                }
            ]
        }""")
        ConnectivityPreservingModelPartController(model, parameters).ImportModelPart()
        self.__CheckModelParts(model_part, model["destination"])
        self.__CheckModelParts(sub_model_part_1, model["destination.sub_source_1"])
        self.__CheckModelParts(sub_sub_model_part, model["destination.sub_source_1.sub_sub_source"])
        self.__CheckModelParts(sub_model_part_2, model["destination.sub_source_2"])

    def test_Duplication8(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("source")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        model_part.CreateNewNode(1, 0, 0, 0)
        model_part.CreateNewNode(2, 1, 0, 0)
        model_part.CreateNewNode(3, 1, 1, 0)
        model_part.CreateNewNode(4, 0, 1, 0)
        model_part.CreateNewNode(5, 2, 0, 0)
        model_part.CreateNewNode(6, 2, 1, 0)

        properties = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("Element3D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewElement("Element3D3N", 2, [1, 3, 4], properties)
        model_part.CreateNewElement("Element3D3N", 3, [2, 5, 6], properties)
        model_part.CreateNewElement("Element3D3N", 4, [2, 6, 3], properties)

        sub_model_part_1 = model_part.CreateSubModelPart("sub_source_1")
        sub_model_part_1.AddElements([2,3,4])
        sub_model_part_1.AddNodes([1, 2, 3, 4, 5, 6])

        sub_sub_model_part = sub_model_part_1.CreateSubModelPart("sub_sub_source")
        sub_sub_model_part.AddElements([3,4])
        sub_sub_model_part.AddNodes([2, 3, 5, 6])

        sub_model_part_2 = model_part.CreateSubModelPart("sub_source_2")
        sub_model_part_2.AddElements([1])
        sub_model_part_2.AddNodes([1, 2, 3])

        parameters = Kratos.Parameters("""{
            "transformation_settings": [
                {
                    "source_model_part_name"     : "source.sub_source_1.sub_sub_source",
                    "destination_model_part_name": "destination.sub_source_1.sub_sub_source",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : ""
                },
                {
                    "source_model_part_name"     : "source.sub_source_2",
                    "destination_model_part_name": "destination.sub_source_2",
                    "destination_element_name"   : "HelmholtzSurfaceElement3D3N",
                    "destination_condition_name" : ""
                }
            ]
        }""")
        ConnectivityPreservingModelPartController(model, parameters).ImportModelPart()
        self.__CheckModelParts(sub_sub_model_part, model["destination.sub_source_1"])
        self.__CheckModelParts(sub_sub_model_part, model["destination.sub_source_1.sub_sub_source"])
        self.__CheckModelParts(sub_model_part_2, model["destination.sub_source_2"])

    def __CheckModelParts(self, model_part_1: Kratos.ModelPart, model_part_2: Kratos.ModelPart):
        self.assertEqual(model_part_1.NumberOfElements(), model_part_2.NumberOfElements())
        self.assertEqual(model_part_1.NumberOfConditions(), model_part_2.NumberOfConditions())

        for element_1 in model_part_1.Elements:
            element_1: Kratos.Element
            self.assertTrue(model_part_2.HasElement(element_1.Id))
            element_2 = model_part_2.GetElement(element_1.Id)
            self.assertEqual(element_1.GetGeometry(), element_2.GetGeometry())

        for condition_1 in model_part_1.Conditions:
            self.assertTrue(model_part_2.HasCondition(condition_1.Id))
            condition_2 = model_part_2.GetCondition(condition_1.Id)
            self.assertEqual(condition_1.GetGeometry(), condition_2.GetGeometry())

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()

