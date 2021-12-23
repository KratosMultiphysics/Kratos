from __future__ import print_function, absolute_import, division
import KratosMultiphysics


import KratosMultiphysics.CableNetApplication as CableNetApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CableNetApplication import edge_cable_element_process

class EdgeCableProcessTests(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def test_sorted_node_list(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("solid_part")
        self._create_nodes(model_part)

        project_settings = KratosMultiphysics.Parameters("""
        {
            "edge_sub_model_part_name"  : "solid_part",
            "element_type"              : "cable",
            "node_id_order"             : [],
            "element_id"                : 1,
            "property_id"               : 1
        }
        """)

        my_process = edge_cable_element_process.EdgeCableElementProcess(current_model,project_settings)
        my_list = my_process.CreateCorrectNodeOrder()

        original_order = [47, 72, 4, 14, 12, 17, 55, 61, 78, 19, 46, 63, 98, 33, 31, 25, 81, 28, 30, 24, 84, 95, 5,
         92, 70, 6, 15, 67, 52, 21, 2, 79, 76, 1, 32, 64, 91, 87, 48, 73, 45, 94, 43, 74, 62, 93, 20, 8, 56, 42, 29,
          49, 37, 13, 11, 99, 39, 26, 89, 53, 7, 96, 35, 66, 10, 77, 83, 3, 40, 90, 59, 51, 57, 80, 27, 68, 16, 34,
           60, 86, 44, 9, 65, 54, 69, 88, 36, 41, 75, 22, 18, 23, 38, 85, 58, 100, 50, 82, 97, 71]

        for node_id_calculated, node_id_result in zip(my_list,original_order):
            self.assertEqual(node_id_calculated, node_id_result)

    def _create_nodes(self,model_part):
        model_part.CreateNewNode(1,132.0,79.2,-343.2)
        model_part.CreateNewNode(2,138.0,82.8,-358.8)
        model_part.CreateNewNode(3,64.0,38.4,-166.4)
        model_part.CreateNewNode(4,194.0,116.39999999999999,-504.40000000000003)
        model_part.CreateNewNode(5,154.0,92.39999999999999,-400.40000000000003)
        model_part.CreateNewNode(6,148.0,88.8,-384.8)
        model_part.CreateNewNode(7,78.0,46.8,-202.8)
        model_part.CreateNewNode(8,104.0,62.4,-270.40000000000003)
        model_part.CreateNewNode(9,36.0,21.599999999999998,-93.60000000000001)
        model_part.CreateNewNode(10,70.0,42.0,-182.0)
        model_part.CreateNewNode(11,90.0,54.0,-234.0)
        model_part.CreateNewNode(12,190.0,114.0,-494.0)
        model_part.CreateNewNode(13,92.0,55.199999999999996,-239.20000000000002)
        model_part.CreateNewNode(14,192.0,115.19999999999999,-499.20000000000005)
        model_part.CreateNewNode(15,146.0,87.6,-379.6)
        model_part.CreateNewNode(16,46.0,27.599999999999998,-119.60000000000001)
        model_part.CreateNewNode(17,188.0,112.8,-488.8)
        model_part.CreateNewNode(18,18.0,10.799999999999999,-46.800000000000004)
        model_part.CreateNewNode(19,180.0,108.0,-468.0)
        model_part.CreateNewNode(20,106.0,63.599999999999994,-275.6)
        model_part.CreateNewNode(21,140.0,84.0,-364.0)
        model_part.CreateNewNode(22,20.0,12.0,-52.0)
        model_part.CreateNewNode(23,16.0,9.6,-41.6)
        model_part.CreateNewNode(24,160.0,96.0,-416.0)
        model_part.CreateNewNode(25,168.0,100.8,-436.8)
        model_part.CreateNewNode(26,84.0,50.4,-218.4)
        model_part.CreateNewNode(27,50.0,30.0,-130.0)
        model_part.CreateNewNode(28,164.0,98.39999999999999,-426.40000000000003)
        model_part.CreateNewNode(29,98.0,58.8,-254.8)
        model_part.CreateNewNode(30,162.0,97.2,-421.2)
        model_part.CreateNewNode(31,170.0,102.0,-442.0)
        model_part.CreateNewNode(32,130.0,78.0,-338.0)
        model_part.CreateNewNode(33,172.0,103.2,-447.2)
        model_part.CreateNewNode(34,44.0,26.4,-114.4)
        model_part.CreateNewNode(35,74.0,44.4,-192.4)
        model_part.CreateNewNode(36,26.0,15.6,-67.60000000000001)
        model_part.CreateNewNode(37,94.0,56.4,-244.4)
        model_part.CreateNewNode(38,14.0,8.4,-36.4)
        model_part.CreateNewNode(39,86.0,51.6,-223.6)
        model_part.CreateNewNode(40,62.0,37.199999999999996,-161.20000000000002)
        model_part.CreateNewNode(41,24.0,14.399999999999999,-62.400000000000006)
        model_part.CreateNewNode(42,100.0,60.0,-260.0)
        model_part.CreateNewNode(43,114.0,68.39999999999999,-296.40000000000003)
        model_part.CreateNewNode(44,38.0,22.8,-98.8)
        model_part.CreateNewNode(45,118.0,70.8,-306.8)
        model_part.CreateNewNode(46,178.0,106.8,-462.8)
        model_part.CreateNewNode(47,198.0,118.8,-514.8000000000001)
        model_part.CreateNewNode(48,122.0,73.2,-317.2)
        model_part.CreateNewNode(49,96.0,57.599999999999994,-249.60000000000002)
        model_part.CreateNewNode(50,6.0,3.5999999999999996,-15.600000000000001)
        model_part.CreateNewNode(51,56.0,33.6,-145.6)
        model_part.CreateNewNode(52,142.0,85.2,-369.2)
        model_part.CreateNewNode(53,80.0,48.0,-208.0)
        model_part.CreateNewNode(54,32.0,19.2,-83.2)
        model_part.CreateNewNode(55,186.0,111.6,-483.6)
        model_part.CreateNewNode(56,102.0,61.199999999999996,-265.2)
        model_part.CreateNewNode(57,54.0,32.4,-140.4)
        model_part.CreateNewNode(58,10.0,6.0,-26.0)
        model_part.CreateNewNode(59,58.0,34.8,-150.8)
        model_part.CreateNewNode(60,42.0,25.2,-109.2)
        model_part.CreateNewNode(61,184.0,110.39999999999999,-478.40000000000003)
        model_part.CreateNewNode(62,110.0,66.0,-286.0)
        model_part.CreateNewNode(63,176.0,105.6,-457.6)
        model_part.CreateNewNode(64,128.0,76.8,-332.8)
        model_part.CreateNewNode(65,34.0,20.4,-88.4)
        model_part.CreateNewNode(66,72.0,43.199999999999996,-187.20000000000002)
        model_part.CreateNewNode(67,144.0,86.39999999999999,-374.40000000000003)
        model_part.CreateNewNode(68,48.0,28.799999999999997,-124.80000000000001)
        model_part.CreateNewNode(69,30.0,18.0,-78.0)
        model_part.CreateNewNode(70,150.0,90.0,-390.0)
        model_part.CreateNewNode(71,0.0,0.0,-0.0)
        model_part.CreateNewNode(72,196.0,117.6,-509.6)
        model_part.CreateNewNode(73,120.0,72.0,-312.0)
        model_part.CreateNewNode(74,112.0,67.2,-291.2)
        model_part.CreateNewNode(75,22.0,13.2,-57.2)
        model_part.CreateNewNode(76,134.0,80.39999999999999,-348.40000000000003)
        model_part.CreateNewNode(77,68.0,40.8,-176.8)
        model_part.CreateNewNode(78,182.0,109.2,-473.2)
        model_part.CreateNewNode(79,136.0,81.6,-353.6)
        model_part.CreateNewNode(80,52.0,31.2,-135.20000000000002)
        model_part.CreateNewNode(81,166.0,99.6,-431.6)
        model_part.CreateNewNode(82,4.0,2.4,-10.4)
        model_part.CreateNewNode(83,66.0,39.6,-171.6)
        model_part.CreateNewNode(84,158.0,94.8,-410.8)
        model_part.CreateNewNode(85,12.0,7.199999999999999,-31.200000000000003)
        model_part.CreateNewNode(86,40.0,24.0,-104.0)
        model_part.CreateNewNode(87,124.0,74.39999999999999,-322.40000000000003)
        model_part.CreateNewNode(88,28.0,16.8,-72.8)
        model_part.CreateNewNode(89,82.0,49.199999999999996,-213.20000000000002)
        model_part.CreateNewNode(90,60.0,36.0,-156.0)
        model_part.CreateNewNode(91,126.0,75.6,-327.6)
        model_part.CreateNewNode(92,152.0,91.2,-395.2)
        model_part.CreateNewNode(93,108.0,64.8,-280.8)
        model_part.CreateNewNode(94,116.0,69.6,-301.6)
        model_part.CreateNewNode(95,156.0,93.6,-405.6)
        model_part.CreateNewNode(96,76.0,45.6,-197.6)
        model_part.CreateNewNode(97,2.0,1.2,-5.2)
        model_part.CreateNewNode(98,174.0,104.39999999999999,-452.40000000000003)
        model_part.CreateNewNode(99,88.0,52.8,-228.8)
        model_part.CreateNewNode(100,8.0,4.8,-20.8)


if __name__ == '__main__':
    KratosUnittest.main()
