# We import the libraries
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestFindConditionsNeighboursProcess(KratosUnittest.TestCase):
    def tearDown(self):
        pass

    def setUp(self):
        # We create the model part
        self.current_model = KratosMultiphysics.Model()
        self.model_part = self.current_model.CreateModelPart("model_part")
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)

        input_mdpa = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_skin")
        model_part_io = KratosMultiphysics.ModelPartIO(input_mdpa)
        model_part_io.ReadModelPart(self.model_part)

        parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"    : "please_provide_model_part_name",
            "average_conditions" : 10
        }
        """)
        parameters["model_part_name"].SetString(self.model_part.Name)
        self.process = KratosMultiphysics.FindConditionsNeighboursProcess(self.current_model, parameters)

    def test_sphere(self):
        self.process.Execute()
        ref_solution = {
            1: [116, 24, 9],
            2: [106, 17, 93],
            3: [4, 21, 6],
            4: [3, 5, 124],
            5: [43, 4, 6],
            6: [33, 5, 3],
            7: [112, 93, 18],
            8: [100, 12, 10],
            9: [1, 100, 10],
            10: [39, 9, 8],
            11: [74, 50, 70],
            12: [117, 14, 8],
            13: [86, 28, 102],
            14: [98, 26, 12],
            15: [122, 46, 43],
            16: [50, 23, 56],
            17: [92, 27, 2],
            18: [19, 65, 7],
            19: [18, 93, 115],
            20: [65, 52, 31],
            21: [3, 118, 82],
            22: [102, 48, 55],
            23: [16, 72, 105],
            24: [1, 81, 30],
            25: [128, 36, 45],
            26: [14, 48, 88],
            27: [17, 72, 95],
            28: [13, 52, 114],
            29: [30, 119, 80],
            30: [24, 29, 84],
            31: [20, 86, 96],
            32: [113, 116, 76],
            33: [34, 44, 6],
            34: [69, 42, 33],
            35: [95, 49, 115],
            36: [25, 118, 57],
            37: [121, 71, 111],
            38: [124, 127, 57],
            39: [40, 71, 10],
            40: [88, 59, 39],
            41: [86, 55, 66],
            42: [107, 34, 87],
            43: [119, 15, 5],
            44: [33, 73, 80],
            45: [127, 25, 57],
            46: [15, 68, 79],
            47: [112, 108, 67],
            48: [26, 22, 62],
            49: [74, 35, 54],
            50: [16, 75, 11],
            51: [115, 77, 60],
            52: [28, 20, 60],
            53: [123, 120, 78],
            54: [49, 95, 72],
            55: [22, 98, 41],
            56: [109, 113, 16],
            57: [38, 45, 36],
            58: [64, 61, 63],
            59: [111, 40, 62],
            60: [51, 114, 52],
            61: [110, 91, 58],
            62: [59, 88, 48],
            63: [109, 105, 58],
            64: [89, 110, 58],
            65: [18, 20, 128],
            66: [41, 78, 120],
            67: [47, 79, 68],
            68: [46, 106, 67],
            69: [34, 94, 87],
            70: [102, 11, 75],
            71: [39, 37, 76],
            72: [27, 23, 54],
            73: [84, 44, 107],
            74: [77, 49, 11],
            75: [50, 111, 70],
            76: [71, 121, 32],
            77: [114, 51, 74],
            78: [90, 53, 66],
            79: [108, 46, 67],
            80: [44, 84, 29],
            81: [24, 91, 85],
            82: [125, 94, 21],
            83: [94, 117, 87],
            84: [73, 30, 80],
            85: [119, 81, 126],
            86: [13, 41, 31],
            87: [42, 69, 83],
            88: [40, 26, 62],
            89: [103, 64, 92],
            90: [99, 78, 98],
            91: [97, 61, 81],
            92: [17, 89, 105],
            93: [2, 19, 7],
            94: [82, 83, 69],
            95: [35, 27, 54],
            96: [31, 120, 104],
            97: [116, 109, 91],
            98: [14, 90, 55],
            99: [117, 125, 90],
            100: [107, 8, 9],
            101: [118, 104, 123],
            102: [70, 22, 13],
            103: [106, 122, 89],
            104: [128, 96, 101],
            105: [23, 92, 63],
            106: [2, 68, 103],
            107: [73, 42, 100],
            108: [79, 47, 124],
            109: [56, 63, 97],
            110: [61, 64, 126],
            111: [37, 59, 75],
            112: [127, 47, 7],
            113: [56, 32, 121],
            114: [77, 28, 60],
            115: [35, 51, 19],
            116: [1, 32, 97],
            117: [12, 83, 99],
            118: [21, 36, 101],
            119: [43, 29, 85],
            120: [53, 96, 66],
            121: [37, 113, 76],
            122: [15, 126, 103],
            123: [125, 101, 53],
            124: [108, 38, 4],
            125: [82, 123, 99],
            126: [122, 85, 110],
            127: [45, 38, 112],
            128: [65, 104, 25]
        }
        solution = KratosMultiphysics.AuxiliarModelPartUtilities(self.model_part).RetrieveConditionsNeighbourConditionsIds()
        for cond in self.model_part.Conditions:
            id = cond.Id
            partial_solution = solution[id]
            partial_ref_solution = ref_solution[id]
            self.assertEqual(partial_solution, partial_ref_solution)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
