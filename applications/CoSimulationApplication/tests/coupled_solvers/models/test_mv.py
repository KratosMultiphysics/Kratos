import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import numpy as np


class TestModelMV(KratosUnittest.TestCase):
    def assertArrayAlmostEqual(self, a1, a2):
        ls1 = list(a1)
        ls2 = list(a2)
        try:
            self.assertEqual(ls1, ls2)
        except AssertionError:
            for i in range(len(ls1)):
                self.assertAlmostEqual(ls1[i], ls2[i])

    def assertArrayEqual(self, a1, a2):
        self.assertEqual(list(a1), list(a2))

    def test_model_mv(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        m = 5
        dz = 2.0
        x = 10.0

        interface_settings = cs_data_structure.Parameters('{"wall": "AREA"}')

        # Create interface
        variable = vars(KM)["AREA"]
        model = cs_data_structure.Model()
        model_part = model.CreateModelPart("wall")
        model_part.AddNodalSolutionStepVariable(variable)
        for i in range(m):
            model_part.CreateNewNode(i, 0.0, 0.0, i * dz)
        step = 0
        for node in model_part.Nodes:
            node.SetSolutionStepValue(variable, step, x)
        interface = CoSimulationInterface(model, interface_settings)

        parameter_file_name = "coupled_solvers/models/test_mv.json"
        with open(parameter_file_name, 'r') as parameter_file:
            settings = cs_data_structure.Parameters(parameter_file.read())

        min_significant = settings["settings"]["min_significant"].GetDouble()

        mv = cs_tools.CreateInstance(settings)
        mv.size = m
        mv.Initialize()
        mv.InitializeSolutionStep()

        r = interface.deepcopy()
        xt = interface.deepcopy()
        r1 = np.array([1, 2, 3, 4, 5])
        xt1 = np.array([5, 4, 3, 2, 1])
        r2 = np.array([8, 5, 5, 5, 8])
        xt2 = np.array([1, 4, 8, 5, 5])
        r3 = np.array([7, 5, 6, 4, 3])
        r4 = np.array([1, 1, 7, 4, 0])
        xt4 = np.array([9, 7, 5, 8, 4])
        r5 = np.array([9, 5, 10, 6, 4])
        xt5 = np.array([5, 1, 2, 3, 9])
        r6 = np.array([7, 8, 1, 2, 3])
        xt6 = np.array([7, 5, 5, 1, 2])
        r7 = np.array([1, 2, 5, 1, 2])
        xt7 = np.array([4, 2, 1, 1, 2])
        r8 = np.array([6, 3, 9, 0, 3])
        xt8 = np.array([3, 1, 2, 3, 9])
        r9 = np.array([1, 3, 5, 0, 8])
        xt9 = np.array([8, 1, 5, 3, 9])
        r10 = np.array([1, 3, -5, 8, 8])
        xt10 = np.array([8, -9, 5, 3, -9])
        r11 = r1
        xt11 = xt1
        r12 = r2
        xt12 = xt2
        r13 = r10 * 0.95
        xt13 = xt10

        is_ready = mv.IsReady()
        self.assertFalse(is_ready)

        r.SetNumpyArray(r1)
        xt.SetNumpyArray(xt1)
        mv.Add(r, xt)

        self.assertTrue(mv.added)
        self.assertArrayEqual(r1, mv.rref)
        self.assertArrayEqual(xt1, mv.xtref)
        self.assertIsNone(mv.ncurr)
        is_ready = mv.IsReady()
        self.assertFalse(is_ready)

        r.SetNumpyArray(r2)
        xt.SetNumpyArray(xt2)
        mv.Add(r, xt)

        self.assertTrue(mv.added)
        self.assertArrayEqual(r2, mv.rref)
        self.assertArrayEqual(xt2, mv.xtref)
        self.assertEqual(mv.v.shape, (m, 1))
        self.assertArrayEqual(mv.v[:, 0], r2 - r1)
        self.assertEqual(mv.w.shape, (m, 1))
        self.assertArrayEqual(mv.w[:, 0], xt2 - xt1)
        n_sol = [-0.388888888888889, -0.166666666666667, -0.111111111111111, -0.0555555555555556, -0.166666666666667,
                0, 0, 0, 0, 0,
                0.486111111111111, 0.208333333333333, 0.138888888888889, 0.0694444444444445, 0.208333333333333,
                0.291666666666667, 0.125000000000000, 0.0833333333333333, 0.0416666666666667, 0.125000000000000,
                0.388888888888889, 0.166666666666667, 0.111111111111111, 0.0555555555555556, 0.166666666666667]
        self.assertArrayAlmostEqual(mv.ncurr.flatten(), n_sol)
        self.assertArrayEqual(mv.nprev.flatten(), np.zeros((m, m)).flatten())
        is_ready = mv.IsReady()
        self.assertTrue(is_ready)
        r.SetNumpyArray(r3)
        dxt = mv.Predict(r)
        dxt_sol = [4.94444444444445, 0, -6.18055555555556, -3.70833333333333, -4.944444444444452]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)

        r.SetNumpyArray(r4)
        xt.SetNumpyArray(xt4)
        mv.Add(r, xt)

        self.assertEqual(mv.v.shape, (m, 2))
        self.assertArrayEqual(mv.v[:, 0], r4 - r2)
        self.assertEqual(mv.w.shape, (m, 2))
        self.assertArrayEqual(mv.w[:, 0], xt4 - xt2)
        n_sol = [-0.306429548563612, -0.216142270861833, 0.251709986320109, -0.0437756497948016, -0.555403556771546,
                0.0718194254445965, -0.0430916552667578, 0.316005471956224, 0.0102599179206566, -0.338577291381669,
                0.550615595075240, 0.169630642954856, 0.422708618331053, 0.0786593707250342, -0.0957592339261285,
                0.445280437756498, 0.0328317373461012, 0.759233926128591, 0.0636114911080711, -0.599179206566348,
                0.474008207934337, 0.115595075239398, 0.485636114911081, 0.0677154582763338, -0.234610123119015]
        self.assertArrayAlmostEqual(mv.ncurr.flatten(), n_sol)
        r.SetNumpyArray(r3)
        dxt = mv.Predict(r)
        dxt_sol = [3.55677154582763, -1.20861833105335, -7.26607387140903, -6.29343365253078, -6.37688098495212]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)

        r.SetNumpyArray(r5)
        xt.SetNumpyArray(xt5)
        mv.Add(r, xt)

        self.assertEqual(mv.v.shape, (m, 2))
        self.assertArrayEqual(mv.v[:, 0], r5 - r4)
        self.assertEqual(mv.w.shape, (m, 2))
        self.assertArrayEqual(mv.w[:, 0], xt5 - xt4)
        n_sol = [-0.258792878853669, -0.180633955709944, 0.376899696048632, 0.0121580547112462, -0.590534085974816,
                -0.460486322188450, -0.200607902735562, -0.446808510638298, -0.159574468085106, 0.0364741641337384,
                -0.266391663048198, -0.0651324359531047, -0.729483282674772, -0.168693009118541, 0.479374728614850,
                -0.379722101606600, -0.171081198436822, -0.316109422492401, -0.123100303951368, -0.0208423795049935,
                0.395788102475033, 0.155449413808076, 0.541033434650456, 0.162613981762918, -0.184107685627443]
        self.assertArrayAlmostEqual(mv.ncurr.flatten(), n_sol)
        r.SetNumpyArray(r3)
        dxt = mv.Predict(r)
        dxt_sol = [2.17629179331307, 7.43617021276596, 5.80395136778116, 5.96504559270517, -6.89209726443769]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        self.assertEqual(mv.v.shape, (m, 2))
        self.assertArrayEqual(mv.v.T.flatten(), np.hstack((r5 - r4, r4 - r2)))
        self.assertEqual(mv.w.shape, (m, 2))
        self.assertArrayEqual(mv.w.T.flatten(), np.hstack((xt5 - xt4, xt4 - xt2)))

        r.SetNumpyArray(r6)
        xt.SetNumpyArray(xt6)
        mv.Add(r, xt)
        r.SetNumpyArray(r7)
        xt.SetNumpyArray(xt7)
        mv.Add(r, xt)
        r.SetNumpyArray(r8)
        xt.SetNumpyArray(xt8)
        mv.Add(r, xt)

        n_sol = [1.59692513368984, -1.67045454545455, -1.19117647058823, 0.612967914438510, -1.93649732620321,
                3.87433155080215, -5.22727272727274, -3.11764705882354, 0.660427807486646, -2.01336898395721,
                4.65909090909091, -6.15909090909093, -3.50000000000001, 0.568181818181834, -1.56818181818181,
                5.05213903743316, -7.02272727272729, -3.32352941176471, 0.736631016042804, -2.20721925133689,
                1.72192513368984, -1.54545454545455, 0.0588235294117627, -0.262032085561494, -0.561497326203206]
        self.assertArrayAlmostEqual(mv.ncurr.flatten(), n_sol)
        r.SetNumpyArray(r3)
        dxt = mv.Predict(r)
        dxt_sol = [7.67847593582868, 21.1203208556145, 21.6136363636358, 23.3649732620315, -1.94652406417126]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        self.assertEqual(mv.v.shape, (m, 5))
        self.assertArrayEqual(mv.v.T.flatten(), np.hstack((r8 - r7, r7 - r6, r6 - r5, r5 - r4, r4 - r2)))
        self.assertEqual(mv.w.shape, (m, 5))
        self.assertArrayEqual(mv.w.T.flatten(), np.hstack((xt8 - xt7, xt7 - xt6, xt6 - xt5, xt5 - xt4, xt4 - xt2)))

        r.SetNumpyArray(r9)
        xt.SetNumpyArray(xt9)
        mv.Add(r, xt)
        r.SetNumpyArray(r10)
        xt.SetNumpyArray(xt10)
        mv.Add(r, xt)

        r.SetNumpyArray(r3)
        dxt = mv.Predict(r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        self.assertEqual(mv.v.shape, (m, 5))
        self.assertArrayEqual(mv.v.T.flatten(), np.hstack((r10 - r9, r9 - r8, r8 - r7, r7 - r6, r6 - r5)))
        self.assertEqual(mv.w.shape, (m, 5))
        self.assertArrayEqual(mv.w.T.flatten(), np.hstack((xt10 - xt9, xt9 - xt8, xt8 - xt7, xt7 - xt6, xt6 - xt5)))

        r.SetNumpyArray(r13)
        xt.SetNumpyArray(xt13)
        mv.Add(r, xt)

        r.SetNumpyArray(r3)
        dxt = mv.Predict(r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        self.assertEqual(mv.v.shape, (m, 5))
        self.assertArrayEqual(mv.v.T.flatten(), np.hstack((r10 - r9, r9 - r8, r8 - r7, r7 - r6, r6 - r5)))
        self.assertEqual(mv.w.shape, (m, 5))
        self.assertArrayEqual(mv.w.T.flatten(), np.hstack((xt10 - xt9, xt9 - xt8, xt8 - xt7, xt7 - xt6, xt6 - xt5)))

        v = mv.v
        w = mv.w
        
        nprev = w @ np.linalg.inv(v.T @ v) @ v.T

        # New solution step
        mv.FinalizeSolutionStep()
        self.assertArrayEqual(mv.nprev.flatten(), nprev.flatten())
        mv.InitializeSolutionStep()
        self.assertIsNone(mv.rref)
        self.assertFalse(mv.added)
        self.assertEqual(mv.v.shape, (m, 0))
        self.assertEqual(mv.w.shape, (m, 0))
        is_ready = mv.IsReady()
        self.assertTrue(is_ready)
        self.assertArrayAlmostEqual(mv.ncurr.flatten(), nprev.flatten())
        r.SetNumpyArray(r3)
        dxt = mv.Predict(r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)

        r.SetNumpyArray(r11)
        xt.SetNumpyArray(xt11)
        mv.Add(r, xt)

        self.assertTrue(mv.added)
        self.assertEqual(mv.v.shape, (m, 0))
        self.assertEqual(mv.w.shape, (m, 0))
        r.SetNumpyArray(r12)
        xt.SetNumpyArray(xt12)
        mv.Add(r, xt)
        n_sol = [-1.07953029089939, 0.852682145716574, -0.00813984520949948, 0.0950093408059783, 0.306645316253002,
                -1.08484565430122, 2.54250066720043, 0.878480562227564, -0.192264923049552, -0.532759540966108,
                0.315863802152833, 0.444989324793168, -0.139022328974290, -0.215427897873855, 0.649152655457700,
                0.519159772262255, -0.566019482252470, -0.0682990837114143, -0.0941975802864517, 0.431578596210302,
                -0.394248732319190, 0.855284227381908, 1.23271950894049, -0.654857219108619, 0.794435548438751]
        self.assertArrayAlmostEqual(mv.ncurr.flatten(), n_sol)

        r.SetNumpyArray(r3)
        dxt = mv.Predict(r)
        dxt_sol = [2.04216706698691, -8.02212881416246, -4.68760564006761, -1.31217195979005, -8.67687483319990]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        self.assertEqual(mv.v.shape, (m, 1))
        self.assertArrayEqual(mv.v, r12 - r11)
        self.assertEqual(mv.w.shape, (m, 1))
        self.assertArrayEqual(mv.w, xt12 - xt11)


if __name__ == '__main__':
    KratosUnittest.main()
