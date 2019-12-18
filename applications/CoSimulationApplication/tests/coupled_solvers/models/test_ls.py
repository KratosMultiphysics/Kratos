import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import numpy as np


class TestModelLS(KratosUnittest.TestCase):
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

    def test_model_ls(self):
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

        parameter_file_name = "coupled_solvers/models/test_ls.json"
        with open(parameter_file_name, 'r') as parameter_file:
            settings = cs_data_structure.Parameters(parameter_file.read())


        # With reuse

        min_significant = settings["setting1"]["settings"]["min_significant"].GetDouble()
        q = settings["setting1"]["settings"]["q"].GetDouble()

        ls = cs_tools.CreateInstance(settings["setting1"])
        ls.size = m
        ls.Initialize()
        ls.InitializeSolutionStep()

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

        is_ready = ls.IsReady()
        self.assertFalse(is_ready)

        r.SetNumpyArray(r1)
        xt.SetNumpyArray(xt1)
        ls.Add(r, xt)

        self.assertTrue(ls.added)
        self.assertArrayEqual(r1, ls.rref)
        self.assertArrayEqual(xt1, ls.xtref)
        is_ready = ls.IsReady()
        self.assertFalse(is_ready)

        r.SetNumpyArray(r2)
        xt.SetNumpyArray(xt2)
        ls.Add(r, xt)

        self.assertTrue(ls.added)
        self.assertArrayEqual(r2, ls.rref)
        self.assertArrayEqual(xt2, ls.xtref)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 1))
        self.assertArrayEqual(v[:, 0], r2 - r1)
        self.assertEqual(w.shape, (m, 1))
        self.assertArrayEqual(w[:, 0], xt2 - xt1)
        is_ready = ls.IsReady()
        self.assertTrue(is_ready)
        r.SetNumpyArray(r3)
        dxt = ls.Predict(r)
        dxt_sol = [4.94444444444445, 0, -6.18055555555556, -3.70833333333333, -4.944444444444452]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)

        r.SetNumpyArray(r4)
        xt.SetNumpyArray(xt4)
        ls.Add(r, xt)

        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 2))
        self.assertArrayEqual(v[:, 0], r4 - r2)
        self.assertEqual(w.shape, (m, 2))
        self.assertArrayEqual(w[:, 0], xt4 - xt2)
        r.SetNumpyArray(r3)
        dxt = ls.Predict(r)
        dxt_sol = [3.55677154582763, -1.20861833105335, -7.26607387140903, -6.29343365253078, -6.37688098495212]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)

        r.SetNumpyArray(r5)
        xt.SetNumpyArray(xt5)
        ls.Add(r, xt)

        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 2))
        self.assertArrayEqual(v[:, 0], r5 - r4)
        self.assertEqual(w.shape, (m, 2))
        self.assertArrayEqual(w[:, 0], xt5 - xt4)
        r.SetNumpyArray(r3)
        dxt = ls.Predict(r)
        dxt_sol = [2.17629179331307, 7.43617021276596, 5.80395136778116, 5.96504559270517, -6.89209726443769]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 2))
        self.assertArrayEqual(v.T.flatten(), np.hstack((r5 - r4, r4 - r2)))
        self.assertEqual(w.shape, (m, 2))
        self.assertArrayEqual(w.T.flatten(), np.hstack((xt5 - xt4, xt4 - xt2)))

        r.SetNumpyArray(r6)
        xt.SetNumpyArray(xt6)
        ls.Add(r, xt)
        r.SetNumpyArray(r7)
        xt.SetNumpyArray(xt7)
        ls.Add(r, xt)
        r.SetNumpyArray(r8)
        xt.SetNumpyArray(xt8)
        ls.Add(r, xt)

        r.SetNumpyArray(r3)
        dxt = ls.Predict(r)
        dxt_sol = [7.67847593582868, 21.1203208556145, 21.6136363636358, 23.3649732620315, -1.94652406417126]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 5))
        self.assertArrayEqual(v.T.flatten(), np.hstack((r8 - r7, r7 - r6, r6 - r5, r5 - r4, r4 - r2)))
        self.assertEqual(w.shape, (m, 5))
        self.assertArrayEqual(w.T.flatten(), np.hstack((xt8 - xt7, xt7 - xt6, xt6 - xt5, xt5 - xt4, xt4 - xt2)))

        r.SetNumpyArray(r9)
        xt.SetNumpyArray(xt9)
        ls.Add(r, xt)
        r.SetNumpyArray(r10)
        xt.SetNumpyArray(xt10)
        ls.Add(r, xt)

        r.SetNumpyArray(r3)
        dxt = ls.Predict(r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 5))
        self.assertArrayEqual(v.T.flatten(), np.hstack((r10 - r9, r9 - r8, r8 - r7, r7 - r6, r6 - r5)))
        self.assertEqual(w.shape, (m, 5))
        self.assertArrayEqual(w.T.flatten(), np.hstack((xt10 - xt9, xt9 - xt8, xt8 - xt7, xt7 - xt6, xt6 - xt5)))

        r.SetNumpyArray(r13)
        xt.SetNumpyArray(xt13)
        ls.Add(r, xt)
        r.SetNumpyArray(r3)

        dxt = ls.Predict(r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 5))
        self.assertArrayEqual(v.T.flatten(), np.hstack((r10 - r9, r9 - r8, r8 - r7, r7 - r6, r6 - r5)))
        self.assertEqual(w.shape, (m, 5))
        self.assertArrayEqual(w.T.flatten(), np.hstack((xt10 - xt9, xt9 - xt8, xt8 - xt7, xt7 - xt6, xt6 - xt5)))

        v1 = ls.vcurr
        w1 = ls.wcurr

        # New solution step
        ls.FinalizeSolutionStep()
        self.assertArrayEqual(ls.vprev[0].flatten(), v1.flatten())
        self.assertArrayEqual(ls.wprev[0].flatten(), w1.flatten())
        ls.InitializeSolutionStep()
        self.assertIsNone(ls.rref)
        self.assertFalse(ls.added)
        self.assertEqual(ls.vcurr.shape, (m, 0))
        self.assertEqual(ls.wcurr.shape, (m, 0))
        is_ready = ls.IsReady()
        self.assertTrue(is_ready)

        r.SetNumpyArray(r11)
        xt.SetNumpyArray(xt11)
        ls.Add(r, xt)

        r.SetNumpyArray(r3)
        dxt = ls.Predict(r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        self.assertTrue(ls.added)
        self.assertEqual(ls.vcurr.shape, (m, 0))
        self.assertEqual(ls.wcurr.shape, (m, 0))

        r.SetNumpyArray(r12)
        xt.SetNumpyArray(xt12)
        ls.Add(r, xt)

        r.SetNumpyArray(r3)
        dxt = ls.Predict(r)
        dxt_sol = [-8.52029914529913, -3.96866096866096, -0.505163817663819, -0.426103988603991, -14.1239316239316]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 5))
        self.assertArrayEqual(v.T.flatten(), np.hstack((r12 - r11, r10 - r9, r9 - r8, r7 - r6, r6 - r5)))
        self.assertEqual(w.shape, (m, 5))
        self.assertArrayEqual(w.T.flatten(), np.hstack((xt12 - xt11, xt10 - xt9, xt9 - xt8, xt7 - xt6, xt6 - xt5)))

        v2 = ls.vcurr
        w2 = ls.wcurr

        # New solution step
        ls.FinalizeSolutionStep()
        self.assertArrayEqual(np.hstack(ls.vprev).flatten(), np.hstack([v2, v1[:, :2], v1[:, 3:]]).flatten())
        self.assertArrayEqual(np.hstack(ls.wprev).flatten(), np.hstack([w2, w1[:, :2], w1[:, 3:]]).flatten())
        ls.InitializeSolutionStep()

        # New solution step
        ls.FinalizeSolutionStep()
        self.assertArrayEqual(np.hstack(ls.vprev).flatten(), np.hstack([np.empty((m, 0)), v2]).flatten())
        self.assertArrayEqual(np.hstack(ls.wprev).flatten(), np.hstack([np.empty((m, 0)), w2]).flatten())
        self.assertEqual(len(ls.vprev), q)
        self.assertEqual(len(ls.wprev), q)
        ls.InitializeSolutionStep()

        # Without reuse

        min_significant = settings["setting2"]["settings"]["min_significant"].GetDouble()
        q = settings["setting2"]["settings"]["q"].GetDouble()

        ls = cs_tools.CreateInstance(settings["setting2"])
        ls.size = m
        ls.Initialize()
        ls.InitializeSolutionStep()

        r.SetNumpyArray(r1)
        xt.SetNumpyArray(xt1)
        ls.Add(r, xt)
        r.SetNumpyArray(r2)
        xt.SetNumpyArray(xt2)
        ls.Add(r, xt)
        r.SetNumpyArray(r4)
        xt.SetNumpyArray(xt4)
        ls.Add(r, xt)
        r.SetNumpyArray(r5)
        xt.SetNumpyArray(xt5)
        ls.Add(r, xt)
        r.SetNumpyArray(r6)
        xt.SetNumpyArray(xt6)
        ls.Add(r, xt)
        r.SetNumpyArray(r7)
        xt.SetNumpyArray(xt7)
        ls.Add(r, xt)
        r.SetNumpyArray(r8)
        xt.SetNumpyArray(xt8)
        ls.Add(r, xt)
        r.SetNumpyArray(r9)
        xt.SetNumpyArray(xt9)
        ls.Add(r, xt)
        r.SetNumpyArray(r10)
        xt.SetNumpyArray(xt10)
        ls.Add(r, xt)

        r.SetNumpyArray(r3)
        dxt = ls.Predict(r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        self.assertArrayAlmostEqual(dxt.GetNumpyArray(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 5))
        self.assertArrayEqual(v.T.flatten(), np.hstack((r10 - r9, r9 - r8, r8 - r7, r7 - r6, r6 - r5)))
        self.assertEqual(w.shape, (m, 5))
        self.assertArrayEqual(w.T.flatten(), np.hstack((xt10 - xt9, xt9 - xt8, xt8 - xt7, xt7 - xt6, xt6 - xt5)))

        # New solution step
        ls.FinalizeSolutionStep()
        self.assertArrayEqual(np.hstack(ls.vprev).flatten(), np.hstack((np.empty((m, 0)))).flatten())
        self.assertArrayEqual(np.hstack(ls.wprev).flatten(), np.hstack((np.empty((m, 0)))).flatten())
        ls.InitializeSolutionStep()
        self.assertIsNone(ls.rref)
        self.assertFalse(ls.added)
        self.assertEqual(ls.vcurr.shape, (m, 0))
        self.assertEqual(ls.wcurr.shape, (m, 0))
        is_ready = ls.IsReady()
        self.assertFalse(is_ready)

        r.SetNumpyArray(r11)
        xt.SetNumpyArray(xt11)
        ls.Add(r, xt)

        is_ready = ls.IsReady()
        self.assertFalse(is_ready)

        r.SetNumpyArray(r12)
        xt.SetNumpyArray(xt12)
        ls.Add(r, xt)

        is_ready = ls.IsReady()
        self.assertTrue(is_ready)


if __name__ == '__main__':
    KratosUnittest.main()
