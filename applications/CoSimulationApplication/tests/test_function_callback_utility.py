import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.function_callback_utility import GenericCallFunction

from math import pi

class TestGenericCallFunction(KratosUnittest.TestCase):

    def test_value(self):
        val = 125.964
        func_str = str(val)
        self.assertAlmostEqual(val, GenericCallFunction(func_str))

    def test_value_manipulation(self):
        val = 125.964
        val2 = 88.787
        func_str = "{}+{}".format(val, val2)
        self.assertAlmostEqual(val+val2, GenericCallFunction(func_str))

        func_str = "{}  +  {}".format(val, val2)
        self.assertAlmostEqual(val+val2, GenericCallFunction(func_str))

    def test_math_constants(self):
        val = 125.964
        func_str = "{}+pi".format(val)
        self.assertAlmostEqual(val+pi, GenericCallFunction(func_str))

    def test_scope(self):
        val = 125.964
        val_t = 1.556
        val_tr = -2.89
        scope = {"t" : val_t, "tr" : val_tr}
        func_str = "{}*t+tr".format(val)
        self.assertAlmostEqual(val*val_t+val_tr, GenericCallFunction(func_str, scope))

    def test_complex_eval(self):
        val_t = 0.6289
        scope = {"t" : val_t}
        func_str = "(sin((t/0.006289)/100*pi/2))**2"
        self.assertAlmostEqual(1.0, GenericCallFunction(func_str, scope))


if __name__ == '__main__':
    KratosUnittest.main()
