from __future__ import print_function, absolute_import, division

import sys
import types
import KratosMultiphysics
import KratosMultiphysics.decorator_injector as DecroatorInjector
import KratosMultiphysics.KratosUnittest as KratosUnittest

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def dummy_decorator_repeat(func):
    def wrapper(*args, **kwargs):
        print("Wrapped Call")
        func(*args, **kwargs)
        func(*args, **kwargs)
    return wrapper

class TestNodalDivergence(KratosUnittest.TestCase):

    def setUp(self):
        def add_once(input):
            input[0] += 1

        @dummy_decorator_repeat
        def add_once_decorated(input):
            input[0] += 1

        virtual_module = types.ModuleType('VirtualModule', 'Module created to provide a context for tests')
        virtual_module.__dict__.update({
            "add_once":add_once,
            "add_once_decorated":add_once_decorated
        })

        sys.modules['VirtualModule'] = virtual_module
        self.virtual_module = virtual_module

    def test_virtual_module(self):
        a = [0]
        self.virtual_module.add_once(a)
        self.assertEqual(a[0], 1)

    def test_virtual_module_decorated(self):
        a = [0]
        self.virtual_module.add_once_decorated(a)
        self.assertEqual(a[0], 2)

    def test_virtual_module_injected(self):
        DecroatorInjector.InjectIntoModule(self.virtual_module, dummy_decorator_repeat, lambda *x: True)

        a = [0]
        self.virtual_module.add_once(a)
        self.assertEqual(a[0], 2)

    def test_virtual_module_decorated__injected(self):
        DecroatorInjector.InjectIntoModule(self.virtual_module, dummy_decorator_repeat, lambda *x: True)

        a = [0]
        self.virtual_module.add_once_decorated(a)
        self.assertEqual(a[0], 4)

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
