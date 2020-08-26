from __future__ import print_function, absolute_import, division

import sys
import types

import KratosMultiphysics
import KratosMultiphysics.profiler as KratosProfiler
import KratosMultiphysics.KratosUnittest as KratosUnittest

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestProfiler(KratosUnittest.TestCase):

    def add_once(input):
            input[0] += 1

    def setUp(self):
        def add_once(input):
            input[0] += 1

        # Generate a virtual module. This should be equivalent to "import VirtualModule"
        virtual_module = types.ModuleType('VirtualModule', 'Module created to provide a context for the decorator injector')
        virtual_module.__dict__.update({
            "add_once":add_once,
        })

        sys.modules['VirtualModule'] = virtual_module
        self.virtual_module = virtual_module

    def test_virtual_module(self):
        p = KratosProfiler.Profiler()
        p.ProfileAllModule(self.virtual_module)

        a = [0]
        self.virtual_module.add_once(a)
        self.assertEqual(a[0], 1)

        p.PrintResults()

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
