from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import blade_beam_mapping_test

class BeamBladeMapping(blade_beam_mapping_test.BladeMappingTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "beam_mapper",
            "echo_level" : 0
        }""")
        super(BeamBladeMapping, cls).setUpMapper(mapper_params)
        cls.print_output = False

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()