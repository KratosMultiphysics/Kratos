import KratosMultiphysics as KM
import test_blade_beam_mapping

class BeamBladeMapping(test_blade_beam_mapping.BladeMappingTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "beam_mapper",
            "search_iterations"        : 30,
            "use_corotation": true,
            "echo_level" : 0
        }""")
        super(BeamBladeMapping, cls).setUpMapper(mapper_params)
        cls.print_output = False

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()