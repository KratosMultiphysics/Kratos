import KratosMultiphysics as Kratos
import spatial_statistics_test_case

class SpatialStatisticsProcessTest(spatial_statistics_test_case.SpatialStatisticsTestCase):
    @classmethod
    def setUpClass(cls):
        super(SpatialStatisticsProcessTest, cls).setUpModelParts("spatial_statistics_process/spatial_statistics_process")

if __name__ == '__main__':
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()