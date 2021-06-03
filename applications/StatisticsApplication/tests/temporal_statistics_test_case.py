from KratosMultiphysics.StatisticsApplication.test_utilities import HistoricalRetrievalMethod
from KratosMultiphysics.StatisticsApplication.test_utilities import NonHistoricalRetrievalMethod

import statistics_test_case


class TemporalStatisticsTestCase(statistics_test_case.StatisticsTestCase):
    @classmethod
    def setUpClass(cls):
        super(TemporalStatisticsTestCase, cls).setUpModelParts(
            "spatial_statistics_process/spatial_statistics_process")

    @staticmethod
    def GetOutputMethod(container_name):
        if (container_name.endswith("_non_historical")):
            return NonHistoricalRetrievalMethod
        elif (container_name.endswith("_historical")):
            return HistoricalRetrievalMethod
        else:
            raise Exception("Unknown container_name. [ container_name = \"" +
                            str(container_name) + "\" ].")

    @staticmethod
    def GetInputMethod(container_name):
        if (container_name.startswith("nodal_historical")):
            return HistoricalRetrievalMethod
        else:
            return NonHistoricalRetrievalMethod

    def GetContainer(self, container_name):
        if (container_name.startswith("nodal")):
            return self.GetModelPart().Nodes
        elif (container_name.startswith("condition")):
            return self.GetModelPart().Conditions
        elif (container_name.startswith("element")):
            return self.GetModelPart().Elements
        else:
            raise Exception("Unknown container_name. [ container_name = \"" +
                            str(container_name) + "\" ].")

    def RunTemporalStatisticsTest(self, norm_type, container_name):
        raise Exception(
            "Calling base class RunTemporalStatisticsTest method in " +
            str(self.__class__.__name__) +
            ". Please implement it in derrived class")


class TemporalStatisticsValueTestCases(TemporalStatisticsTestCase):
    def testHistoricalHistoricalValueMethod(self):
        self.RunTemporalStatisticsTest("none", "nodal_historical_historical")

    def testHistoricalNonHistoricalValueMethod(self):
        self.RunTemporalStatisticsTest("none",
                                       "nodal_historical_non_historical")

    def testNodalNonHistoricalValueMethod(self):
        self.RunTemporalStatisticsTest("none", "nodal_non_historical")

    def testConditionNonHistoricalValueMethod(self):
        self.RunTemporalStatisticsTest("none", "condition_non_historical")

    def testElementNonHistoricalValueMethod(self):
        self.RunTemporalStatisticsTest("none", "element_non_historical")


class TemporalStatisticsNormTestCases(TemporalStatisticsTestCase):
    def testHistoricalHistoricalNormMethod(self):
        self.RunTemporalStatisticsTest("magnitude",
                                       "nodal_historical_historical")

    def testHistoricalNonHistoricalNormMethod(self):
        self.RunTemporalStatisticsTest("magnitude",
                                       "nodal_historical_non_historical")

    def testNodalNonHistoricalNormMethod(self):
        self.RunTemporalStatisticsTest("magnitude", "nodal_non_historical")

    def testConditionNonHistoricalNormMethod(self):
        self.RunTemporalStatisticsTest("magnitude", "condition_non_historical")

    def testElementNonHistoricalNormMethod(self):
        self.RunTemporalStatisticsTest("magnitude", "element_non_historical")
