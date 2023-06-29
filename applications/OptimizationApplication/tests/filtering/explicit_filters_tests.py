import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestExplicitVertexMorphingFilter(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        with kratos_unittest.WorkFolderScope(".", __file__):
            ReadModelPart("solid", cls.model_part)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id + 1, node.Id + 2]))

        for condition in cls.model_part.Conditions:
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([condition.Id, condition.Id + 1, condition.Id + 2]))

        for element in cls.model_part.Elements:
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([element.Id, element.Id + 1, element.Id + 2]))

    def test_NodalExplicitVertexMorphingFilterCondition(self):
        radius = 2.0
        model_part = self.model_part.GetSubModelPart("design")
        vm_filter = KratosOA.NodalExplicitVertexMorphingFilter(model_part, "linear", 1000)

        filter_radius = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, radius)
        vm_filter.SetFilterRadius(filter_radius)

        unfiltered_field = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.VariableExpressionIO.Read(unfiltered_field, Kratos.VELOCITY, True)
        filtered_data = vm_filter.FilterField(unfiltered_field)

        ref_values = numpy.array([
                [ 7.020893792141365,   8.020893792141365,   9.020893792141365 ],
                [ 7.151659899691424,   8.151659899691424,   9.151659899691424 ],
                [ 7.40711302540754,    8.40711302540754,    9.407113025407538 ],
                [ 7.537879132957599,   8.537879132957599,   9.537879132957599 ],
                [ 7.6542761484549775,  8.654276148454978,   9.654276148454978 ],
                [ 7.790712885995792,   8.790712885995795,   9.790712885995791 ],
                [ 8.057243653947024,   9.057243653947022,  10.057243653947024 ],
                [ 8.193680391487838,   9.193680391487836,  10.193680391487836 ],
                [ 7.493458985567711,   8.493458985567711,   9.49345898556771  ],
                [ 7.636340089562771,   8.636340089562772,   9.636340089562772 ],
                [ 7.975658423143609,   8.975658423143608,   9.97565842314361  ],
                [ 7.896005112693083,   8.896005112693082,   9.896005112693082 ],
                [ 7.825347527844546,   8.825347527844546,   9.825347527844546 ]
            ])

        self.assertAlmostEqual(numpy.linalg.norm(filtered_data.Evaluate() - ref_values), 0.0, 12)

        filtered_data = vm_filter.FilterIntegratedField(unfiltered_field)
        ref_values = numpy.array([
                [14.041787584282735, 16.041787584282737, 18.041787584282737],
                [14.303319799382853, 16.30331979938285,  18.303319799382855],
                [14.814226050815082, 16.814226050815087, 18.814226050815083],
                [15.075758265915203, 17.0757582659152,   19.075758265915205],
                [22.96282844536494,  25.96282844536494,  28.962828445364945],
                [23.372138657987385, 26.372138657987385, 29.372138657987385],
                [24.171730961841078, 27.171730961841085, 30.17173096184108 ],
                [24.581041174463518, 27.581041174463525, 30.581041174463525],
                [22.480376956703143, 25.480376956703143, 28.480376956703147],
                [22.909020268688323, 25.90902026868833,  28.909020268688323],
                [23.92697526943084,  26.92697526943084,  29.92697526943084 ],
                [23.688015338079257, 26.68801533807926,  29.688015338079257],
                [23.476042583533644, 26.47604258353364,  29.476042583533648]
            ])
        self.assertAlmostEqual(numpy.linalg.norm(filtered_data.Evaluate() - ref_values), 0.0, 12)

    def test_NodalExplicitVertexMorphingFilterElement(self):
        radius = 2.0
        model_part = self.model_part.GetSubModelPart("structure")
        vm_filter = KratosOA.NodalExplicitVertexMorphingFilter(model_part, "linear", 1000)

        filter_radius = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, radius)
        vm_filter.SetFilterRadius(filter_radius)

        unfiltered_field = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.VariableExpressionIO.Read(unfiltered_field, Kratos.VELOCITY, True)
        filtered_data = vm_filter.FilterField(unfiltered_field)

        ref_values = numpy.array([
                [ 7.195467517797387,  8.195467517797386,  9.195467517797386],
                [ 7.318570814762067,  8.318570814762067,  9.318570814762067],
                [ 7.559054550072059,  8.559054550072057,  9.559054550072059],
                [ 7.682157847036736,  8.682157847036736,  9.682157847036738],
                [ 7.883514769564958,  8.883514769564957,  9.883514769564957],
                [ 8.006618066529636,  9.006618066529636, 10.006618066529636],
                [ 8.24710180183963,   9.24710180183963,  10.24710180183963 ],
                [ 8.370205098804309,  9.370205098804307, 10.370205098804309],
                [ 7.645886355496233,  8.645886355496232,  9.645886355496232],
                [ 8.210354586583765,  9.210354586583765, 10.210354586583765],
                [ 7.822178630754171,  8.822178630754172,  9.822178630754172],
                [ 8.134818666059648,  9.134818666059648, 10.134818666059646],
                [ 8.061427954957626,  9.061427954957626, 10.061427954957626],
                [ 7.996325696590015,  8.996325696590015,  9.996325696590015]
            ])

        self.assertAlmostEqual(numpy.linalg.norm(filtered_data.Evaluate() - ref_values), 0.0, 12)

        filtered_data = vm_filter.FilterIntegratedField(unfiltered_field)
        ref_values = numpy.array([
                [172.69122042713727,  196.69122042713727,  220.69122042713727 ],
                [175.6456995542896,   199.6456995542896,   223.64569955428962 ],
                [181.4173092017294,   205.4173092017294,   229.4173092017294  ],
                [184.37178832888168,  208.37178832888168,  232.37178832888168 ],
                [189.20435446955898,  213.20435446955898,  237.20435446955904 ],
                [192.15883359671125,  216.15883359671128,  240.15883359671125 ],
                [197.93044324415112,  221.93044324415112,  245.93044324415115 ],
                [200.8849223713034,   224.8849223713034,   248.88492237130342 ],
                [ 73.40050901276382,   83.00050901276381,   92.60050901276382 ],
                [ 78.81940403120414,   88.41940403120415,   98.01940403120413 ],
                [ 62.577429046033366,  70.57742904603337,   78.57742904603337 ],
                [ 65.07854932847718,   73.07854932847718,   81.07854932847717 ],
                [ 77.38970836759322,   86.98970836759321,   96.5897083675932  ],
                [ 76.76472668726413,   86.36472668726414,   95.96472668726415 ]
            ])
        self.assertAlmostEqual(numpy.linalg.norm(filtered_data.Evaluate() - ref_values), 0.0, 12)

    def test_ConditionExplicitVertexMorphingFilter(self):
        radius = 1.0
        vm_filter = KratosOA.ConditionExplicitVertexMorphingFilter(self.model_part, "linear", 1000)

        filter_radius = Kratos.Expression.ConditionExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, radius)
        vm_filter.SetFilterRadius(filter_radius)

        unfiltered_field = Kratos.Expression.ConditionExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(unfiltered_field, Kratos.VELOCITY)
        filtered_data = vm_filter.FilterField(unfiltered_field)

        ref_values = numpy.array([
            [ 6.758330624355697,  7.758330624355698,  8.758330624355697],
            [ 8.563643579842978,  9.563643579842978, 10.563643579842978],
            [ 7.694047957950673,  8.694047957950673,  9.694047957950675],
            [ 8.292863410975192,  9.292863410975192, 10.292863410975192],
            [ 7.268144993577022,  8.268144993577023,  9.268144993577023],
            [ 9.966790169575301, 10.9667901695753,   11.966790169575303],
            [ 8.919588772181047,  9.919588772181049, 10.919588772181045],
            [10.194739547987977, 11.194739547987977, 12.194739547987975],
            [ 9.497741144178006, 10.497741144178008, 11.497741144178008],
            [12.302092912247636, 13.302092912247637, 14.302092912247636],
            [11.254891514853382, 12.254891514853382, 13.254891514853384],
            [12.930118217431987, 13.930118217431989, 14.930118217431987],
            [ 9.092184727756681, 10.092184727756681, 11.092184727756681],
            [ 9.919698370675327, 10.919698370675325, 11.919698370675327],
            [11.454569130368588, 12.454569130368586, 13.454569130368588],
            [12.486120423462921, 13.486120423462921, 14.486120423462921],
            [10.595147944188897, 11.595147944188895, 12.595147944188895],
            [11.636465499406407, 12.63646549940641,  13.636465499406409],
            [13.171336259099668, 14.17133625909967,  15.171336259099666],
            [14.686329600585715, 15.686329600585713, 16.686329600585715]
        ])
        self.assertAlmostEqual(numpy.linalg.norm(filtered_data.Evaluate() - ref_values), 0.0, 12)

        filtered_data = vm_filter.FilterIntegratedField(unfiltered_field)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_data.Evaluate() - ref_values / 0.25), 0.0, 12)

    def test_ElementExplicitVertexMorphingFilter(self):
        radius = 1.0
        vm_filter = KratosOA.ElementExplicitVertexMorphingFilter(self.model_part, "linear", 1000)

        filter_radius = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, radius)
        vm_filter.SetFilterRadius(filter_radius)

        unfiltered_field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(unfiltered_field, Kratos.VELOCITY)
        filtered_data = vm_filter.FilterField(unfiltered_field)

        ref_values = numpy.array([
            [11.55444743868707,  12.554447438687067, 13.554447438687067],
            [11.918041958443546, 12.918041958443546, 13.918041958443544],
            [11.928193300444164, 12.928193300444164, 13.928193300444162],
            [12.282886864969507, 13.282886864969505, 14.282886864969505],
            [11.5471098313716,   12.547109831371602, 13.547109831371603],
            [11.246916707604543, 12.24691670760454,  13.24691670760454 ],
            [12.36654479019958,  13.36654479019958,  14.366544790199583],
            [10.896867880093462, 11.896867880093463, 12.896867880093463],
            [11.655383133574905, 12.655383133574908, 13.655383133574906],
            [11.455040010196624, 12.455040010196623, 13.45504001019662 ],
            [12.511738531311057, 13.511738531311059, 14.511738531311057],
            [12.231281395258035, 13.231281395258033, 14.231281395258035],
            [12.221777211365726, 13.221777211365728, 14.221777211365726],
            [13.052716422133567, 14.052716422133567, 15.052716422133567],
            [13.436252830633203, 14.436252830633203, 15.436252830633203],
            [13.168849831216496, 14.1688498312165,   15.168849831216498],
            [11.369765135635856, 12.369765135635854, 13.369765135635857],
            [12.375621221163348, 13.375621221163348, 14.375621221163346],
            [11.584368222275419, 12.584368222275417, 13.584368222275417],
            [12.41708255783569,  13.417082557835688, 14.41708255783569 ],
            [12.22618075312927,  13.226180753129269, 14.226180753129274],
            [12.88371497110807,  13.883714971108073, 14.883714971108072],
            [12.180052031010879, 13.180052031010879, 14.180052031010879],
            [13.302694861592688, 14.302694861592686, 15.302694861592688]
        ])

        self.assertAlmostEqual(numpy.linalg.norm(filtered_data.Evaluate() - ref_values), 0.0, 12)

        filtered_data = vm_filter.FilterIntegratedField(unfiltered_field)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_data.Evaluate() - ref_values * 24.0), 0.0, 10)


if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()