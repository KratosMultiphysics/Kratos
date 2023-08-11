import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication


def run():
    KM.Tester.SetVerbosity(KM.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    KM.Tester.RunTestCases("*LocalElasticAnisotropicVariables*")
    KM.Tester.RunTestCases("*LocalElasticAnisotropicDamage*")
    #KM.Tester.RunTestSuite("KratosConstitutiveLawsFastSuite")


if __name__ == '__main__':
    run()
