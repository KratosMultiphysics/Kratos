# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Check external dependencies
try:
  import KratosMultiphysics 
  import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
  missing_external_dependencies = False
  missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1)


# Tests for elements:

# Small displacement elements SDE
class PatchTest2D(TF.TestFactory):
    file_name = "contact_patch_test/ContactPatchTest_EP_2D"
    file_parameters = "contact_patch_test/ProjectParameters.json"
    

def SetTestSuite(suites):
    small_suite = suites['small']
    
    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            #SDE
            PatchTest2D,
        ])
    )

    return small_suite
