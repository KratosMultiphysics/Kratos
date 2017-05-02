from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
KratosMultiphysics.CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CompareTwoPointOutputFilesMPIProcess(Model, settings["Parameters"])

class CompareTwoPointOutputFilesMPIProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):

    def __init__(self,model_part,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "file_name_1"            : "PLEASE_SPECIFY_FILENAME",
            "file_name_2"            : "PLEASE_SPECIFY_FILENAME",
            "decimal_places"         : 5
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
        self.file_name_1 = self.params["file_name_1"].GetString()
        self.file_name_2 = self.params["file_name_2"].GetString()
        self.decimal_places = self.params["decimal_places"].GetInt()

    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass
            
    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        KratosMPI.mpi.world.barrier()
        if KratosMPI.mpi.rank == 0:
            with open(self.file_name_1,'r') as f1:
                lines1 = f1.readlines()
            with open(self.file_name_2,'r') as f2:
                lines2 = f2.readlines()

            # assert headers are the same
            while lines1[0].lstrip()[0] == '#' or lines2[0].lstrip()[0] == '#':
                self.assertTrue(lines1.pop(0) == lines2.pop(0))

            # assert values are equal up to given tolerance
            for line1, line2 in zip(lines1, lines2):
                for v1, v2 in zip(line1.split(), line2.split()):
                    self.assertAlmostEqual(float(v1), float(v2), self.decimal_places)
