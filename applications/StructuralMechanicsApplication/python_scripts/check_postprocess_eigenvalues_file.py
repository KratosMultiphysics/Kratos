from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Other imports
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CheckPostprocessEigenvaluesFile(Model, settings["Parameters"])

class CheckPostprocessEigenvaluesFile(KratosMultiphysics.Process, KratosUnittest.TestCase):
    """
    This process compares the *.post.res file written by the PostprocessEigenvaluesProcess
    """
    def __init__(self, Model, settings):
        
        ## Settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "reference_file_name"   : "",
            "output_file_name"      : ""
        }
        """)
        
        ## Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
    def ExecuteInitialize(self):
        self.reference_file_name = os.path.join(os.getcwd(), self.settings["reference_file_name"].GetString())
        self.output_file_name = os.path.join(os.getcwd(), self.settings["output_file_name"].GetString())
        
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
        self._compare_post_res_file(self.reference_file_name, self.output_file_name)


    def _compare_post_res_file(self, ReferenceFileName, OutputFileName):
        ref_file = open(ReferenceFileName,'r')
        out_file = open(OutputFileName,'r')
        
        lines1 = ref_file.readlines()
        lines2 = out_file.readlines()
        num_lines_1 = len(lines1)

        if num_lines_1 != len(lines2):
            raise Exception("Files have different number of lines!")

        results_start_index = -1
        results_available = False

        # comparing the header
        for i in range(num_lines_1):
            if lines1[i].startswith("Result "):
                results_start_index = i
                results_available = True
                break

            tmp1 = lines1[i].split()
            tmp2 = lines2[i].split()

            if len(tmp1) != len(tmp2):
                raise Exception("Lines have different length!")

            for j in range(len(tmp1)):
                self.assertTrue(tmp1[j] == tmp2[j], msg=tmp1[j] + " != " + tmp2[j])

        # comparing the results
        if results_available:
            while results_start_index < num_lines_1:
                results_start_index = self._compare_results_block(lines1, lines2, results_start_index)
        
        ref_file.close()
        out_file.close()

    def _compare_results_block(self, lines1, lines2, current_index):
        # comparing result labels
        tmp1 = lines1[current_index].split()
        tmp2 = lines2[current_index].split()

        if len(tmp1) != len(tmp2):
            raise Exception("Result labels have different length!")

        for j in range(len(tmp1)):
            self.assertTrue(tmp1[j] == tmp2[j], msg=tmp1[j] + " != " + tmp2[j])

        current_index += 1 # skipping "Values"-line

        # comparing results
        while lines1[current_index+1] != "End Values\n":
            current_index += 1
            tmp1 = lines1[current_index].split()
            tmp2 = lines2[current_index].split()
            if len(tmp1) != len(tmp2):
                raise Exception("Different number of results!")

            for j in range(len(tmp1)):
                # comparing the values
                # abs bcs the sign does not matter for eigenvalues
                # tolerance is small bcs values in file are also small
                self.assertAlmostEqual(abs(float(tmp1[j])),
                                       abs(float(tmp2[j])), 15)

        return current_index+2 # directly incrementing to get the new result label