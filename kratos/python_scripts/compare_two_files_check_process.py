from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

import filecmp
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CompareTwoFilesCheckProcess(Model, settings["Parameters"])

class CompareTwoFilesCheckProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):

    def __init__(self, Model, params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "reference_file_name"   : "",
            "output_file_name"      : "",
            "remove_output_file"    : true,
            "comparison_type"       : "deterministic",
            "tolerance"             : 1.0e-6,
            "dimension"             : 3
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
        self.reference_file_name = os.path.join(os.getcwd(), self.params["reference_file_name"].GetString())
        self.output_file_name = os.path.join(os.getcwd(), self.params["output_file_name"].GetString())
        self.remove_output_file = self.params["remove_output_file"].GetBool()
        self.comparison_type = self.params["comparison_type"].GetString()
        self.tolerance = self.params["tolerance"].GetDouble()
        self.dimension = self.params["dimension"].GetInt()

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
        if (self.comparison_type == "deterministic"):
            value = filecmp.cmp(self.reference_file_name, self.output_file_name)
            self.assertTrue(value)
        elif (self.comparison_type == "mesh_file"):
            error = _ReadVertices(self.reference_file_name, self.output_file_name, self.dimension)
            self.assertTrue(error < self.tolerance)
        elif (self.comparison_type == "sol_file"):
            error = _ReadMetric(self.reference_file_name, self.output_file_name, self.dimension)
            self.assertTrue(error < self.tolerance)
        elif (self.comparison_type == "post_res_file"):
            self.__ComparePostResFile()
        else:
            raise NameError('Requested comparision type "' + self.comparison_type + '" not implemented yet')

        if self.remove_output_file == True:
            kratos_utils.DeleteFileIfExisting(self.output_file_name)


    def __ComparePostResFile(self):
        with open(self.reference_file_name,'r') as ref_file, open(self.output_file_name,'r') as out_file:
            lines_ref = ref_file.readlines()
            lines_out = out_file.readlines()
            num_lines_1 = len(lines_ref)

            if num_lines_1 != len(lines_out):
                self.assertTrue(False, msg="Files have different number of lines!")

            results_start_index = -1
            results_found = False

            # comparing the header
            for i in range(num_lines_1):
                if lines_ref[i].startswith("Result "):
                    results_start_index = i
                    results_found = True
                    break

                lines_ref_splitted = lines_ref[i].split()
                lines_out_splitted = lines_out[i].split()

                if len(lines_ref_splitted) != len(lines_out_splitted):
                    self.assertTrue(False, msg="Lines have different length!")

                for ref_value, out_value in zip(lines_ref_splitted, lines_out_splitted):
                    self.assertTrue(ref_value == out_value,
                                    msg=ref_value + " != " + out_value)

            # comparing the results
            if results_found:
                while results_start_index < num_lines_1:
                    results_start_index = self.__CompareResultsBlock(lines_ref, lines_out, results_start_index)


    def __CompareResultsBlock(self, lines1, lines2, current_index):
        # comparing result labels
        lines_1_splitted = lines1[current_index].split()
        lines_2_splitted = lines2[current_index].split()

        if len(lines_1_splitted) != len(lines_2_splitted):
            self.assertTrue(False, msg="Result labels have different length!")

        for val_1, val_2 in zip(lines_1_splitted, lines_2_splitted):
            self.assertTrue(val_1 == val_2,
                            msg=val_1 + " != " + val_2)

        current_index += 1 # skipping "Values"-line

        # comparing results
        while lines1[current_index+1] != "End Values\n":
            current_index += 1
            lines_1_splitted = lines1[current_index].split()
            lines_2_splitted = lines2[current_index].split()
            if len(lines_1_splitted) != len(lines_2_splitted):
                self.assertTrue(False, msg="Different number of results!")

            for val_1, val_2 in zip(lines_1_splitted, lines_2_splitted):
                self.assertAlmostEqual(float(val_1),
                                       float(val_2),
                                       self.tolerance)

        return current_index+2 # directly incrementing to get the new result label

def _ConvertStringToListFloat(line, space = " ", endline = ""):
    list_values = []
    string_values = (line.replace(endline,"")).split(space)
    for string in string_values:
        list_values.append(float(string))
    #value = ""
    #for i in line:
        #if ((i == " ") or (i == "  ")):
            #num_value = float(value)
            #list_values.append(num_value)
            #value = ""
        #else:
            #value += i
    return list_values

def _ReadVertices(input_file1, input_file2, dimension):
    with open(input_file1,'r') as f1, open(input_file2,'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

        numline = 0
        for line1 in lines1:
            numline += 1

            if("Vertices" in line1):
                line = lines1[numline]
                nvertices = int(line)
                numline += 1
                break

        error = 0.0
        for i in range(numline, nvertices + numline):
            tmp1 = _ConvertStringToListFloat(lines1[i], "", "\n")
            tmp2 = _ConvertStringToListFloat(lines2[i], "", "\n")
            if (dimension == 2):
                error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0)**(0.5)
            else:
                error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0 + (tmp1[2] - tmp2[2])**2.0)**(0.5)

    return (error/nvertices)

def _ReadMetric(input_file1, input_file2, dimension):
    with open(input_file1,'r') as f1, open(input_file2,'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

        numline = 0
        for line1 in lines1:
            numline += 1

            if("SolAtVertices" in line1):
                line = lines1[numline]
                nvertices = int(line)
                numline += 2
                break

        error = 0.0
        for i in range(numline, nvertices + numline):

            if dimension == 2:
                space = " "
                end_line = " \n"
            else:
                space = "  "
                end_line = "  \n"

            if (lines1[i][0] == " "):
                lines1[i] = lines1[i][1:]
            if (lines2[i][0] == " "):
                lines1[i][0] = lines2[i][1:]
            tmp1 = _ConvertStringToListFloat(lines1[i], space, end_line)
            tmp2 = _ConvertStringToListFloat(lines2[i], space, end_line)

            if (dimension == 2):
                error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0 + (tmp1[2] - tmp2[2])**2.0)**(0.5)
            else:
                error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0 + (tmp1[2] - tmp2[2])**2.0 + (tmp1[3] - tmp2[3])**2.0 + (tmp1[4] - tmp2[4])**2.0 + (tmp1[5] - tmp2[5])**2.0)**(0.5)

    return (error/nvertices)
