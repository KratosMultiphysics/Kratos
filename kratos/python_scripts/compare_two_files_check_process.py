from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics

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

    def __init__(self, model, params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "reference_file_name"   : "",
            "output_file_name"      : "",
            "remove_output_file"    : true,
            "comparison_type"       : "deterministic",
            "decimal_places"        : 6,
            "dimension"             : 3
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        params.ValidateAndAssignDefaults(default_parameters)
        self.reference_file_name = os.path.join(os.getcwd(), params["reference_file_name"].GetString())
        self.output_file_name = os.path.join(os.getcwd(), params["output_file_name"].GetString())
        self.remove_output_file = params["remove_output_file"].GetBool()
        self.comparison_type = params["comparison_type"].GetString()
        self.decimal_places = params["decimal_places"].GetInt()
        self.dimension = params["dimension"].GetInt()

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
        if kratos_utils.IsRankZero():
            if (self.comparison_type == "deterministic"):
                value = filecmp.cmp(self.reference_file_name, self.output_file_name)
                self.assertTrue(value)
            elif (self.comparison_type == "mesh_file"):
                self.__CompareMeshVertivesFile()
            elif (self.comparison_type == "sol_file"):
                self.__CompareSolMetricFile()
            elif (self.comparison_type == "post_res_file"):
                self.__ComparePostResFile()
            elif (self.comparison_type == "dat_file"):
                self.__CompareDatFile()
            else:
                raise NameError('Requested comparision type "' + self.comparison_type + '" not implemented yet')

        if self.remove_output_file == True:
            kratos_utils.DeleteFileIfExisting(self.output_file_name) # this checks internally if it is rank 0

    def __GetFileLines(self):
        with open(self.reference_file_name,'r') as ref_file:
            lines_ref = ref_file.readlines()
        with open(self.output_file_name,'r') as out_file:
            lines_out = out_file.readlines()

        num_lines_ref = len(lines_ref)
        num_lines_out = len(lines_out)

        err_msg  = "Files have different number of lines!"
        err_msg += "\nNum Lines Reference File: " + str(num_lines_ref)
        err_msg += "\nNum Lines Output File: " + str(num_lines_out)
        self.assertTrue(num_lines_ref == num_lines_out, msg=err_msg)

        return lines_ref, lines_out

    def __ComparePostResFile(self):
        lines_ref, lines_out = self.__GetFileLines()

        results_start_index = -1
        results_found = False

        num_lines = len(lines_ref)

        # comparing the header
        for i in range(num_lines):
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
            while results_start_index < num_lines:
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
                                       self.decimal_places)

        return current_index+2 # directly incrementing to get the new result label

    def __CompareDatFile(self):
        lines_ref, lines_out = self.__GetFileLines()

        # assert headers are the same
        while lines_ref[0].lstrip()[0] == '#' or lines_out[0].lstrip()[0] == '#':
            self.assertTrue(lines_ref.pop(0) == lines_out.pop(0))

        # assert values are equal up to given tolerance
        for line_ref, line_out in zip(lines_ref, lines_out):
            for v1, v2 in zip(line_ref.split(), line_out.split()):
                self.assertAlmostEqual(float(v1),
                                       float(v2),
                                       self.decimal_places)

    def __CompareMeshVertivesFile(self):
        lines_ref, lines_out = self.__GetFileLines()

        numline = 0
        for line1 in lines_ref:
            numline += 1

            if("Vertices" in line1):
                line = lines_ref[numline]
                nvertices = int(line)
                numline += 1
                break

        error = 0.0
        for i in range(numline, nvertices + numline):
            tmp1 = ConvertStringToListFloat(lines_ref[i], "", "\n")
            tmp2 = ConvertStringToListFloat(lines_out[i], "", "\n")
            if (dimension == 2):
                error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0)**(0.5)
            else:
                error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0 + (tmp1[2] - tmp2[2])**2.0)**(0.5)

        error /= nvertices
        self.assertTrue(error < GetTolerance(self.decimal_places))

    def __CompareSolMetricFile(self):
        lines_ref, lines_out = self.__GetFileLines()

        numline = 0
        for line1 in lines_ref:
            numline += 1

            if("SolAtVertices" in line1):
                line = lines_ref[numline]
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

            if (lines_ref[i][0] == " "):
                lines_ref[i] = lines_ref[i][1:]
            if (lines_out[i][0] == " "):
                lines_ref[i][0] = lines_out[i][1:]
            tmp1 = ConvertStringToListFloat(lines_ref[i], space, end_line)
            tmp2 = ConvertStringToListFloat(lines_out[i], space, end_line)

            if (dimension == 2):
                error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0 + (tmp1[2] - tmp2[2])**2.0)**(0.5)
            else:
                error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0 + (tmp1[2] - tmp2[2])**2.0 + (tmp1[3] - tmp2[3])**2.0 + (tmp1[4] - tmp2[4])**2.0 + (tmp1[5] - tmp2[5])**2.0)**(0.5)

        error /= nvertices
        self.assertTrue(error < GetTolerance(self.decimal_places))


def ConvertStringToListFloat(line, space = " ", endline = ""):
    list_values = []
    string_values = (line.replace(endline,"")).split(space)
    for string in string_values:
        list_values.append(float(string))

    return list_values

def GetTolerance(decimal_places):
    # convert e.g. 5 to 1e-5
    tolerance = 0.1**decimal_places
    tolerance = round(tolerance, decimal_places) # to remove rounding errors
    return(tolerance)
