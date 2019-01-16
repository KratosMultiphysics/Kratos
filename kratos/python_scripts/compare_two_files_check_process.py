from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.KratosUnittest import isclose as t_isclose

# Other imports
import filecmp
import os
import math

def Factory(settings, current_model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CompareTwoFilesCheckProcess(settings["Parameters"])

class CompareTwoFilesCheckProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    """This process compares files that are written during a simulation
    against reference files.
    Please see the "ExecuteFinalize" functions for details about the
    available file-formats
    """
    def __init__(self, params):
        KratosMultiphysics.Process.__init__(self)
        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "help"                  : "This process checks that two files are the same. This can be used in order to create tests, where a given solution is expected",
            "reference_file_name"   : "",
            "output_file_name"      : "",
            "remove_output_file"    : true,
            "comparison_type"       : "deterministic",
            "tolerance"             : 1e-6,
            "relative_tolerance"    : 1e-9,
            "dimension"             : 3
        }
        """)

        # backwards compatibility
        if params.Has("decimal_places"):
            if params.Has("tolerance") or params.Has("relative_tolerance"):
                raise Exception('Conflicting settings specified, please remove "decimal_places"')
            decimal_places = params["decimal_places"].GetInt()
            params.RemoveValue("decimal_places")
            warning =  'W-A-R-N-I-N-G: You have specified "decimal_places", '
            warning += 'which is deprecated and will be removed soon.\n'
            warning += 'Please specify "tolerance" and "relative_tolerance" instead!'
            KratosMultiphysics.Logger.PrintWarning("CompareTwoFilesCheckProcess", warning)
            tol = 0.1**decimal_places
            params.AddEmptyValue("tolerance").SetDouble(tol)

        ## Overwrite the default settings with user-provided parameters
        params.ValidateAndAssignDefaults(default_parameters)

        # abspath to make paths os-independent
        ref_file_name = os.path.abspath(params["reference_file_name"].GetString())
        out_file_name = os.path.abspath(params["output_file_name"].GetString())

        self.reference_file_name = os.path.join(os.getcwd(), ref_file_name)
        self.output_file_name = os.path.join(os.getcwd(), out_file_name)

        self.remove_output_file = params["remove_output_file"].GetBool()
        self.comparison_type = params["comparison_type"].GetString()
        self.tol = params["tolerance"].GetDouble()
        self.reltol = params["relative_tolerance"].GetDouble()
        self.dimension = params["dimension"].GetInt()

        self.info_msg = "".join([  "\n[%s]: Failed with following parameters:\n" % self.__class__.__name__,
                                    params.PrettyPrintJsonString()
                                ])

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
        """The files are compared in this function
        Please see the respective files for details on the format of the files
        """
        if (self.comparison_type == "deterministic"):
            value = filecmp.cmp(self.reference_file_name, self.output_file_name)
            self.assertTrue(value, msg = self.info_msg)
        elif (self.comparison_type == "mesh_file"):
            self.__CompareMeshVerticesFile()
        elif (self.comparison_type == "sol_file"):
            self.__CompareSolMetricFile()
        elif (self.comparison_type == "post_res_file"):
            self.__ComparePostResFile()
        elif (self.comparison_type == "dat_file"):
            self.__CompareDatFile()
        elif (self.comparison_type == "dat_file_variables_time_history"):
            self.__CompareDatFileVariablesTimeHistory()
        else:
            raise NameError('Requested comparision type "' + self.comparison_type + '" not implemented yet')

        if self.remove_output_file == True:
            kratos_utils.DeleteFileIfExisting(self.output_file_name)

    def __GetFileLines(self):
        """This function reads the reference and the output file
        It returns the lines read from both files and also compares
        if they contain the same numer of lines
        """
        # check if files are valid
        if not os.path.isfile(self.reference_file_name):
            err_msg  = 'The specified reference file name "'
            err_msg += self.reference_file_name
            err_msg += '" is not valid!'
            raise Exception(err_msg)
        if not os.path.isfile(self.output_file_name):
            err_msg  = 'The specified output file name "'
            err_msg += self.output_file_name
            err_msg += '" is not valid!'
            raise Exception(err_msg)

        # "readlines" adds a newline at the end of the line,
        # which will be removed with rstrip afterwards
        with open(self.reference_file_name,'r') as ref_file:
            lines_ref = ref_file.readlines()
        with open(self.output_file_name,'r') as out_file:
            lines_out = out_file.readlines()

        # removing trailing newline AND whitespaces than can mess with the comparison
        lines_ref = [line.rstrip() for line in lines_ref]
        lines_out = [line.rstrip() for line in lines_out]

        num_lines_ref = len(lines_ref)
        num_lines_out = len(lines_out)

        err_msg  = "Files have different number of lines!"
        err_msg += "\nNum Lines Reference File: " + str(num_lines_ref)
        err_msg += "\nNum Lines Output File: " + str(num_lines_out)
        self.assertTrue(num_lines_ref == num_lines_out, msg=err_msg + self.info_msg)

        return lines_ref, lines_out

    def __ComparePostResFile(self):
        """Comparing output files from GiD containing results in ASCII format
        => *.post.res
        see https://www.gidhome.com/documents/customizationmanual/POSTPROCESS%20DATA%20FILES
        """

        lines_ref, lines_out = self.__GetFileLines()

        results_found, results_start_index = self.__ComparePostResFileHeader(lines_ref, lines_out)

        # comparing the results (if there are any)
        if results_found:
            while results_start_index < len(lines_ref):
                results_start_index = self.__ComparePostResFileResultsBlock(lines_ref, lines_out, results_start_index)

    def __ComparePostResFileHeader(self, lines_ref, lines_out):
        """This function compares the headers of *.post.res files
        It loops the lines until Results are found
        """
        results_start_index = -1
        results_found = False

        for i in range(len(lines_ref)):
            if lines_ref[i].startswith("Result "):
                # Now the end of the header is found
                results_start_index = i
                results_found = True
                break

            lines_ref_splitted = lines_ref[i].split()
            lines_out_splitted = lines_out[i].split()

            self.assertTrue(len(lines_ref_splitted) == len(lines_out_splitted),
                            msg="Lines have different length!" + self.info_msg)

            for ref_value, out_value in zip(lines_ref_splitted, lines_out_splitted):
                self.assertTrue(ref_value == out_value,
                                msg=ref_value + " != " + out_value + self.info_msg)

        return results_found, results_start_index

    def __ComparePostResFileResultsBlock(self, lines1, lines2, current_index):
        """This function compares results blocks of *.post.res files
        """
        # comparing result labels
        lines_1_splitted = lines1[current_index].split()
        lines_2_splitted = lines2[current_index].split()

        if len(lines_1_splitted) != len(lines_2_splitted):
            self.assertTrue(False, msg="Result labels have different length!" + self.info_msg)

        for val_1, val_2 in zip(lines_1_splitted, lines_2_splitted):
            self.assertTrue(val_1 == val_2,
                            msg=val_1 + " != " + val_2 + self.info_msg)

        current_index += 1 # skipping "Values"-line

        # comparing results
        while lines1[current_index+1] != "End Values":
            current_index += 1
            lines_1_splitted = lines1[current_index].split()
            lines_2_splitted = lines2[current_index].split()
            if len(lines_1_splitted) != len(lines_2_splitted):
                self.assertTrue(False, msg="Different number of results!" + self.info_msg)

            for val_1, val_2 in zip(lines_1_splitted, lines_2_splitted):
                self.__CheckCloseValues(float(val_1), float(val_2))

        return current_index+2 # directly incrementing to get the new result label

    def __CompareDatFile(self):
        """This function compares files with tabular data.
        => *.dat
        Lines starting with "#" are comments and therefore compared for equality
        Other lines are compared to be almost equal with the specified tolerance
        """
        lines_ref, lines_out = self.__GetFileLines()

        # assert headers are the same
        lines_ref, lines_out = self.__CompareDatFileComments(lines_ref, lines_out)

        # assert values are equal up to given tolerance
        self.__CompareDatFileResults(lines_ref, lines_out)

    def __CompareDatFileVariablesTimeHistory(self):
        """This function compares files with tabular data.
        => *.dat
        Lines starting with "#" are comments and therefore compared for equality
        Other lines are compared to be almost equal with the specified tolerance
        If the comparison fails, it prints the location of failure with details
        The expected format is the one written by the PointOutputProcess:

        # some basic file information
        # time var_name_1 var_name_2
        0.1 1.2345 2.852
        0.2 0.889 -89.444
        .
        .
        .
        """
        lines_ref, lines_out = self.__GetFileLines()

        # extracting the names of output variables eg: time, VELOCITY_X, VELOCITY_Y, VELOCITY_Z
        variable_names = lines_ref[1].split()[2:]

        # assert headers are the same
        lines_ref, lines_out = self.__CompareDatFileComments(lines_ref, lines_out)

        # assert values are equal up to given tolerance
        self.__CompareDatFileResultsWithLocation(lines_ref, lines_out, variable_names)

    def __CompareDatFileComments(self, lines_ref, lines_out):
        """This function compares the comments of files with tabular data
        The lines starting with "#" are being compared
        These lines are removed from the list of lines
        """
        for line_ref, line_out in zip(lines_ref, lines_out):
            if line_ref.lstrip()[0] == '#' or line_out.lstrip()[0] == '#':
                self.assertTrue(line_ref == line_out, msg = self.info_msg)

        lines_ref = [line for line in lines_ref if not(line.lstrip()[0] == '#')]
        lines_out = [line for line in lines_out if not(line.lstrip()[0] == '#')]

        return lines_ref, lines_out

    def __CompareDatFileResults(self, lines_ref, lines_out):
        """This function compares the data of files with tabular data
        The comment lines were removed beforehand
        """
        for line_ref, line_out in zip(lines_ref, lines_out):
            for v1, v2 in zip(line_ref.split(), line_out.split()):
                self.__CheckCloseValues(float(v1), float(v2))

    def __CompareDatFileResultsWithLocation(self, lines_ref, lines_out, variable_names):
        """This function compares the data of files with tabular data
        It also prints the exact location where data doesnt match each other
        The comment lines were removed beforehand
        """
        for line_ref, line_out in zip(lines_ref, lines_out):
            for i_var, (v1, v2) in enumerate(zip(line_ref.split(), line_out.split())):
                if i_var == 0: # comparing time:
                    additional_info = 'Different time found!'
                    self.__CheckCloseValues(float(v1), float(v2), additional_info)
                    current_time = v1
                else: # comparing variables
                    additional_info  = 'Failed for variable ' + variable_names[i_var-1]
                    additional_info += ' at time: ' + current_time
                    self.__CheckCloseValues(float(v1), float(v2), additional_info)

    def __CompareMeshVerticesFile(self):
        """This function compares the output of the MMG meshing library
        => *.mesh
        see https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-tutorials
        """
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
            if (self.dimension == 2):
                error += math.sqrt((tmp1[0] - tmp2[0])**2 + (tmp1[1] - tmp2[1])**2)
            else:
                error += math.sqrt((tmp1[0] - tmp2[0])**2 + (tmp1[1] - tmp2[1])**2 + (tmp1[2] - tmp2[2])**2)

        error /= nvertices
        self.assertTrue(error < self.tol, msg = self.info_msg)

    def __CompareSolMetricFile(self):
        """This function compares the output of the MMG meshing library
        => *.sol
        see https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options/mmg-remesher-option-sol
        """
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
            if self.dimension == 2:
                space = " "
                end_line = " \n"
            else:
                space = "  "
                end_line = "  \n"

            if (lines_ref[i][0] == " "):
                tmp1 = ConvertStringToListFloat(lines_ref[i][1:], space, end_line)
            else:
                tmp1 = ConvertStringToListFloat(lines_ref[i], space, end_line)
            if (lines_out[i][0] == " "):
                tmp2 = ConvertStringToListFloat(lines_out[i][1:], space, end_line)
            else:
                tmp2 = ConvertStringToListFloat(lines_out[i], space, end_line)

            if (self.dimension == 2):
                error += math.sqrt((tmp1[0] - tmp2[0])**2 + (tmp1[1] - tmp2[1])**2 + (tmp1[2] - tmp2[2])**2)
            else:
                error += math.sqrt((tmp1[0] - tmp2[0])**2 + (tmp1[1] - tmp2[1])**2 + (tmp1[2] - tmp2[2])**2 + (tmp1[3] - tmp2[3])**2 + (tmp1[4] - tmp2[4])**2 + (tmp1[5] - tmp2[5])**2)

        error /= nvertices
        self.assertTrue(error < self.tol, msg = self.info_msg)

    def __CheckCloseValues(self, val_a, val_b, additional_info=""):
        isclosethis = t_isclose(val_a, val_b, rel_tol=self.reltol, abs_tol=self.tol)
        full_msg =  self.info_msg + "\n"
        full_msg += str(val_a) + " != " + str(val_b) + ", rel_tol = "
        full_msg += str(self.reltol) + ", abs_tol = " + str(self.tol)
        if additional_info != "":
            full_msg += "\n" + additional_info
        self.assertTrue(isclosethis, msg=full_msg)


def ConvertStringToListFloat(line, space = " ", endline = ""):
    """This function converts a string into a list of floats
    """
    list_values = []
    string_values = (line.replace(endline,"")).split(space)
    for string in string_values:
        if (string != ""):
            list_values.append(float(string))

    return list_values
