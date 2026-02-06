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
    if not isinstance(settings, KratosMultiphysics.Parameters):
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

    def Execute(self):
        """Executes all functions required to compare the files
        Intended to be directly used within python-scripts
        """
        self.ExecuteFinalize()

    def ExecuteFinalize(self):
        """The files are compared in this function
        Please see the respective files for details on the format of the files
        """

        KratosMultiphysics.Testing.GetDefaultDataCommunicator().Barrier()
        ## only do the check in ranks zero, otherwise this can experience in race condition
        if KratosMultiphysics.Testing.GetDefaultDataCommunicator().Rank() != 0:
            return

        if self.comparison_type == "deterministic":
            value = filecmp.cmp(self.reference_file_name, self.output_file_name)
            self.assertTrue(value, msg = self.info_msg)
        elif self.comparison_type == "mesh_file":
            self.__CompareMeshVerticesFile()
        elif self.comparison_type == "sol_file":
            self.__CompareSolMetricFile()
        elif self.comparison_type == "post_res_file":
            self.__ComparePostResFile()
        elif self.comparison_type == "dat_file":
            self.__CompareDatFile()
        elif self.comparison_type == "csv_file":
            self.__CompareCSVFile()
        elif self.comparison_type == "dat_file_variables_time_history":
            self.__CompareDatFileVariablesTimeHistory()
        elif self.comparison_type == "vtk":
            self.__CompareVtkFile()
        elif self.comparison_type == "case":
            self.__CompareCaseFile()
        elif self.comparison_type == "geo":
            self.__CompareGeoFile()
        elif self.comparison_type == "ensight_solution":
            self.__CompareEnsightSolutionFile()
        else:
            raise NameError('Requested comparison type "' + self.comparison_type + '" not implemented yet')

        if self.remove_output_file:
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
        self.__CompareDelimittedFileResults(lines_ref, lines_out, None)

    def __CompareCSVFile(self):
        """This function compares files with tabular data.
        => *.csv
        Lines starting with "#" are comments and therefore compared for equality
        Other lines are compared to be almost equal with the specified tolerance
        """
        lines_ref, lines_out = self.__GetFileLines()

        # assert headers are the same
        lines_ref, lines_out = self.__CompareDatFileComments(lines_ref, lines_out)

        # assert values are equal up to given tolerance
        self.__CompareDelimittedFileResults(lines_ref, lines_out, ",")

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

    def __CompareDelimittedFileResults(self, lines_ref, lines_out, delimiter):
        """This function compares the data of files with tabular data
        The comment lines were removed beforehand
        """
        for line_ref, line_out in zip(lines_ref, lines_out):
            for v1, v2 in zip(line_ref.split(delimiter), line_out.split(delimiter)):
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

            if "Vertices" in line1:
                line = lines_ref[numline]
                nvertices = int(line)
                numline += 1
                break

        error = 0.0
        for i in range(numline, nvertices + numline):
            tmp1 = ConvertStringToListFloat(lines_ref[i], "", "\n")
            tmp2 = ConvertStringToListFloat(lines_out[i], "", "\n")
            if self.dimension == 2:
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

            if "SolAtVertices" in line1:
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
                space = " "
                end_line = "  \n"

            if lines_ref[i][0] == " ":
                tmp1 = ConvertStringToListFloat(lines_ref[i][1:], space, end_line)
            else:
                tmp1 = ConvertStringToListFloat(lines_ref[i], space, end_line)
            if lines_out[i][0] == " ":
                tmp2 = ConvertStringToListFloat(lines_out[i][1:], space, end_line)
            else:
                tmp2 = ConvertStringToListFloat(lines_out[i], space, end_line)

            if len(tmp1) == 1:
                error += tmp1[0] - tmp2[0]
            else:
                if self.dimension == 2:
                    error += math.sqrt((tmp1[0] - tmp2[0])**2 + (tmp1[1] - tmp2[1])**2 + (tmp1[2] - tmp2[2])**2)
                else:
                    error += math.sqrt((tmp1[0] - tmp2[0])**2 + (tmp1[1] - tmp2[1])**2 + (tmp1[2] - tmp2[2])**2 + (tmp1[3] - tmp2[3])**2 + (tmp1[4] - tmp2[4])**2 + (tmp1[5] - tmp2[5])**2)

        error /= nvertices
        self.assertTrue(error < self.tol, msg = self.info_msg)

    def __CompareVtkFile(self):
        """This function compares vtk files in ASCII format
        """

        def ReadVectorFromLine(line, conversion_fct, delimiter=' '):
            return [conversion_fct(word) for word in line.split(delimiter)]

        def CheckHeader(first_lines_file):
            # Expected header:
            '''
            # vtk DataFile Version 4.0
            vtk output
            ASCII
            DATASET UNSTRUCTURED_GRID
            '''
            if first_lines_file[2] != "ASCII":
                raise Exception("Only acsii files can be compared!")
            if first_lines_file[3] != "DATASET UNSTRUCTURED_GRID":
                raise Exception("unknown dataset")

        def ComparePoints(lines_ref, lines_out, line_counter):
            if not lines_out[line_counter].startswith("POINTS"):
                raise Exception('output-file is missing "POINTS" in the same location (expected line: {})'.format(line_counter))

            num_points_ref = int(lines_ref[line_counter].split(" ")[1])
            num_points_out = int(lines_out[line_counter].split(" ")[1])
            if not num_points_ref == num_points_out:
                raise Exception('output-file has wrong number of points: ref: {}, out: {}'.format(num_points_ref, num_points_out))

            # Comparing point coordinates
            for i in range(1, num_points_ref+1):
                ref_coords = ReadVectorFromLine(lines_ref[line_counter+i], float)
                out_coords = ReadVectorFromLine(lines_out[line_counter+i], float)

                for val_1, val_2 in zip(ref_coords, out_coords):
                    self.__CheckCloseValues(val_1, val_2)

        def CompareCells(lines_ref, lines_out, line_counter):
            if not lines_out[line_counter].startswith("CELLS"):
                raise Exception('output-file is missing "CELLS" in the same location (expected line: {})'.format(line_counter))

            num_cells_ref = int(lines_ref[line_counter].split(" ")[1])
            num_cells_out = int(lines_out[line_counter].split(" ")[1])
            if not num_cells_ref == num_cells_out:
                raise Exception('output-file has wrong number of cells: ref: {}, out: {}'.format(num_cells_ref, num_cells_out))

            num_connectivities_ref = int(lines_ref[line_counter].split(" ")[2])
            num_connectivities_out = int(lines_out[line_counter].split(" ")[2])
            if not num_connectivities_ref == num_connectivities_out:
                raise Exception('output-file has wrong number of connectivities: ref: {}, out: {}'.format(num_connectivities_ref, num_connectivities_out))

            # Comparing connectivities
            for i in range(1, num_cells_ref+1):
                ref_connectivities = ReadVectorFromLine(lines_ref[line_counter+i], int)
                out_connectivities = ReadVectorFromLine(lines_out[line_counter+i], int)

                for val_1, val_2 in zip(ref_connectivities, out_connectivities):
                    self.assertTrue(val_1 == val_2, msg='Wrong connectivity in line {}: ref: {}, out: {}'.format(line_counter+i+1, val_1, val_2))

        def CompareCellTypes(lines_ref, lines_out, line_counter):
            if not lines_out[line_counter].startswith("CELL_TYPES"):
                raise Exception('output-file is missing "CELL_TYPES" in the same location (expected line: {})'.format(line_counter))

            num_cells_ref = int(lines_ref[line_counter].split(" ")[1])
            num_cells_out = int(lines_out[line_counter].split(" ")[1])
            if not num_cells_ref == num_cells_out:
                raise Exception('output-file has wrong number of cells: ref: {}, out: {}'.format(num_cells_ref, num_cells_out))

            # Comparing cell types
            for i in range(1, num_cells_ref+1):
                ref_cell_type = ReadVectorFromLine(lines_ref[line_counter+i], int)
                out_cell_type = ReadVectorFromLine(lines_out[line_counter+i], int)
                self.assertTrue(ref_cell_type[0] == out_cell_type[0], msg='Wrong cell type in line {}: ref: {}, out: {}'.format(line_counter+i+1, ref_cell_type[0], out_cell_type[0]))

        def CompareData(lines_ref, lines_out, line_counter):
            def CompareFieldData(lines_ref, lines_out, line_counter, num_entities):
                ref_line_splitted = lines_ref[line_counter].split(" ")
                name_ref = ref_line_splitted[0]
                dim_ref = int(ref_line_splitted[1])
                num_entities_ref = int(ref_line_splitted[2])

                out_line_splitted = lines_out[line_counter].split(" ")
                name_out = out_line_splitted[0]
                dim_out = int(out_line_splitted[1])
                num_entities_out = int(out_line_splitted[2])

                if not num_entities_ref == num_entities:
                    raise Exception('Num entities is wrong in ref: expected: {}, got: {}'.format(num_entities_ref, num_entities))
                if not num_entities_out == num_entities:
                    raise Exception('Num entities is wrong in out: expected: {}, got: {}'.format(num_entities_out, num_entities))

                if not name_ref == name_out:
                    raise Exception('name of field is not matching: ref: {}, out: {}'.format(name_ref, name_out))

                if not dim_ref == dim_out:
                    raise Exception('dimension of field is not matching: ref: {}, out: {}'.format(dim_ref, dim_out))

                # Compare the values
                for i in range(1, num_entities+1):
                    ref_vals = ReadVectorFromLine(lines_ref[line_counter+i], float)
                    out_vals = ReadVectorFromLine(lines_out[line_counter+i], float)

                    for val_1, val_2 in zip(ref_vals, out_vals):
                        self.__CheckCloseValues(val_1, val_2)


            # Check if POINT_DATA or CELL_DATA
            data_type_ref = lines_ref[line_counter].split(" ")[0]
            data_type_out = lines_out[line_counter].split(" ")[0]
            if data_type_ref != data_type_out:
                raise Exception('data type is not matching: ref: {}, out: {}'.format(data_type_ref, data_type_out))

            num_entities_ref = int(lines_ref[line_counter].split(" ")[1])
            num_entities_out = int(lines_out[line_counter].split(" ")[1])
            if not num_entities_ref == num_entities_out:
                raise Exception('output-file has wrong number of entities for: {}: ref: {}, out: {}'.format(data_type_ref, num_entities_ref, num_entities_out))

            if not lines_ref[line_counter+1].startswith("FIELD"):
                raise Exception('reference-file is missing "FIELD" (expected in line: {})'.format(line_counter))
            if not lines_out[line_counter+1].startswith("FIELD"):
                raise Exception('output-file is missing "FIELD" (expected in line: {})'.format(line_counter))

            num_fields_ref = int(lines_ref[line_counter+1].split(" ")[2])
            num_fields_out = int(lines_out[line_counter+1].split(" ")[2])
            if not num_fields_ref == num_fields_out:
                raise Exception('output-file has wrong number of fields for: {}: ref: {}, out: {}'.format(data_type_ref, num_fields_ref, num_fields_out))

            for i in range(num_fields_ref):
                CompareFieldData(lines_ref, lines_out, line_counter+2 + i*(num_entities_ref+1), num_entities_ref)

        ######
        lines_ref, lines_out = self.__GetFileLines()

        CheckHeader(lines_ref[0:4])
        CheckHeader(lines_out[0:4])

        for line_counter, line_ref in enumerate(lines_ref):
            if line_ref.startswith("POINTS"):
                ComparePoints(lines_ref, lines_out, line_counter)
            if line_ref.startswith("CELLS"):
                CompareCells(lines_ref, lines_out, line_counter)
            if line_ref.startswith("CELL_TYPES"):
                CompareCellTypes(lines_ref, lines_out, line_counter)
            if line_ref.startswith("POINT_DATA") or line_ref.startswith("CELL_DATA"):
                CompareData(lines_ref, lines_out, line_counter)

    def __CompareCaseFile(self):
        """
        This function compares two EnSight Case (.case) files.
        It checks the format, geometry, time, and variable definitions.
        """

        def ParseCaseFile(lines):
            """
            Parses the lines of a .case file into a structured dictionary.
            """
            data = {
                "VARIABLE": [] # Use a list to store multiple variable definitions
            }
            current_section = None

            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                # Check for section headers
                if line in ["FORMAT", "GEOMETRY", "VARIABLE", "TIME"]:
                    current_section = line
                    continue

                # Parse content based on the current section
                if current_section == "FORMAT":
                    key, val = line.split(':')
                    data.setdefault("FORMAT", {})[key.strip()] = val.strip()

                elif current_section == "GEOMETRY":
                    key, val = line.split(':', 1) # Split only on the first colon
                    data.setdefault("GEOMETRY", {})[key.strip()] = val.strip()

                elif current_section == "TIME":
                    key, val = line.split(':', 1)
                    key = key.strip()
                    val = val.strip()
                    if key == "time values":
                        # Parse list of floats
                        parsed_val = [float(v) for v in val.split()]
                    else:
                        # Parse integers
                        parsed_val = int(val)
                    data.setdefault("TIME", {})[key] = parsed_val

                elif current_section == "VARIABLE":
                    # E.g.: "scalar per node: 1 PRESSURE Main.PRESSURE.****.node.scl"
                    parts = line.split(':')
                    var_type = parts[0].strip()

                    # Further split the remaining part by whitespace
                    details = parts[1].strip().split()
                    var_count = int(details[0])
                    var_name = details[1]
                    var_file = details[2]

                    variable_info = {
                        "description": var_type,
                        "count": var_count,
                        "name": var_name,
                        "file": var_file
                    }
                    data["VARIABLE"].append(variable_info)

            # Sort variables by name to make comparison order-independent
            data["VARIABLE"].sort(key=lambda v: v['name'])
            return data

        ###### Main execution starts here ######
        lines_ref, lines_out = self.__GetFileLines()

        # Parse both files into structured dictionaries
        data_ref = ParseCaseFile(lines_ref)
        data_out = ParseCaseFile(lines_out)

        # List to collect all error messages
        error_messages = []

        # 1. Compare FORMAT section
        format_ref = data_ref.get("FORMAT", {})
        format_out = data_out.get("FORMAT", {})
        if format_ref != format_out:
            error_messages.append(
                f"FORMAT sections do not match.\n  Ref: {format_ref}\n  Out: {format_out}"
            )

        # 2. Compare GEOMETRY section
        geometry_ref = data_ref.get("GEOMETRY", {})
        geometry_out = data_out.get("GEOMETRY", {})
        if geometry_ref != geometry_out:
            error_messages.append(
                f"GEOMETRY sections do not match.\n  Ref: {geometry_ref}\n  Out: {geometry_out}"
            )

        # 3. Compare TIME section
        # Use copies to safely modify them by popping "time values"
        time_ref = data_ref.get("TIME", {}).copy()
        time_out = data_out.get("TIME", {}).copy()

        if set(time_ref.keys()) != set(time_out.keys()):
            error_messages.append(
                f"TIME section keys do not match.\n  Ref keys: {set(time_ref.keys())}\n  Out keys: {set(time_out.keys())}"
            )
        else:
            # Handle the special case of "time values", which is a list of floats
            if "time values" in time_ref:
                ref_time_values = time_ref.pop("time values")
                out_time_values = time_out.pop("time values")

                if len(ref_time_values) != len(out_time_values):
                    error_messages.append(
                        f"Number of 'time values' entries does not match. Ref: {len(ref_time_values)}, Out: {len(out_time_values)}"
                    )
                else:
                    # Compare each float value with tolerance
                    for i, (ref_val, out_val) in enumerate(zip(ref_time_values, out_time_values)):
                        try:
                            # Assuming __CheckCloseValues raises an AssertionError on failure
                            self.__CheckCloseValues(ref_val, out_val)
                        except AssertionError as e:
                            error_messages.append(
                                f"Mismatch in 'time values' at index {i}: Ref={ref_val}, Out={out_val}. Details: {e}"
                            )

            # Compare all other remaining key-value pairs (expected to be integers)
            for key, ref_val in time_ref.items():
                out_val = time_out[key]
                if int(ref_val) != int(out_val):
                    error_messages.append(
                        f"TIME section mismatch for key '{key}'. Ref: {ref_val}, Out: {out_val}"
                    )

        # 4. Compare VARIABLE section
        vars_ref = data_ref.get("VARIABLE", [])
        vars_out = data_out.get("VARIABLE", [])

        if len(vars_ref) != len(vars_out):
            error_messages.append(
                f"Number of defined VARIABLEs does not match. Ref: {len(vars_ref)}, Out: {len(vars_out)}"
            )
        else:
            for var_ref, var_out in zip(vars_ref, vars_out):
                if var_ref != var_out:
                    var_name = var_ref.get('name', 'UNKNOWN')
                    error_messages.append(
                        f"Mismatch in VARIABLE definition for '{var_name}'.\n  Ref: {var_ref}\n  Out: {var_out}"
                    )

        # A single assertion at the end. It passes if the error_messages list is empty.
        if error_messages:
            print("Case file comparison failed:\n\n" + "\n\n".join(error_messages))
        self.assertTrue(not error_messages)

    def __CompareGeoFile(self):
        """
        This function compares two EnSight 6 Geometry (.geo) files in ASCII format.
        It checks nodes, coordinates, parts, and element connectivities.
        """
        def ReadVectorFromLine(line, conversion_fct):
            """
            Reads a vector of values from a single line string, 
            using any whitespace (spaces, tabs, etc.) as a delimiter.
            """
            # Use split() without arguments to handle any amount of whitespace as a delimiter
            # and automatically discard empty strings.
            return [conversion_fct(word) for word in line.split()]

        def CompareHeader(lines_ref, lines_out):
            """
            Checks the header section of the geo file. It expects the first 4 lines
            to define the file type and node/element id settings.
            """
            header_lines_to_check = [
                "EnSight 6 Geometry File",
                "Written by Kratos Multi-Physics",
                "node id given",
                "element id given"
            ]
            for i, line_content in enumerate(header_lines_to_check):
                if line_content not in lines_ref[i] or line_content not in lines_out[i]:
                    raise Exception(f'Header mismatch in line {i+1}. Expected "{line_content}"')

            # We can start parsing from after the header
            return len(header_lines_to_check)

        def CompareCoordinates(lines_ref, lines_out, line_counter, error_messages):
            """
            Compares the 'coordinates' block of the two files.
            """
            # Verify the "coordinates" keyword
            if "coordinates" not in lines_ref[line_counter]:
                raise Exception(f'Reference file is missing "coordinates" keyword at line {line_counter+1}')
            if "coordinates" not in lines_out[line_counter]:
                raise Exception(f'Output file is missing "coordinates" keyword at line {line_counter+1}')
            line_counter += 1

            # Compare number of nodes
            num_nodes_ref = int(lines_ref[line_counter])
            num_nodes_out = int(lines_out[line_counter])
            if num_nodes_ref != num_nodes_out:
                error_messages.append(f'Mismatch in number of nodes. Ref: {num_nodes_ref}, Out: {num_nodes_out}')
                # If counts differ, we cannot safely compare them. Skip to the end of this block.
                return line_counter + 1 + num_nodes_ref
            line_counter += 1

            # Compare node IDs and coordinates for each node
            for i in range(num_nodes_ref):
                current_line_idx = line_counter + i
                ref_node_data = ReadVectorFromLine(lines_ref[current_line_idx], float)
                out_node_data = ReadVectorFromLine(lines_out[current_line_idx], float)

                # Node ID (integer)
                if int(ref_node_data[0]) != int(out_node_data[0]):
                    error_messages.append(f'Node ID mismatch in line {current_line_idx+1}. Ref: {int(ref_node_data[0])}, Out: {int(out_node_data[0])}')

                # Coordinates (floats)
                for i_coord, (ref_coord, out_coord) in enumerate(zip(ref_node_data[1:], out_node_data[1:])):
                    try:
                        self.__CheckCloseValues(ref_coord, out_coord)
                    except AssertionError as e:
                        error_messages.append(f'Coordinate {i_coord+1} mismatch in line {current_line_idx+1}. Ref: {ref_coord}, Out: {out_coord}. Details: {e}')

            return line_counter + num_nodes_ref

        def ComparePart(lines_ref, lines_out, line_counter, error_messages):
            """
            Compares a single 'part' block, including its elements and connectivities.
            """
            list_of_element_types = [
                "point", "bar2", "bar3", "tria3", "tria6", "quad4", "quad8",
                "tetra4", "tetra10", "pyramid5", "pyramid13", "penta6",
                "penta15", "hexa8", "hexa20"
            ]
            if not lines_ref[line_counter].strip().startswith("part"):
                raise Exception(f'Reference file is missing "part" keyword at line {line_counter+1}')
            if not lines_out[line_counter].strip().startswith("part"):
                raise Exception(f'Output file is missing "part" keyword at line {line_counter+1}')
            line_counter += 1

            # Compare part number and description
            if lines_ref[line_counter].strip() != lines_out[line_counter].strip():
                error_messages.append(f'Part number mismatch at line {line_counter+1}.')
            line_counter += 1

            # Iterate through all element types in this part
            while line_counter < len(lines_ref) and lines_ref[line_counter] in list_of_element_types:
                # Compare element type
                elem_type_ref = lines_ref[line_counter].strip()
                elem_type_out = lines_out[line_counter].strip()
                if elem_type_ref != elem_type_out:
                    error_messages.append(f'Element type mismatch for part at line {line_counter-2}. Ref: "{elem_type_ref}", Out: "{elem_type_out}"')
                line_counter += 1

                # Compare number of elements
                num_elements_ref = int(lines_ref[line_counter])
                num_elements_out = int(lines_out[line_counter])
                if num_elements_ref != num_elements_out:
                    error_messages.append(f'Mismatch in number of elements for part "{elem_type_ref}". Ref: {num_elements_ref}, Out: {num_elements_out}')
                    # Skip to end of this part block based on reference file
                    return line_counter + 1 + num_elements_ref
                line_counter += 1

                # Compare element IDs and connectivities
                for i in range(num_elements_ref):
                    current_line_idx = line_counter + i
                    ref_connectivity = ReadVectorFromLine(lines_ref[current_line_idx], int)
                    out_connectivity = ReadVectorFromLine(lines_out[current_line_idx], int)
                    if ref_connectivity != out_connectivity:
                        error_messages.append(f'Element connectivity mismatch in line {current_line_idx+1}.\n  Ref: {ref_connectivity}\n  Out: {out_connectivity}')
                line_counter += num_elements_ref

            return line_counter # Return the updated line counter

        ###### Main execution starts here ######
        lines_ref, lines_out = self.__GetFileLines()

        # Pre-process by removing leading/trailing whitespace from all lines
        lines_ref = [line.strip() for line in lines_ref]
        lines_out = [line.strip() for line in lines_out]

        error_messages = []

        try:
            # Compare headers and get starting index for data
            line_counter = CompareHeader(lines_ref, lines_out)

            # The next block must be coordinates
            line_counter = CompareCoordinates(lines_ref, lines_out, line_counter, error_messages)

            # Process all subsequent 'part' blocks until the end of the file
            while line_counter < len(lines_ref):
                # If it's not an empty line, it must be a part
                line_counter = ComparePart(lines_ref, lines_out, line_counter, error_messages)

        except IndexError:
            error_messages.append("Comparison failed because files have different structures or lengths, leading to an out-of-bounds error.")
        except Exception as e:
            error_messages.append(f"A critical error occurred during comparison: {e}")

        # Final assertion on the collected results
        if error_messages:
            print("Geo file comparison failed:\n\n" + "\n\n".join(error_messages))
        self.assertTrue(not error_messages)

    def __CompareEnsightSolutionFile(self):
        """
        This function compares two EnSight 6 Solution files (scalar or vector).
        It checks the header for equality and then compares all numerical values
        using a hardcoded tolerance, while preserving the class method structure.
        Supported formats: .scl, .vec
        """
        # A hardcoded tolerance is defined for the check, ignoring class parameters.
        HARDCODED_TOLERANCE = 1e-9

        # Use the existing class method to read file lines.
        lines_ref, lines_out = self.__GetFileLines()

        # 1. Perform a direct comparison of the header line.
        header_ref = lines_ref[0]
        header_out = lines_out[0]
        if header_ref != header_out:
            # Explicitly fail the test if the headers do not match.
            error_msg = "EnSight solution file headers do not match.\n"
            error_msg += f"  Reference: '{header_ref.strip()}'\n"
            error_msg += f"  Output:    '{header_out.strip()}'"
            self.fail(error_msg + self.info_msg)

        number_of_lines = len(lines_ref)
        line_counter = 1 # Start after the header line

        # 2. Iterate through the numerical data lines (skipping the header and parts).
        while line_counter < number_of_lines:
            line_ref = lines_ref[line_counter]
            line_out = lines_out[line_counter]

            # Skip "part" lines if present (only in element and condition file types).
            if not "Per_node" in header_ref and "part" in line_ref:
                # Check it is the same part number, split and the number of the part is the second word
                part_number_ref = int(line_ref.split()[1])
                part_number_out = int(line_out.split()[1])
                if part_number_ref != part_number_out:
                    error_msg = f"Part numbers do not match at line {line_counter+1}.\n"
                    error_msg += f"  Reference: {part_number_ref}\n"
                    error_msg += f"  Output:    {part_number_out}"
                    self.fail(error_msg + self.info_msg)

                # Check the same type of element follows
                line_counter += 1
                line_ref = lines_ref[line_counter]
                line_out = lines_out[line_counter]
                element_type_ref = line_ref.strip()
                element_type_out = line_out.strip()
                if element_type_ref != element_type_out:
                    error_msg = f"Element types do not match after part {part_number_ref} at line {line_counter+1}.\n"
                    error_msg += f"  Reference: '{element_type_ref}'\n"
                    error_msg += f"  Output:    '{element_type_out}'"
                    self.fail(error_msg + self.info_msg)

                # Update the counter to move to the next line after "part" and element type
                line_counter += 1
                line_ref = lines_ref[line_counter]
                line_out = lines_out[line_counter]

            # Split lines by any whitespace to get numerical values as strings.
            values_ref_str = line_ref.split()
            values_out_str = line_out.split()

            # Ensure both lines contain the same number of values using a direct comparison.
            len_ref = len(values_ref_str)
            len_out = len(values_out_str)
            error_msg = f"Line {line_counter+1} has a different number of values.\n"
            error_msg += f"  Reference count: {len_ref}\n"
            error_msg += f"  Output count:    {len_out}"
            self.assertTrue(len_ref == len_out, msg=error_msg + self.info_msg)

            # 3. Compare each pair of numerical values from the line.
            for val_ref_str, val_out_str in zip(values_ref_str, values_out_str):
                try:
                    val_ref = float(val_ref_str)
                    val_out = float(val_out_str)

                    # Perform the hardcoded comparison.
                    is_close = math.isclose(val_ref, val_out, abs_tol=HARDCODED_TOLERANCE)

                    # Use the test framework's assertion to check the result and report failure.
                    full_msg = self.info_msg + "\n"
                    full_msg += f"Mismatch on data line {line_counter}: {val_ref} != {val_out}"
                    full_msg += f" (hardcoded tolerance = {HARDCODED_TOLERANCE})"
                    self.assertTrue(is_close, msg=full_msg)

                except ValueError:
                    # If conversion to float fails, fail the test immediately.
                    self.fail(f"Could not convert value to float on line {i+1}. "
                              f"Ref: '{val_ref_str}', Out: '{val_out_str}'" + self.info_msg)

            # Update counter to move to the next line.
            line_counter += 1

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
        if string != "":
            list_values.append(float(string))

    return list_values
