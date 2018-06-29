from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import os

class TimeBasedAsciiFileWriterUtility(object):
    """This utility handles a file to which results are to be written.
    It internally handles restarts and the closing/opening of files.

    It should be used by e.g. processes that write data to files
    """
    def __init__(self, model_part, params, file_header):

        default_settings = KratosMultiphysics.Parameters('''{
            "output_file_name"  : "",
            "save_in_folder"    : true,
            "output_folder"     : "TimeBasedAsciiResults",
            "write_buffer_size" : -1
        }''')
        # write_buffer_size: -1 means we use the system default

        self.model_part = model_part

        params.ValidateAndAssignDefaults(default_settings)

        # file name and folder path specifications and check
        self.output_file_name = params["output_file_name"].GetString()
        self.save_in_folder = params["save_in_folder"].GetBool()
        self.output_folder = params["output_folder"].GetString()
        self.__ValidateAndAssignOutputFolderPath()

        # size of the buffer in bytes. Set to "0" for flushing always
        self.write_buffer_size = params["write_buffer_size"].GetInt()

        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.file = self.__InitializeOutputFile(file_header)
        else:
            restart_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            out_file = self.__AddToExistingOutputFile(file_header, restart_time)

            if out_file is not None:
                self.file = out_file
            else:
                warn_msg  = "No valid data file was found after restarting,\n"
                warn_msg += "writing to a new file"
                KratosMultiphysics.Logger.PrintWarning("TimeBasedAsciiFileWriterUtility", warn_msg)
                self.file = self.__InitializeOutputFile(file_header)

    def __OpenOutputFile(self):
        return open(self.output_file_name,"w", buffering=self.write_buffer_size)

    def __InitializeOutputFile(self, file_header):
        output_file = self.__OpenOutputFile()
        output_file.write(file_header)
        return output_file

    def __AddToExistingOutputFile(self, file_header, restart_time):
        if not os.path.isfile(self.output_file_name):
            return None

        try: # We try to open the file and transfer the info
            with open(self.output_file_name,'r') as out_file:
                lines_existing_file = out_file.readlines()

            # search for time, return false if it was not found
            # copy corresponding lines to new file and open it
            is_found = False

            # comparing the header
            old_header = ""
            while lines_existing_file[0].lstrip()[0] == '#':
                old_header += lines_existing_file.pop(0)

            if file_header != old_header:
                warn_msg  = "Headers in " + self.output_file_name + " after restarting do not match, \n"
                warn_msg += "appending results after restart from time " + str(restart_time) + " not possible"
                KratosMultiphysics.Logger.PrintWarning("TimeBasedAsciiFileWriterUtility", warn_msg)
                return None

            # comparing the data
            new_lines =""
            for line in lines_existing_file:
                new_lines+=line
                time_in_file = float(line.split()[0])
                if abs(time_in_file - restart_time) < 1e-20:
                    is_found = True
                    break

            if not(is_found):
                warn_msg  = "No line was found in " + self.output_file_name + " after restarting containing indicated restart time, \n"
                warn_msg += "appending results after restart from time " + str(restart_time) + " not possible"
                KratosMultiphysics.Logger.PrintWarning("TimeBasedAsciiFileWriterUtility", warn_msg)
                return None

            output_file = self.__OpenOutputFile() # this overwrites the old file
            output_file.write(file_header)
            output_file.write(new_lines)

            return output_file
        except:
            return None

    def __ValidateAndAssignOutputFolderPath(self):
        # check and correct file extension
        if not self.output_file_name.endswith(".dat"):
            self.output_file_name += ".dat"

        # check if relative path was erroneously specified in file name
        raw_path, raw_output_file_name = os.path.split(self.output_file_name)

        if self.output_folder != "":
            if raw_path != "":
                # assign the default value
                self.output_folder = "TimeBasedAsciiResults"

                warn_msg  = 'Relative path "'+ raw_path +'" contained wrongly in "output_file_name": "'+ self.output_file_name +'"\n'
                warn_msg += 'Use parameter "output_folder" to specify correctly\n'
                warn_msg += 'Using the default relative path "' + self.output_folder + '" instead'
                KratosMultiphysics.Logger.PrintWarning("TimeBasedAsciiFileWriteUtility", warn_msg)

            absolute_folder_path = os.path.join(os.getcwd(), self.output_folder)
            self.output_file_name = os.path.join(self.output_folder, raw_output_file_name)

        else:
            if raw_path != "":
                warn_msg  = 'Relative path "'+ raw_path +'" contained wrongly in "output_file_name": "'+ self.output_file_name +'"\n'
                warn_msg += 'Use the parameters "save_in_folder" and "output_folder" to specify correctly\n'
                warn_msg += 'Using the current directory instead'
                KratosMultiphysics.Logger.PrintWarning("TimeBasedAsciiFileWriteUtility", warn_msg)

            absolute_folder_path = os.getcwd()
            self.output_file_name = raw_output_file_name

        # make sure that the absolute path to the desired output folder exists
        if not os.path.isdir(absolute_folder_path):
            os.makedirs(absolute_folder_path)

