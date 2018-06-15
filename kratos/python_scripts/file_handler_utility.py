from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import os

class FileHandlerUtility(object):
    """This utility handles a file to which results are to be written.
    It internally handles restarts and the closing/opening of files.

    It should be used by e.g. processes that write data to files
    """
    def __init__(self, model_part, params, file_header):

        default_settings = KratosMultiphysics.Parameters('''{
            "output_file_name"     : "",
            "flush_frequency"      : ""
        }''')

        self.model_part = model_part

        params.ValidateAndAssignDefaults(default_settings)

        output_file_name = params["output_file_name"].GetString()

        flush_frequency = params["flush_frequency"].GetString()
        if flush_frequency != "": # user specified a flush frequency
            self.flush_output = True
            self.flush_frequency = float(flush_frequency)
            self.next_flush = self.flush_frequency
        else:
            self.flush_output = False

        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.file = InitializeOutputFile(output_file_name, file_header)
        else:
            restart_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            existing_file_is_valid, out_file = AddToExistingOutputFile(output_file_name, file_header, restart_time)

            if existing_file_is_valid:
                self.file = out_file
            else:
                warn_msg  = "No valid data file was found after restarting,\n"
                warn_msg += "writing to a new file"
                KratosMultiphysics.Logger.PrintWarning("FileHandlerUtility", warn_msg)
                self.file = InitializeOutputFile(output_file_name, file_header)

    def write(self, string_to_write):
        self.file.write(string_to_write)

        if self.flush_output:
            time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            if time > self.next_flush:
                self.flush()
            if self.flush_frequency > 0.0: # Note: if == 0, we'll always flush
                while self.next_flush <= time:
                    self.next_flush += self.flush_frequency

    def flush(self):
        self.file.flush()

    def close(self):
        self.__CloseOutputFiles()

    def __del__(self):
         # in case "close" is not called
         # this can happen if a simulation is forcefully stopped on a cluster
        self.__CloseOutputFiles()

    def __CloseOutputFiles(self):
        '''Close output file.'''
        if not self.file.closed:
            self.flush()
            self.file.close()


def InitializeOutputFile(output_file_name, file_header):
    output_file = open(output_file_name,"w")
    output_file.write(file_header)
    return output_file

def AddToExistingOutputFile(output_file_name, file_header, restart_time):
    if not os.path.isfile(output_file_name):
        return False, None

    try: # We try to open the file and transfer the info
        with open(output_file_name,'r') as out_file:
            lines_existing_file = out_file.readlines()


        # search for time, return false if it was not found
        # copy corresponding lines to new file and open it
        is_found = False

        # comparing the header
        old_header = ""
        while lines_existing_file[0].lstrip()[0] == '#':
            old_header += lines_existing_file.pop(0)

        if file_header != old_header:
            warn_msg  = "Headers in " + output_file_name + " after restarting do not match, \n"
            warn_msg += "appending results after restart from time " + str(restart_time) + " not possible"
            KratosMultiphysics.Logger.PrintWarning("FileHandlerUtility", warn_msg)
            return is_found, None

        # comparing the data
        new_lines =""
        for line in lines_existing_file:
            new_lines+=line
            time_in_file = float(line.split()[0])
            if abs(time_in_file - restart_time) < 1e-20:
                is_found = True
                break

        if not(is_found):
            warn_msg  = "No line was found in " + output_file_name + " after restarting containing indicated restart time, \n"
            warn_msg += "appending results after restart from time " + str(restart_time) + " not possible"
            KratosMultiphysics.Logger.PrintWarning("FileHandlerUtility", warn_msg)
            return is_found, None

        output_file = open(output_file_name,"w") # this overwrites the old file
        output_file.write(file_header)
        output_file.write(new_lines)

        return is_found, output_file
    except:
        return False, None

