# Importing the Kratos Library
import KratosMultiphysics

# Other imports
from pathlib import Path

class TimeBasedAsciiFileWriterUtility:
    """This utility handles a file to which results are to be written.
    It internally handles restarts and the closing/opening of files.

    It should be used by e.g. processes that write data to files
    """
    def __init__(self, model_part, params, file_header):

        default_settings = KratosMultiphysics.Parameters('''{
            "file_name"  : "",
            "output_path": "",
            "write_buffer_size" : -1
        }''')
        # write_buffer_size: -1 means we use the system default
        # write_buffer_size:  0 means no buffering is done. IMPORTANT : Only for binary output.
        # write_buffer_size > 0 means value specified is the size of buffer

        self.model_part = model_part
        has_initial_write_buffer_size = params.Has("write_buffer_size")

        if params.Has("folder_name"):
            params.AddValue("output_path",params["folder_name"])
            params.RemoveValue("folder_name")
            KratosMultiphysics.Logger.PrintWarning('TimeBasedAsciiFileWriterUtility', '"folder_name" key is deprecated. Use "output_path" instead.')
        params.ValidateAndAssignDefaults(default_settings)

        # file name and folder path specifications and check
        self.file_name = Path(params["file_name"].GetString())
        self.output_path = Path(params["output_path"].GetString())
        self.__ValidateAndAssignOutputFolderPath()

        # size of the buffer in bytes. Set to "0" for flushing always
        is_mpi_execution = (model_part.GetCommunicator().TotalProcesses() > 1)
        if not is_mpi_execution and not has_initial_write_buffer_size:
            info_msg  = "File output buffer size set to 1 \n"
            info_msg += "for TimeBasedAsciiFileWriterUtility output file "+ str(self.file_name)
            KratosMultiphysics.Logger.PrintInfo("TimeBasedAsciiFileWriterUtility", info_msg)
            self.write_buffer_size = 1
        else:
            self.write_buffer_size = params["write_buffer_size"].GetInt()
            if (self.write_buffer_size == 0):
                err_msg  = "Buffer size of 0 not possible for ASCII output. \n"
                err_msg  += "\t\tPlease choose a number greater than 1 or set -1 for default size. \n"
                raise Exception(err_msg)

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
        return open(str(self.file_name),"w", buffering=self.write_buffer_size)

    def __InitializeOutputFile(self, file_header):
        output_file = self.__OpenOutputFile()
        output_file.write(file_header)
        return output_file

    def __AddToExistingOutputFile(self, file_header, restart_time):
        if not self.file_name.is_file():
            return None

        try: # We try to open the file and transfer the info
            with open(str(self.file_name),'r') as out_file:
                lines_existing_file = out_file.readlines()

            # search for time, return false if it was not found
            # copy corresponding lines to new file and open it
            is_found = False

            # comparing the header
            old_header = ""
            while lines_existing_file[0].lstrip()[0] == '#':
                old_header += lines_existing_file.pop(0)

            if file_header != old_header:
                warn_msg  = "Headers in " + self.file_name + " after restarting do not match, \n"
                warn_msg += "appending results after restart from time " + str(restart_time) + " not possible"
                KratosMultiphysics.Logger.PrintWarning("TimeBasedAsciiFileWriterUtility", warn_msg)
                return None

            # comparing the data
            new_lines =""
            for line in lines_existing_file:
                new_lines+=line
                time_in_file = float(line.split()[0])
                if abs(time_in_file - restart_time) < 1e-12:
                    is_found = True
                    break

            if is_found:
                info_msg = 'Restart for file "{}" successful'.format(self.file_name)
                KratosMultiphysics.Logger.PrintInfo("TimeBasedAsciiFileWriterUtility", info_msg)
            else:
                warn_msg  = "No line was found in " + str(self.file_name) + " after restarting containing indicated restart time, \n"
                warn_msg += "appending results after restart from time " + str(restart_time) + " not possible.\n"
                warn_msg += "To avoid loss of data continuing writing from the end of file\n"
                KratosMultiphysics.Logger.PrintWarning("TimeBasedAsciiFileWriterUtility", warn_msg)

                new_lines+= "# indicated restart time " + str(restart_time) + " not found (up to the given tolerance)\n"
                new_lines+= "# continuing writing after restart from the previous end of file\n"

            output_file = self.__OpenOutputFile() # this overwrites the old file
            output_file.write(file_header)
            output_file.write(new_lines)

            return output_file
        except:
            return None

    def __ValidateAndAssignOutputFolderPath(self):
        # a file name must be specified
        if self.file_name == "":
            raise Exception('No "file_name" was specified!')

        # check and correct file extension
        if self.file_name.suffix != ".dat":
            self.file_name = self.file_name.with_suffix(".dat")

        if self.file_name.parent != Path(): # check if folder was specified in "file_name"
            warn_msg  = 'Path contained wrongly in "file_name" "{}", this will be ignored in the future\n'.format(self.file_name)
            warn_msg += 'Use parameter "output_path" to specify it correctly.'
            KratosMultiphysics.Logger.PrintWarning("TimeBasedAsciiFileWriterUtility", warn_msg)

        self.file_name = self.output_path / self.file_name

        # make sure that the path to the desired output folder exists
        KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories(str(self.output_path))
