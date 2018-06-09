import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import python_process

# other imports
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeDragProcess(Model, settings["Parameters"])


class ComputeDragProcess(python_process.PythonProcess):
    def __init__(self, Model, settings ):

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"                   : 0,
                "model_part_name"           : "please_specify_model_part_name",
                "interval"                  : [0.0, 1e30],
                "write_drag_output_file"    : true,
                "print_drag_to_screen"      : false
            }
            """)

        # Detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString()):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = settings["interval"][0].GetDouble()
        self.interval[1] = settings["interval"][1].GetDouble()
        self.print_drag_to_screen = settings["print_drag_to_screen"].GetBool()
        self.write_drag_output_file = settings["write_drag_output_file"].GetBool()

        if (self.write_drag_output_file) and (self.model_part.GetCommunicator().MyPID() == 0):
          
            output_file_name = settings["model_part_name"].GetString() + ".drag"

            if self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                # if the simulation is restarted we search for a  file
                # existing from the previous run to append the data
                restart_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
                existing_file_is_valid, out_file = AddToExistingOutputFile(output_file_name,
                                                                        restart_time)

                if existing_file_is_valid:
                    self.output_file = out_file
                # if no valid file can be found we create a new one
                # and issue a warning
                else:
                    warn_msg  = "No data file was found after restarting,\n"
                    warn_msg += "writing to a new file"
                    KratosMultiphysics.Logger.PrintWarning("ComputeDragProcess", warn_msg)

                    self.output_file = InitializeOutputFile(output_file_name, settings["model_part_name"].GetString())

            else: # no restart, regular simulation
                self.output_file = InitializeOutputFile(output_file_name, settings["model_part_name"].GetString())

    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if((current_time >= self.interval[0]) and  (current_time < self.interval[1])):
            # Compute the drag force
            drag_force = KratosFluid.DragUtilities().CalculateBodyFittedDrag(self.model_part)

            # Write the drag force values
            if (self.model_part.GetCommunicator().MyPID() == 0):
                # Print drag values to screen
                if (self.print_drag_to_screen):
                    print("DRAG RESULTS:")
                    print("Current time: " + str(current_time) + " x-drag: " + str(drag_force[0]) + " y-drag: " + str(drag_force[1]) + " z-drag: " + str(drag_force[2]))
                
                # Print drag values to file
                if (self.write_drag_output_file):
                    self.output_file.write(str(current_time)+"   "+str(drag_force[0])+"   "+str(drag_force[1])+"   "+str(drag_force[2])+"\n")

    def ExecuteFinalize(self):
        self.__CloseOutputFile()

    def __del__(self):
         # in case "ExecuteFinalize" is not called
         # this can happen if a simulation is forcefully stopped on a cluster
        self.__CloseOutputFile()

    def __CloseOutputFile(self):
        '''Close output file.'''
        self.output_file.close()

def InitializeOutputFile(output_file_name, model_part_name):
    output_file = open(output_file_name,"w")
    out  = '# ' + model_part_name + ' drag \n'
    out += '# Time Fx Fy Fz \n'
    
    output_file.write(out)

    return output_file

def AddToExistingOutputFile(output_file_name, restart_time):
    if not os.path.isfile(output_file_name):
        return False, None

    try: # We try to open the file and transfer the info
        with open(output_file_name,'r') as out_file:
            lines_existing_file = out_file.readlines()

        output_file = open(output_file_name,"w") # this overwrites the old file

        # search for time, return false if it was not found
        # copy corresponding lines to new file and open it
        is_found = False

        for line in lines_existing_file:
            output_file.write(line)
            if line.startswith(str(restart_time)):
                is_found = True
                break

        if not(is_found):
            warn_msg  = "No line was found in " + output_file_name + " after restarting containing indicated restart time, \n"
            warn_msg += "appending results after restart from time " + str(restart_time) + " not possible"
            KratosMultiphysics.Logger.PrintWarning("ComputeDragProcess", warn_msg)

        return True, output_file
    except:
        return False, None