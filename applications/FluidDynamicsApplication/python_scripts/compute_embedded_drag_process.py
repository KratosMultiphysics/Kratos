import KratosMultiphysics
import python_process

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeEmbeddedDragProcess(Model, settings["Parameters"])


class ComputeEmbeddedDragProcess(python_process.PythonProcess):
    def __init__(self, Model, settings ):

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"                   : 0,
                "model_part_name"           : "",
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

        self.Model = Model
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = settings["interval"][0].GetDouble()
        self.interval[1] = settings["interval"][1].GetDouble()
        self.print_drag_to_screen = settings["print_drag_to_screen"].GetBool()
        self.write_drag_output_file = settings["write_drag_output_file"].GetBool()

        if (self.write_drag_output_file) and (self.fluid_model_part.GetCommunicator().MyPID()==0):
            # Set drag output file name
            self.drag_filename = settings["model_part_name"].GetString() + ".drag"

            # File creation to store the drag evolution
            with open(self.drag_filename, 'w') as file:
                file.write("Embedded objects total drag")
                file.write("\n")
                file.write("Time   Fx   Fy   Fz \n")
                file.close()


    def ExecuteFinalizeSolutionStep(self):

        current_time = self.fluid_model_part.ProcessInfo[KratosMultiphysics.TIME]

        if((current_time >= self.interval[0]) and (current_time < self.interval[1])):

            drag_x = KratosMultiphysics.VariableUtils().SumHistoricalNodeScalarVariable(KratosMultiphysics.REACTION_X, self.fluid_model_part, 0)
            drag_y = KratosMultiphysics.VariableUtils().SumHistoricalNodeScalarVariable(KratosMultiphysics.REACTION_Y, self.fluid_model_part, 0)
            drag_z = KratosMultiphysics.VariableUtils().SumHistoricalNodeScalarVariable(KratosMultiphysics.REACTION_Z, self.fluid_model_part, 0)

            if (self.print_drag_to_screen) and (self.fluid_model_part.GetCommunicator().MyPID()==0):
                print("DRAG RESULTS:")
                print("Current time: " + str(current_time) + " x-drag: " + str(drag_x) + " y-drag: " + str(drag_y) + " z-drag: " + str(drag_z))

            if (self.write_drag_output_file) and (self.fluid_model_part.GetCommunicator().MyPID()==0):
                with open(self.drag_filename, 'a') as file:
                    file.write(str(current_time)+"   "+str(drag_x)+"   "+str(drag_y)+"   "+str(drag_z)+"\n")
                    file.close()
