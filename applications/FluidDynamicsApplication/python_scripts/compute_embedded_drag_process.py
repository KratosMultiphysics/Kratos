# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# other imports
from time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeEmbeddedDragProcess(model, settings["Parameters"])


class ComputeEmbeddedDragProcess(KratosMultiphysics.Process):
    def __init__(self, model, params ):
        """ 
        Auxiliary class to output total flow forces over obstacles 
        in fluid dynamics problems for an embedded model part.
        """
        super(ComputeEmbeddedDragProcess,self).__init__()

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"           : "",
                "interval"                  : [0.0, 1e30],
                "write_drag_output_file"    : true,
                "print_drag_to_screen"      : false,
                "write_buffer_size"         : -1,
                "print_format"              : ""
            }
            """)

        self.params = params
        # Detect "End" as a tag and replace it by a large number
        if(self.params.Has("interval")):
            if(self.params["interval"][1].IsString()):
                if(self.params["interval"][1].GetString() == "End"):
                    self.params["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+self.params["interval"].PrettyPrintJsonString())

        self.params.ValidateAndAssignDefaults(default_settings)

        self.format = self.params["print_format"].GetString()

        # getting the ModelPart from the Model
        model_part_name = self.params["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        else:
            self.model_part = model[self.params["model_part_name"].GetString()]

    def ExecuteInitialize(self):
    
        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = self.params["interval"][0].GetDouble()
        self.interval[1] = self.params["interval"][1].GetDouble()
        self.print_drag_to_screen = self.params["print_drag_to_screen"].GetBool()
        self.write_drag_output_file = self.params["write_drag_output_file"].GetBool()

        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_drag_output_file):
                
                output_file_name = self.params["model_part_name"].GetString() + "_drag.dat"

                file_handler_params = KratosMultiphysics.Parameters('''{ "output_file_name" : "" }''')
                
                file_handler_params["output_file_name"].SetString(output_file_name)
                file_handler_params.AddValue("write_buffer_size", self.params["write_buffer_size"])

                file_header = GetFileHeader(self.params["model_part_name"].GetString())
                self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part, 
                    file_handler_params, file_header).file

    def ExecuteFinalizeSolutionStep(self):
    
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if((current_time >= self.interval[0]) and  (current_time < self.interval[1])):
            # Compute the drag force
            drag_force = KratosCFD.DragUtilities().CalculateEmbeddedDrag(self.model_part)

            # Write the drag force values
            if (self.model_part.GetCommunicator().MyPID() == 0):
                if (self.print_drag_to_screen):
                    KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedDragProcess", "DRAG RESULTS:")
                    KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedDragProcess","Current time: " + str(current_time) + " x-drag: " + format(drag_force[0],self.format) + " y-drag: " + format(drag_force[1],self.format) + " z-drag: " + format(drag_force[2],self.format))

                if (self.write_drag_output_file):
                    self.output_file.write(str(current_time)+" "+format(drag_force[0],self.format)+" "+format(drag_force[1],self.format)+" "+format(drag_force[2],self.format)+"\n")

    def ExecuteFinalize(self):
        if (self.model_part.GetCommunicator().MyPID() == 0):
            self.output_file.close()

def GetFileHeader(model_part_name):

    header  = '# Embedded drag for model part ' + model_part_name + '\n'
    header += '# Time Fx Fy Fz \n'

    return header