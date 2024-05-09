# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import math

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return MPMPointOutputProcess(model, settings["Parameters"])


class MPMPointOutputProcess(KratosMultiphysics.Process):
    """
    Auxiliary base class to output total flow forces
    over obstacles in fluid dynamics problems.
    A derived class needs to be implemented to be able to use
    this functionality, as calling the base class alone is not enough.
    """
    def __init__(self, model, params ):
        """
        Auxiliary class to output total flow forces over obstacles
        in fluid dynamics problems for a body fitted model part.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters('''{
            "help"                 : "This process writes results from a geometrical position (point) in the model to a file. It first searches the entity containing the requested output location and then interpolates the requested variable(s). The output can be requested for elements, conditions and nodes. For nodes no geometrical interpolation is performed, the exact coordinates have to be specified. This process works in MPI as well as with restarts. It can serve as a basis for other processes (e.g. MultiplePointsOutputProcess). Furthermore it can be used for testing in MPI where the node numbers can change",
            "model_part_name"      : "",
            "interval"             : [0.0, 1e30],
            "position"             : [0, 0, 0],
            "print_format"         : ".8f",
            "write_tracking_output_file" : true,
            "output_pressure"      : false,
            "output_file_settings": {}
            }''')


        #default_settings = KratosMultiphysics.Parameters("""
        #    {
        #        "model_part_name"           : "",
        #        "interval"                  : [0.0, 1e30],
        #        "print_drag_to_screen"      : false,
        #        "print_format"              : ".8f",
        #        "write_drag_output_file"    : true,
        #        "output_file_settings": {}
        #    }
        #    """)

        self.params = params
        # Detect "End" as a tag and replace it by a large number
        if(self.params.Has("interval")):
            if(self.params["interval"][1].IsString()):
                if(self.params["interval"][1].GetString() == "End"):
                    self.params["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+self.params["interval"].PrettyPrintJsonString())

        self.params.ValidateAndAssignDefaults(default_settings)

        # variables for search of given position in model part
        #self.search_done = False

        self.format = self.params["print_format"].GetString()

        # getting the ModelPart from the Model
        model_part_name = self.params["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        else:
            self.model_part = model[self.params["model_part_name"].GetString()]

    def ExecuteBeforeSolutionLoop(self):

        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = self.params["interval"][0].GetDouble()
        self.interval[1] = self.params["interval"][1].GetDouble()
        #self.print_tracking_to_screen = self.params["print_tracking_to_screen"].GetBool()
        self.write_drag_output_file = self.params["write_tracking_output_file"].GetBool()
        self.output_press = self.params["output_pressure"].GetBool()

        #Read position point
        point_position = self.params["position"].GetVector()
        if point_position.Size() != 3:
            raise Exception('The position has to be provided with 3 coordinates!')
        self.initial_point = KratosMultiphysics.Point(point_position[0],
                                         point_position[1],
                                         point_position[2])

        if (self.write_drag_output_file):
            #if (self.model_part.GetCommunicator().MyPID() == 0):
                #output_file_name = self.params["model_part_name"].GetString() + "_tracking.dat"
                file_handler_params = KratosMultiphysics.Parameters(self.params["output_file_settings"])
                #file_handler_params.AddEmptyValue("file_name")
                #file_handler_params["file_name"].SetString(output_file_name)
                file_header = "" #self._GetFileHeader()
                self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part,
                    file_handler_params, file_header).file


        self.particle_id = 0
        minor_distance = 0
        for  mp in self.model_part.Elements:
            mp_coord = mp.CalculateOnIntegrationPoints(KratosParticle.MP_COORD,self.model_part.ProcessInfo)[0]
            distance=ComputeEucledianDistance(mp_coord, self.initial_point)
            if (minor_distance ==0 or minor_distance > distance):
                minor_distance = distance
                self.particle_id = mp.Id
                mp_coord_minor_distance = mp_coord
        if (minor_distance !=0) :
            print("----------------------------------------------")
            print("----Tracking of a material particle in time---")
            print("----------------------------------------------")
            print("Coordinates introduced: " + str(self.initial_point))
            print("Particle found!")
            print("Particle_id: " + str(self.particle_id))
            print("Coordinates of the particle: " + str(mp_coord_minor_distance))
            print("Please see the _tracking.dat file to see the data.")
            print("-----------------------------------------------")
        if (minor_distance==0):
            raise Exception("No particle found to do the tracking of time.")
        self.output_file.write("#Tracking of a material particle in time" + "\n")
        self.output_file.write("#Particle id: " + str(self.particle_id) + " with coordinates: " + str(mp_coord_minor_distance)+ " is the closest to the selected point: " + str(self.initial_point) + "\n")
        self.output_file.write("#The tracking time (s) for the mp is from " + str(self.interval[0]) + " to " + str(self.interval[1]) + " seconds." + "\n")
        self.output_file.write("#TIME" + " " + "COORD_X" + " " + "COORD_Y"+ " " + "COORD_Z" + " "  + "VELOC_X"+ " " + "VELOC_Y" + "" + " " + "VELOC_Z"+ " "+ "DISPL_X"+ " " + "DISPL_Y" + " " + "DISPL_Z")
        if (self.output_press==True):
            self.output_file.write(" " + "PRESSURE"+ "\n")
        else:
            self.output_file.write("\n")


    def ExecuteFinalizeSolutionStep(self):
        current_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        #self.initial_point = KratosMultiphysics.Point(4.8, 6, 0)
        variable_disp = KratosMultiphysics.KratosGlobals.GetVariable( "MP_DISPLACEMENT" )
        variable_vel = KratosMultiphysics.KratosGlobals.GetVariable( "MP_VELOCITY" )
        if (self.output_press==True):
            variable_press = KratosMultiphysics.KratosGlobals.GetVariable( "MP_PRESSURE" )


        mp=self.model_part.Elements[self.particle_id]
        mp_coord = mp.CalculateOnIntegrationPoints(KratosParticle.MP_COORD,self.model_part.ProcessInfo)[0]
        displacement = mp.CalculateOnIntegrationPoints(variable_disp,self.model_part.ProcessInfo)[0]
        velocity = mp.CalculateOnIntegrationPoints(variable_vel,self.model_part.ProcessInfo)[0]
        if (self.output_press==True):
            pressure = mp.CalculateOnIntegrationPoints(variable_press,self.model_part.ProcessInfo)[0]

        # Write the tracking of points values
        #if (self.model_part.GetCommunicator().MyPID() == 0):
        if (current_time >= self.interval[0] and current_time <=self.interval[1]):
            self.output_file.write(str(current_time)+" "+format(mp_coord[0],self.format)+" "+format(mp_coord[1],self.format)+" "+format(mp_coord[2],self.format)+" "+format(velocity[0],self.format)+" "+format(velocity[1],self.format)+" "+format(velocity[2],self.format)+" "+format(displacement[0],self.format)+" "+format(displacement[1],self.format)+" "+format(displacement[2],self.format)+" ")
            if (self.output_press==True):
                self.output_file.write(format(pressure,self.format)+"\n")
            else:
                self.output_file.write("\n")

    def ExecuteFinalize(self):
        #if (self.write_tracking_output_file):
            #if (self.model_part.GetCommunicator().MyPID() == 0):
        self.output_file.close()

def ComputeEucledianDistance(particle, point):
    distance = math.sqrt(math.pow(particle[0] - point[0],2) + math.pow(particle[1] - point[1],2) + math.pow(particle[2] - point[2],2))
    return distance
