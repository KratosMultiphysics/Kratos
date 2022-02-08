import KratosMultiphysics
import KratosMultiphysics.PfemMeltingApplication as PfemM

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ApplyLaserProcess(Model, settings["Parameters"])

class ApplyLaserProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)
        print("here")
        # Check the default values
        default_settings = KratosMultiphysics.Parameters( """
        {
            "model_part_name" : "CHOOSE_FLUID_MODELPART_NAME",
            "filename"        : "provide_the_name_of_the laser_file"
        }  """ )


        # Get the fluid model part from the Model container
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]


        #print(settings["filename"])

        with open("LaserSettings.json",'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        #print(project_parameters)
        #print(project_parameters["laser_settings"]["laser_profile"])

        self.radius=project_parameters["laser_settings"]["laser_profile"]["radius"].GetDouble()
        self.power=project_parameters["laser_settings"]["laser_profile"]["power"].GetDouble()
        self.shape=project_parameters["laser_settings"]["laser_profile"]["shape"].GetString()

        self.path=project_parameters["laser_settings"]["path"]
        #print(self.path)

        laser_path = project_parameters["laser_settings"]["path"]
        #table1=mat["Tables"]["Table1"]
        print(laser_path.size())
        i=0
        #taux1= table1["Table1"]
        self.new_table_x = KratosMultiphysics.PiecewiseLinearTable()
        self.new_table_y = KratosMultiphysics.PiecewiseLinearTable()
        self.new_table_z = KratosMultiphysics.PiecewiseLinearTable()

        self.maximum_time= laser_path[laser_path.size()-1]["time"].GetDouble()
        
        while(i < laser_path.size()):
            time=laser_path[0]["time"].GetDouble()
            self.new_table_x.AddRow(laser_path[i]["time"].GetDouble(), laser_path[i]["x"].GetDouble())
            self.new_table_y.AddRow(laser_path[i]["time"].GetDouble(), laser_path[i]["y"].GetDouble())
            self.new_table_z.AddRow(laser_path[i]["time"].GetDouble(), laser_path[i]["z"].GetDouble())
            i = i + 1

        # Set the Boussinesq force process
        self.ApplyLaserProcess = PfemM.ApplyLaserProcess(self.fluid_model_part, project_parameters["laser_settings"])

    def ExecuteInitialize(self):
        self.ApplyLaserProcess.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.ApplyLaserProcess.ExecuteInitializeSolutionStep()
        current_time = self.fluid_model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(current_time < self.maximum_time):
            x=self.new_table_x.GetValue(current_time)
            y=self.new_table_y.GetValue(current_time)
            z=self.new_table_z.GetValue(current_time)
        else:
            x=laser_path[laser_path.size()-1]["x"].GetDouble()
            y=laser_path[laser_path.size()-1]["y"].GetDouble()
            z=laser_path[laser_path.size()-1]["z"].GetDouble()
	
        
        self.ApplyLaserProcess.ApplyLaser(x, y, z)

