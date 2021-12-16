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
          
        mat = project_parameters["laser_settings"]["path"]
        #table1=mat["Tables"]["Table1"]
        print(mat.size())
        i=0
        #taux1= table1["Table1"]
        self.new_table_x = KratosMultiphysics.PiecewiseLinearTable()
        self.new_table_y = KratosMultiphysics.PiecewiseLinearTable()
        self.new_table_z = KratosMultiphysics.PiecewiseLinearTable()
        while(i < mat.size()):
            time=mat[0]["time"].GetDouble()
            self.new_table_x.AddRow(mat[i]["time"].GetDouble(), mat[i]["x"].GetDouble())
            self.new_table_y.AddRow(mat[i]["time"].GetDouble(), mat[i]["y"].GetDouble())
            self.new_table_z.AddRow(mat[i]["time"].GetDouble(), mat[i]["z"].GetDouble())
            i = i + 1


       
        # Set the Boussinesq force process
        self.ApplyLaserProcess = PfemM.ApplyLaserProcess(self.fluid_model_part, project_parameters)

    def ExecuteInitialize(self):
        self.ApplyLaserProcess.ExecuteInitialize()
        

        

    def ExecuteInitializeSolutionStep(self):
        self.ApplyLaserProcess.ExecuteInitializeSolutionStep()
        current_time = self.fluid_model_part.ProcessInfo[KratosMultiphysics.TIME]

        x=self.new_table_x.GetValue(current_time)
        y=self.new_table_y.GetValue(current_time)
        z=self.new_table_z.GetValue(current_time)

        self.ApplyLaserProcess.laserapplication(self.radius,self.power, x, y, z)
         
