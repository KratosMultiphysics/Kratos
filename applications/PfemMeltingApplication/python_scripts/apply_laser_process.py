import KratosMultiphysics
import KratosMultiphysics.PfemMeltingApplication as PfemM

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ApplyLaserProcess(Model, settings["Parameters"])

class ApplyLaserProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)
        # Check the default values
        default_settings = KratosMultiphysics.Parameters( """
        {
            "model_part_name" : "CHOOSE_FLUID_MODELPART_NAME",
            "filename"        : "provide_the_name_of_the laser_file"
        }  """ )


        # Get the fluid model part from the Model container
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]

        #first laser
        with open("LaserSettings.json",'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        self.radius=project_parameters["laser_settings"]["laser_profile"]["radius"].GetDouble()
        self.shape=project_parameters["laser_settings"]["laser_profile"]["shape"].GetString()

        self.path=project_parameters["laser_settings"]["path"]

        list_of_coordinates = []
        list_of_power = []
        if(self.shape=="custom"):

            len = project_parameters["laser_settings"]["laser_profile"]["values"]
            i=0
            while(i < len.size()):
                coordinates = float(len[i]["distance"].GetDouble())
                #power = float(len[i]["power_per_unit_area"].GetDouble())
                power = float(len[i]["power_deviation_from_flat"].GetDouble())

                list_of_coordinates.append(coordinates)
                list_of_power.append(power)
                i = i + 1

            sum=0.0
            for i in range(len.size()-1):
                c_i_sq = list_of_coordinates[i] * list_of_coordinates[i]
                c_i1_sq = list_of_coordinates[i+1] * list_of_coordinates[i+1]
                c_i_cu = c_i_sq * list_of_coordinates[i]
                c_i1_cu = c_i1_sq * list_of_coordinates[i+1]
                slope = (list_of_power[i+1] - list_of_power[i]) / (list_of_coordinates[i+1] - list_of_coordinates[i])
                sum += 0.5 * list_of_power[i] * (c_i1_sq-c_i_sq) + slope * 0.3333333333333333 * (c_i1_cu - c_i_cu)


            total_heat=2.0 * 3.1416 * sum

        self.laser_path = project_parameters["laser_settings"]["path"]

        i=0

        self.new_table_x = KratosMultiphysics.PiecewiseLinearTable()
        self.new_table_y = KratosMultiphysics.PiecewiseLinearTable()
        self.new_table_z = KratosMultiphysics.PiecewiseLinearTable()
        self.new_table_Q = KratosMultiphysics.PiecewiseLinearTable()

        self.maximum_time= self.laser_path[self.laser_path.size()-1]["time"].GetDouble()

        while(i < self.laser_path.size()):
            time=self.laser_path[0]["time"].GetDouble()
            self.new_table_x.AddRow(self.laser_path[i]["time"].GetDouble(), self.laser_path[i]["x"].GetDouble())
            self.new_table_y.AddRow(self.laser_path[i]["time"].GetDouble(), self.laser_path[i]["y"].GetDouble())
            self.new_table_z.AddRow(self.laser_path[i]["time"].GetDouble(), self.laser_path[i]["z"].GetDouble())
            self.new_table_Q.AddRow(self.laser_path[i]["time"].GetDouble(), self.laser_path[i]["power"].GetDouble())
            i = i + 1

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
            Q=self.new_table_Q.GetValue(current_time)

        else:
            x=self.laser_path[self.laser_path.size()-1]["x"].GetDouble()
            y=self.laser_path[self.laser_path.size()-1]["y"].GetDouble()
            z=self.laser_path[self.laser_path.size()-1]["z"].GetDouble()
            Q=self.laser_path[self.laser_path.size()-1]["power"].GetDouble()

        self.ApplyLaserProcess.ApplyLaser(x, y, z, Q)



