import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import math
import pre_calculated_fluid_analysis
BaseAnalysis = pre_calculated_fluid_analysis.PreCalculatedFluidAnalysis

class RotatingAleAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)
        self.SetRotator()
        self.time_full_flow = 0.05
        self.inlet_group_number = 4
        self.inlet_velocity = - 9.06

    def SetRotator(self):
        self.rotator = SDEM.MeshRotationUtility(self.project_parameters)

    def SetBetaParameters(self):
        BaseAnalysis.SetBetaParameters(self)
        self.project_parameters.AddEmptyValue("ALE_option").SetBool(True)
        Add = self.project_parameters.AddEmptyValue
        Add("steps_per_average_step").SetInt(1)
        Add("initial_averaging_time").SetDouble(0.0)
        Add("rotated_stationary_flow_option").SetBool(False)
        Add("averaging_has_already_been_done").SetBool(False)
        Add("stationary_start_time").SetDouble(0.0)

    def UpdateALEMeshMovement(self, time):
        if self.project_parameters["custom_fluid"]["ALE_option"].GetBool():
            self.rotator.RotateMesh(self.fluid_model_part, time)
            self.projection_module.UpdateDatabase(self.h_min)

    def AssessStationarity(self):
        # BaseAnalysis.AssessStationarity(self)
        if self.time > self.project_parameters["stationary_start_time"].GetDouble():
           self.stationarity = True

        if self.stationarity:
            self.rotator.SetStationaryField(self.fluid_model_part, self.time)

    def SetFluidLoader(self):
        if self.project_parameters["rotated_stationary_flow_option"].GetBool():
            import KratosMultiphysics.SwimmingDEMApplication.hdf5_io_tools_PID as hdf5_io_tools_PID
            import KratosMultiphysics.SwimmingDEMApplication.average_field as average_field

            rotation_axis_initial_point = self.project_parameters['frame_of_reference']["frame_rotation_axis_initial_point"].GetVector()
            rotation_axis_final_point = self.project_parameters['frame_of_reference']["frame_rotation_axis_final_point"].GetVector()
            angular_velocity_module = self.project_parameters['frame_of_reference']["angular_velocity_magnitude"].GetDouble()
            dataset_name = 'stationary_field'
            original_file_name = self.project_parameters["prerun_fluid_file_name"].GetString()
            initial_averaging_time = self.project_parameters["initial_averaging_time"].GetDouble()
            steps_per_average_step = self.project_parameters["steps_per_average_step"].GetInt()

            averager = average_field.Averager(rotation_axis_initial_point = rotation_axis_initial_point,
                                              rotation_axis_final_point = rotation_axis_final_point,
                                              angular_velocity_module = angular_velocity_module,
                                              dataset_name = dataset_name,
                                              original_file_name = original_file_name,
                                              original_file_path = self.main_path,
                                              initial_time = initial_averaging_time,
                                              steps_per_average_step = steps_per_average_step)

            self.fluid_loader = hdf5_io_tools_PID.FluidHDF5LoaderPID(self.project_parameters,
                                                                     self.all_model_parts.Get('FluidPart'),
                                                                     self.all_model_parts.Get('SpheresPart'),
                                                                     self.main_path,
                                                                     averager)
        else:
            BaseAnalysis.SetFluidLoader(self)

    def FluidSolve(self, time='None', solve_system=True):
        self.SetInletVelocity(time)
        rotated_stationary_flow_option = self.project_parameters["rotated_stationary_flow_option"].GetBool()
        averaging_has_already_been_done = self.project_parameters["averaging_has_already_been_done"].GetBool()

        if rotated_stationary_flow_option and not averaging_has_already_been_done:
            self.fluid_loader.averager.PerformAverage(reference_time = time)
            self.project_parameters["averaging_has_already_been_done"].SetBool(True)

        BaseAnalysis.FluidSolve(self, time, solve_system)

        if rotated_stationary_flow_option and not solve_system:
            self.rotator.RotateFluidVelocities(time)

    def SetInletVelocity(self, time):
        if time <= self.time_full_flow:
            alpha = math.sin(0.5 * math.pi * time / self.time_full_flow)
            self.fluid_model_part.GetProperties()[self.inlet_group_number][Kratos.IMPOSED_VELOCITY_Z_VALUE] = alpha * self.inlet_velocity

            for node in self.fluid_model_part.GetMesh(self.inlet_group_number).Nodes:
                node.SetSolutionStepValue(Kratos.VELOCITY_Z, alpha * self.inlet_velocity)
        else:
            pass
