import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ShallowWaterGiDOutputProcess(model, settings["Parameters"])

from KratosMultiphysics.gid_output_process import GiDOutputProcess

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class ShallowWaterGiDOutputProcess(GiDOutputProcess):
    def __init__(self, model, settings):
        first_level_defaults = KM.Parameters("""
        {
            "output_name"               : "file",
            "model_part_name"           : "unknown",
            "result_file_configuration" : {},
            "point_data_configuration"  : []
        }
        """)
        settings.ValidateAndAssignDefaults(first_level_defaults)
        self.param = settings

        self.base_file_name = self.param["output_name"].GetString()
        self.model_part = model[self.param["model_part_name"].GetString()]
        self.body_io = None
        self.volume_list_files = []

        # Initializing the topographic post process tool
        self.topographic_tool = SW.PostProcessUtilities(self.model_part)
        self.topographic_tool.DefineAuxiliaryProperties()

        # This output does not support cuts. Initializing empty instances
        self.cut_model_part = None
        self.cut_io = None
        self.output_surface_index = 0
        self.cut_list_files = []

        # I am not interested on point output, but initialize that stuff
        point_data_configuration = self.param["point_data_configuration"]
        if point_data_configuration.size() > 0:
            import point_output_process
            self.point_output_process = point_output_process.PointOutputProcess(self.model_part,point_data_configuration)
        else:
            self.point_output_process = None

        self.step_count = 0
        self.printed_step_count = 0
        self.next_output = 0.0

    def ExecuteInitialize(self):
        super(ShallowWaterGiDOutputProcess, self).ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        super(ShallowWaterGiDOutputProcess, self).ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        super(ShallowWaterGiDOutputProcess, self).ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        super(ShallowWaterGiDOutputProcess, self).ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        pass

    def IsOutputStep(self):
        return super(ShallowWaterGiDOutputProcess, self).IsOutputStep()

    def PrintOutput(self):
        super(ShallowWaterGiDOutputProcess, self).PrintOutput()

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        super(ShallowWaterGiDOutputProcess, self).ExecuteFinalize()

    def Check(self):
        return 1

    def _GiDOutputProcess__write_mesh(self, label):
        if self.body_io is not None:
            self.body_io.InitializeMesh(label)
            if self.body_output:
                # Write the standard mesh
                self.body_io.WriteMesh(self.model_part.GetMesh())

                # Write the topographic mesh
                self.__prepare_topographic_mesh()
                self.body_io.WriteMesh(self.model_part.GetMesh())
                self.__reset_standard_mesh()

            if self.node_output:
                self.body_io.WriteNodeMesh(self.model_part.GetMesh())
            self.body_io.FinalizeMesh()

        if self.cut_io is not None:
            self.cut_io.InitializeMesh(label)
            self.cut_io.WriteMesh(self.cut_model_part.GetMesh())
            self.cut_io.FinalizeMesh()

    def _GiDOutputProcess__write_nodal_results(self, label):
        if self.body_io is not None:
            # Write the standard results
            self.__prepare_standard_mesh()
            for variable in self.nodal_variables:
                self.body_io.WriteNodalResults(variable, self.model_part.GetCommunicator().LocalMesh().Nodes, label, 0)

            # Write the topographic results
            self.__prepare_topographic_mesh()
            self.body_io.WriteNodalResults(SW.BATHYMETRY, self.model_part.GetCommunicator().LocalMesh().Nodes, label, 0)
            self.__reset_standard_mesh()

        if self.cut_io is not None:
            self.cut_manager.UpdateCutData(self.cut_model_part, self.model_part)
            for variable in self.nodal_variables:
                self.cut_io.WriteNodalResults(variable, self.cut_model_part.GetCommunicator().LocalMesh().Nodes, label, 0)

    def __prepare_standard_mesh(self):
        self.topographic_tool.SetFreeSurfaceMeshPosition()
        self.topographic_tool.AssignDryWetProperties(KM.FLUID)

    def __prepare_topographic_mesh(self):
        self.topographic_tool.RestoreDryWetProperties()
        self.topographic_tool.ApplyIdOffsetOnNodes()
        self.topographic_tool.ApplyIdOffsetOnElements()
        self.topographic_tool.ApplyIdOffsetOnConditions()
        self.topographic_tool.AssignSolidFluidProperties()
        self.topographic_tool.SetBathymetryMeshPosition()
    
    def __reset_standard_mesh(self):
        self.topographic_tool.UndoIdOffsetOnNodes()
        self.topographic_tool.UndoIdOffsetOnElements()
        self.topographic_tool.UndoIdOffsetOnConditions()
        self.topographic_tool.RestoreSolidFluidProperties()
        self.topographic_tool.ResetMeshPosition()
