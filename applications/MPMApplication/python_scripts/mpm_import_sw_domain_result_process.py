import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.MPMApplication as KratosMPM
# import KratosMultiphysics.MappingApplication as Mapping
from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput
from KratosMultiphysics.HDF5Application import import_model_part_from_hdf5_process
from KratosMultiphysics.HDF5Application import single_mesh_temporal_input_process
from os.path import commonprefix
from pathlib import Path
import re

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MpmImportSWDomainResultProcess(model, settings["Parameters"])

class MpmImportSWDomainResultProcess(KM.Process):
    """MpmImportSWDomainResultProcess

    Read the Shallow Water nodal values from an HDF5 file. Then generate and set them as initial conditions for the material points.
    """

    def GetDefaultParameters(self):
        default_parameters = KM.Parameters("""{
            "input_model_part_name"       : "sw_input_model_part",
            "read_historical_database"    : true,
            "interval"                    : [0.0,"End"],
            "minimum_water_depth"         : 0.0,
            "file_settings"               : {}
        }""")

        return default_parameters

    def __init__(self, model, settings):
        """The constructor of the MpmImportSWDomainResultProcess."""

        KM.Process.__init__(self)
        self.model = model
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.input_model_part = model.CreateModelPart(self.settings["input_model_part_name"].GetString())


        self.interval = KM.IntervalUtility(self.settings)

        self.hdf5_import = import_model_part_from_hdf5_process.Factory(self._CreateHDF5ParametersImport(), model)
        self.hdf5_process = single_mesh_temporal_input_process.Factory(self._CreateHDF5Parameters(), model)
        # self._GetInputTimes(self.settings['file_settings'])

        # temporary hacks to get mpm model parts
        self.grid_model_part = self.model["Background_Grid"]
        self.material_point_model_part = self.model["MPM_Material"]


    def Check(self): #
        '''Check the processes.'''
        self.hdf5_import.Check()
        self.hdf5_process.Check()
        # free_surface_is_present = False
        # height_is_present = False
        # for variable in self.variables:
        #     if variable == SW.FREE_SURFACE_ELEVATION:
        #         free_surface_is_present = True
        #     if variable == SW.HEIGHT:
        #         height_is_present = True
        # if free_surface_is_present and not height_is_present:
        #     self.variables.append(SW.HEIGHT)


    def ExecuteInitialize(self):#
        '''Read the input_model_part and set the variables.'''
        self.hdf5_import.ExecuteInitialize()
        print("andihehe", self.model)

        # self.hdf5_process.ExecuteInitialize()
        self.input_model_part.ProcessInfo = self.grid_model_part.ProcessInfo
        self.input_model_part.SetProperties(self.model["Initial_MPM_Material"].Properties)
        print("andihehe", self.input_model_part.GetProperties()[1])

        # for element in self.input_model_part.Elements:
        #     element.SetProperties(self.model["Initial_MPM_Material"].GetProperties(1))
        self.hdf5_process.ExecuteInitializeSolutionStep()
        KratosMPM.GenerateMaterialPointElementFromSwModel(self.grid_model_part,
                                                          self.input_model_part,
                                                          self.material_point_model_part,
                                                          False,
                                                          self.settings["minimum_water_depth"].GetDouble() )


    # def _GetInputTimes(self, file_settings):
    #     # Get all the file names
    #     file = Path(file_settings["file_name"].GetString())
    #     directory = file.parent
    #     file_names = [str(f) for f in directory.glob("*.h5")]

    #     if len(file_names) == 0:
    #         msg = self.__class__.__name__ + ": The specified path is empty or does not exist: '{}'"
    #         raise Exception(msg.format(directory))

    #     # Extract the time stamp from the end of each file name
    #     self.times = []
    #     pattern = re.compile(r'(\d+(?:\.\d+)?)\.h5$')

    #     for f in file_names:
    #         match = pattern.search(Path(f).name)

    #         # Skip files without a time stamp (e.g. sw_output_model_part.h5)
    #         if match is None:
    #             continue

    #         self.times.append(float(match.group(1)))

    #     if len(self.times) == 0:
    #         raise Exception(
    #             self.__class__.__name__ + ": No time-stamped .h5 files were found."
    #         )

    #     self.times.sort()

    def _CreateHDF5ParametersImport(self): #
        hdf5_settings = KM.Parameters()
        hdf5_settings.AddValue("model_part_name", self.settings["input_model_part_name"])
        import_file_settings = KM.Parameters("""{"file_name" : "hdf5_output/<model_part_name>-0.0000.h5"}""")
        hdf5_settings.AddValue("file_settings", import_file_settings)
        data_settings = KM.Parameters("""{"list_of_variables" : []}""")
        if self.settings["read_historical_database"].GetBool():
            hdf5_settings.AddValue("nodal_solution_step_data_settings", data_settings)
        else:
            hdf5_settings.AddValue("nodal_data_value_settings", data_settings)
        hdf5_process_settings = KM.Parameters()
        hdf5_process_settings.AddValue("Parameters", hdf5_settings)
        return hdf5_process_settings

    def _CreateHDF5Parameters(self): #
        hdf5_settings = KM.Parameters()
        hdf5_settings.AddValue("model_part_name", self.settings["input_model_part_name"])
        hdf5_settings.AddValue("file_settings", self.settings["file_settings"])
        data_settings = KM.Parameters("""{"list_of_variables" : ["MOMENTUM","VELOCITY","HEIGHT","VERTICAL_VELOCITY","TOPOGRAPHY"]}""")
        if self.settings["read_historical_database"].GetBool():
            hdf5_settings.AddValue("nodal_solution_step_data_settings", data_settings)
        else:
            hdf5_settings.AddValue("nodal_data_value_settings", data_settings)
        hdf5_process_settings = KM.Parameters()
        hdf5_process_settings.AddValue("Parameters", hdf5_settings)
        return hdf5_process_settings
