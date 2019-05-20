from __future__ import print_function, absolute_import, division

# Other imports
import KratosMultiphysics.CoSimulationApplication.custom_solver_interfaces.co_simulation_io_factory as cs_io_factory
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure
import collections

##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change any thing in this class.
#
# This Class servers as a base class for all the Solver interfaces which will be implemented
class CoSimulationBaseSolver(object):
    ## Constructor :  The base class for the CoSimulation Solver interfaces
    #                  Constructor of the Base-Solver interface
    #                  Deriving classes should call it in their constructors
    #
    #  @param self                      The object pointer.
    #  @param cs_solver_settings     python dictionary : with the solver settings.
    def __init__(self, solver_name, cs_solver_settings):
        default_settings = cs_data_structure.Parameters("""
        {
            "solver_type" : "",
            "io_settings": {},
            "settings" : {},
            "data" : [],
            "echo_level" : 0
        }
        """)
        self.name = solver_name
        self.cs_solver_settings = cs_solver_settings
        self.cs_solver_settings.ValidateAndAssignDefaults(default_settings)
        self.SetEchoLevel(self.cs_solver_settings["echo_level"].IsInt())
        self.data_map = self._GetDataMap()
        self.geo_names = self._GetGeometryNames()
        self.model = cs_data_structure.Model() # Where all the co-simulation meshes are stored.
        # This is the map of all the geometries that a solver can have
        self.geometry_map = {}
        self.io_is_initialized = False


    ## Initialize : Initialize function for the solver class. Necessary
    #               initialization of the variables and objects to be done here.
    #  @param self                      The object pointer. especially the geometries
    def Initialize(self):
        # Initialize IO
        if not self.io_is_initialized:
            self.InitializeIO(self.echo_level)

        # Initialize data (and geometries)

    ## InitializeIO : Initialize the IO class for this solver.
    #                   usually a particular type of solver has a particular default IO type
    #  @param self                      The object pointer.
    #  @param io_echo_level             int : echo level for the io to be initialized.
    def InitializeIO(self, io_echo_level=0):
        solver_name = self.name
        if self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is already initialized!')

        self.io = cs_io_factory.CreateIO(self.model, self.cs_solver_settings["io_settings"])
        self.io.SetEchoLevel(io_echo_level)
        self.io_is_initialized = True

    ## Finalize :  Initialize function for the solver class.
    #               finalization of the variables and objects to be done here.
    #  @param self                      The object pointer.
    def Finalize(self):
        pass

    ## AdvanceInTime :  This function will advance the current solver to the given current_time
    #
    #  @param self                      The object pointer.
    #  @current_time                    The time to which the current solver should advance
    def AdvanceInTime(self, current_time):
        pass

    ## Predict : Predict the solution of the next solution step.
    #
    #  @param self            The object pointer.
    def Predict(self):
        pass

    ## InitializeSolutionStep : Called once in the beginning of the solution step
    #
    #  @param self            The object pointer.
    def InitializeSolutionStep(self):
        for name, data in self.data_map.items():
            for filter in data.filters:
                filter.InitializeSolutionStep()

    ## FinalizeSolutionStep : Called once at the end of the solution step
    #
    #  @param self            The object pointer.
    def FinalizeSolutionStep(self):
        for name, data in self.data_map.items():
            for filter in data.filters:
                filter.FinalizeSolutionStep()

    ## InitializeCouplingIteration : Called once in the beginning of the coupled iteration
    #
    #  @param self            The object pointer.
    def InitializeCouplingIteration(self):
        for name, data in self.data_map.items():
            for filter in data.filters:
                filter.InitializeCouplingIteration()

    ## FinalizeCouplingIteration : Called once at the end of the coupled iteration
    #
    #  @param self            The object pointer.
    def FinalizeCouplingIteration(self):
        for name, data in self.data_map.items():
            for filter in data.filters:
                filter.FinalizeCouplingIteration()

    ## OutputSolutionStep : Called once at the end of the solution step.
    #                       The output of the solvers and / or cosimulation output
    #                       can be performed in this function.
    #
    #  @param self            The object pointer.
    def OutputSolutionStep(self):
        pass

    ## SolveSolutionStep : This function implements the coupling workflow
    #
    #  @param self            The object pointer.
    def SolveSolutionStep(self):
        pass

    ## GetDataConfig : This function gets the definition of data in this solver
    #
    #  @param self            The object pointer.
    def GetInterfaceData(self, data_name):
        if data_name in self.data_map.keys():
            return self.data_map[data_name]
        else:
            raise Exception(cs_tools.bcolors.FAIL+ "Requested data field " + data_name + " does not exist in the solver "+self.name+cs_tools.bcolors.ENDC)

    ## ImportCouplingInterfaceData : This function imports the requested data from
    #               from_client
    #
    #  @param self            The object pointer.
    #  @param data_name       string : Name of the data to be imported from from_client
    #  @param from_client     python obj : The client from which data_name has to be imported
    #                         Default is None that means, the solver imports the mesh from itself,
    #                         that is the actual solver for which this acts as an alias
    def ImportCouplingInterfaceData(self, data_object, from_client=None):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ImportCouplingInterfaceData(data_object, from_client) ## Mapping happens inside here
        ## Once the mapping is done here, we apply all the filters before making the data usable
        data_object.ApplyFilters()

    ## ImportCouplingInterface : This function imports the requested surface/volume
    #               mesh from from_client
    #
    #  @param self            The object pointer.
    #  @param mesh_name       string : Name of the mesh to be imported from from_client
    #  @param from_client     python obj : The client from which data_name has to be imported.
    #                         Default is None that means, the solver imports the mesh from itself,
    #                         that is the actual solver for which this acts as an alias
    def ImportCouplingInterface(self, mesh_name, from_client=None):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ImportCouplingInterface(mesh_name, from_client)

    ## ExportCouplingInterfaceData : This function exports the requested data to
    #               to_client
    #
    #  @param self            The object pointer.
    #  @param to_client       The client to which the data is to be exported
    #                         Default is None that means, the solver imports the mesh from itself,
    #                         that is the actual solver for which this acts as an alias
    def ExportCouplingInterfaceData(self, data_object, to_client=None):
        ## Before we export the data, we apply all the filters so the solver importing can readily use it.
        data_object.ApplyFilters()
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ExportCouplingInterfaceData(data_object, to_client)

    ## ExportCouplingInterface : This function exports the requested surface/volume
    #               to to_client
    #
    #  @param self            The object pointer.
    #  @param to_client       The client to which the Mesh is to be exported
    #                         Default is None that means, the solver imports the mesh from itself,
    #                         that is the actual solver for which this acts as an alias
    def ExportCouplingInterface(self, mesh_name, to_client=None):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ExportCouplingInterface(mesh_name, to_client)

    ## GetDeltaTime : Function to obtain the time step of this solver
    #
    #  @param self            The object pointer.
    def GetDeltaTime(self):
        raise NotImplementedError(cs_tools.bcolors.FAIL + "CoSimulationBaseSolver : Calling GetDeltaTime function from base Co-Simulation Solver class!" + cs_tools.bcolors.ENDC)

    ## PrintInfo : Function to display information about the
    #              specifics of the coupled solver.
    #
    #  @param self            The object pointer.
    def PrintInfo(self):
        raise NotImplementedError(cs_tools.bcolors.FAIL + "CoSimulationBaseSolver : Calling PrintInfo function from base Co-Simulation Solver class!" + cs_tools.bcolors.ENDC)

    ## SetEchoLevel : Function to set the Echo level for this solver
    #
    #  @param self            The object pointer.
    def SetEchoLevel(self, level):
        self.echo_level = level

    ## Check : Called once at the beginning of simulation.
    #          Necessary setup of the solver can be verified here.
    #
    #  @param self            The object pointer.
    def Check(self):
        raise NotImplementedError(cs_tools.bcolors.FAIL + "CoSimulationBaseSolver : Calling Check function from base Co-Simulation Solver class!" + cs_tools.bcolors.ENDC)

    ## _GetIOType : Private Function to obtain the type of IO this solver has as default
    #
    #  @param self            The object pointer.
    def _GetIOType(self):
        raise NotImplementedError(cs_tools.bcolors.FAIL + "CoSimulationBaseSolver : Calling _GetIOName function from base Co-Simulation Solver class!" + cs_tools.bcolors.ENDC)

    ## _GetDataMap : Private Function to obtain the map of data objects
    #
    #  @param self            The object pointer.
    def _GetDataMap(self):
        data_map = collections.OrderedDict()
        for data_conf in self.cs_solver_settings["data"]:
            data_name = data_conf["name"].GetString()
            data_object = cs_tools.CouplingInterfaceData(data_conf, self)
            data_map[data_name] = data_object
        return data_map

    def _GetGeometryNames(self):
        geometry_name_list = []
        for name, data in self.data_map.items():
            geometry_name = data.geometry_name
            if geometry_name not in geometry_name_list:
                geometry_name_list.append(geometry_name)
        return geometry_name_list
