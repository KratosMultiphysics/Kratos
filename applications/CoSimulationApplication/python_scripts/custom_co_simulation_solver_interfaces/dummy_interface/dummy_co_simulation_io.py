# co simulation imports
from co_simulation_base_io import CoSimulationBaseIO
import co_simulation_tools as tools


def Create(settings):
    return DummyCoSimulationIO(settings)


class DummyCoSimulationIO(CoSimulationBaseIO):
    ## The constructor
    #
    #  @param self            The object pointer.
    def __init__(self, settings):
        self.echo_level = 0
        self.settings = settings

    ## ImportData :  used to import data from other clients
    #
    #  IMPORTANT :: Import and Export by default will import the mesh and data from the remote solver.
    #
    #  @param self            The object pointer.
    #  @param data_config     python dictionary : configuration of the data to be imported.
    #                                             data will be imported into this dictionary
    #  @param from_client     python object : The client from which data is to be imported.
    def ImportData(self, data_config, from_client=None):
        if(from_client != None):
            # exchange data from python cosim solver
            pass
        else:
            # import data from remote solver
            pass


    ## ImportMesh :  used to import mesh from other clients
    #
    #  @param self            The object pointer.
    #  @param mesh_conig      python dictionary : configuration of the mesh to be imported.
    #                                             mesh will be imported into this dictionary.
    #  @param from_client     python object : The client from which mesh is to be imported.
    def ImportMesh(self, mesh_conig, from_client=None):
        if(from_client != None):
            # exchange mesh from python cosim solver
            pass
        else:
            # import mesh from remote solver
            pass

    ## ExportData :  used to export data to other clients
    #
    #  @param self            The object pointer.
    #  @param data_config     python dictionary : configuration of the mesh to be exported.
    #                                             also contains the data to export.
    #  @param to_client       python object : The client to which mesh is to be exported.
    def ExportData(self, data_config, to_client=None):
        if(to_client != None): # IMPORTANT : exchanging data between python cosim solvers should be avoided here.
            # put data on to python cosim solver.
            pass
        else:
            # export the given data to the remote solver
            pass

    ## ExportMesh :  used to export mesh to other clients
    #
    #  @param self            The object pointer.
    #  @param mesh_conig      python dictionary : configuration of the mesh to be exported.
    #                                             also contains the mesh data to export.
    #  @param to_client       python object : The client to which mesh is to be exported.
    def ExportMesh(self, mesh_conig, to_client=None):
        if(to_client != None): # IMPORTANT : exchanging mesh between python cosim solvers should be avoided here.
            # put mesh on to python cosim solver.
            pass
        else:
            # export the given mesh to the remote solver
            pass
        pass

    ## Sets the echo level of for this IO object. Used for output of information during CoSimulation
    #
    #  @param self            The object pointer.
    #  @param level           int : echo level
    def SetEchoLevel(self, level):
        self.echo_level = level

    ## Prints the information about the current IO. The derived classes should implement it.
    #
    #  @param self            The object pointer.
    def PrintInfo(self):
        print(tools.bcolors.WARNING + "!!!WARNING!!! Calling PrintInfo from base IO class!" + tools.bcolors.ENDC)

    ## Check : Checks the current IO against faulty settings. Derived classes may implement this function.
    #
    #  @param self            The object pointer.
    def Check(self):
        print(tools.bcolors.WARNING + "!!!WARNING!!! Calling Check from base IO class!" + tools.bcolors.ENDC)
