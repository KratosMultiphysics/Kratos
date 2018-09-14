from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
import co_simulation_tools as tools

##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change any thing in this class.
#
# This Class servers as a base class for all the Input-output methods to be implemented
class CoSimulationBaseIO(object):
    ## The constructor
    #
    #  @param self            The object pointer.
    def __init__(self):
        self.echo_level = 0

    ## ImportData :  used to import data from other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param self            The object pointer.
    #  @param data_config     python dictionary : configuration of the data to be imported.
    #                                             data will be imported into this dictionary
    #  @param from_solver     python object : The solver from which data is to be imported.
    def ImportData(self, data_config, from_solver):
        """
        if(from_client not None):
            # exchange data from python cosim solver
        else:
            # import data from remote solver
        """
        raise NotImplementedError(tools.bcolors.FAIL + "From BaseIO : The method ImportData is not implemented in the IO class!" + tools.bcolors.ENDC)

    ## ImportMesh :  used to import mesh from other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param self            The object pointer.
    #  @param mesh_conig      python dictionary : configuration of the mesh to be imported.
    #                                             mesh will be imported into this dictionary.
    #  @param from_solver     python object : The solver from which mesh is to be imported.
    def ImportMesh(self, mesh_conig, from_solver):
        """
        if(from_client not None):
            # exchange mesh from python cosim solver
        else:
            # import mesh from remote solver
        """
        raise NotImplementedError(tools.bcolors.FAIL + "From BaseIO : The method ImportMesh is not implemented in the IO class!" + tools.bcolors.ENDC)

    ## ExportData :  used to export data to other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param self            The object pointer.
    #  @param data_config     python dictionary : configuration of the mesh to be exported.
    #                                             also contains the data to export.
    #  @param to_solver       python object : The solver to which mesh is to be exported.
    def ExportData(self, data_config, to_solver):
        """
        if(to_client not None): # IMPORTANT : exchangig data between python cosim solvers should be avoided here.
            # put data on to python cosim solver.
        else:
            # export the given data to the remote solver
        """
        raise NotImplementedError(tools.bcolors.FAIL + "From BaseIO : The method ExportData is not implemented in the IO class!" + tools.bcolors.ENDC)

    ## ExportMesh :  used to export mesh to other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param self            The object pointer.
    #  @param mesh_conig      python dictionary : configuration of the mesh to be exported.
    #                                             also contains the mesh data to export.
    #  @param to_solver       python object : The solver to which mesh is to be exported.
    def ExportMesh(self, mesh_conig, to_solver):
        """
        if(to_client not None): # IMPORTANT : exchangig mesh between python cosim solvers should be avoided here.
            # put mesh on to python cosim solver.
        else:
            # export the given mesh to the remote solver
        """
        raise NotImplementedError(tools.bcolors.FAIL + "From BaseIO : The method ExportMesh is not implemented in the IO class!" + tools.bcolors.ENDC)

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
