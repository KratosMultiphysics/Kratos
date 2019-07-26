from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change any thing in this class.
#
# This Class servers as a base class for all the Input-output methods to be implemented
class CoSimulationIO(object):
    ## The constructor
    #
    def __init__(self, model, settings):
        self.echo_level = 0
        self.settings = settings
        self.model = model

    ## ImportCouplingInterfaceData :  used to import data from other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param data_object     python dictionary : object of the data to be imported.
    #                                             data will be imported into this dictionary
    #  @param from_solver     python object : The solver from which data is to be imported.
    def ImportCouplingInterfaceData(self, data_object, from_solver=None):
        """
        if(from_client not None):
            # exchange data from python cosim solver
        else:
            # import data from remote solver
        """
        raise NotImplementedError(cs_tools.darkred("From BaseIO : The method ImportCouplingInterfaceData is not implemented in the IO class!"))

    ## ImportCouplingInterface :  used to import mesh from other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param mesh_config      python dictionary : configuration of the mesh to be imported.
    #  @param from_solver     python object : The solver from which mesh is to be imported.
    def ImportCouplingInterface(self, mesh_config, from_solver=None):
        """
        if(from_client not None):
            # exchange mesh from python cosim solver
        else:
            # import mesh from remote solver
        """
        raise NotImplementedError(cs_tools.darkred("From BaseIO : The method ImportCouplingInterface is not implemented in the IO class!"))

    ## ExportCouplingInterfaceData :  used to export data to other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param data_object     python dictionary : object of the mesh to be exported.
    #                                             also contains the data to export.
    #  @param to_solver       python object : The solver to which mesh is to be exported.
    def ExportCouplingInterfaceData(self, data_object, to_solver=None):
        """
        if(to_client not None): # IMPORTANT : exchanging data between python cosim solvers should be avoided here.
            # put data on to python cosim solver.
        else:
            # export the given data to the remote solver
        """
        raise NotImplementedError(cs_tools.darkred("From BaseIO : The method ExportCouplingInterfaceData is not implemented in the IO class!"))

    ## ExportCouplingInterface :  used to export mesh to other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param mesh_config      python dictionary : configuration of the mesh to be exported.
    #                                             also contains the mesh data to export.
    #  @param to_solver       python object : The solver to which mesh is to be exported.
    def ExportCouplingInterface(self, mesh_config, to_solver=None):
        """
        if(to_client not None): # IMPORTANT : exchanging mesh between python cosim solvers should be avoided here.
            # put mesh on to python cosim solver.
        else:
            # export the given mesh to the remote solver
        """
        raise NotImplementedError(cs_tools.darkred("From BaseIO : The method ExportCouplingInterface is not implemented in the IO class!"))

    ## Prints the information about the current IO. The derived classes should implement it.
    #
    def PrintInfo(self):
        cs_tools.PrintWarning(cs_tools.red("WARNING", "Calling PrintInfo from base IO class!"))

    ## Check : Checks the current IO against faulty settings. Derived classes may implement this function.
    #
    def Check(self):
        cs_tools.PrintWarning(cs_tools.red("WARNING", "Calling Check from base IO class!"))
