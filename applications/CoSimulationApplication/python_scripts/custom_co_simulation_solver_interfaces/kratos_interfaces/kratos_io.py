try:
    import KratosMultiphysics
    # Check that applications were imported in the main script
    KratosMultiphysics.CheckRegisteredApplications("CoSimulationApplication")
    import KratosMultiphysics.CoSimulationApplication as CoSimApp
except ModuleNotFoundError:
    print(tools.bcolors.FAIL + 'KRATOS is not available ! Please ensure that Kratos is available for usage !'+ tools.bcolors.ENDC)
    exit()

from base_co_simulation_classes.co_simulation_base_io import CoSimulationBaseIO
import co_simulation_tools as tools

# Other imports
import os

def Create(custom_settings):
    return KratosIo(custom_settings)

class KratosIo(CoSimulationBaseIO):
    def __init__(self, custom_settings):
        self.settings = custom_settings
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
            {
                "io_type":"kratos_io",
                "settings":{
                }
            }
        """)
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.mappers = {}

        super(KratosIo, self).__init__(self.settings)
        ### Constructing the IO for this solver

    ## ImportData :  used to import data from other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param self            The object pointer.
    #  @param data_config     python dictionary : configuration of the data to be imported.
    #                                             data will be imported into this dictionary
    #  @param from_solver     python object : The solver from which data is to be imported.
    def ImportData(self, data_config, from_solver):
        if(from_solver):
            # exchange data from python cosim solver
            if(data_config.Has("mapper_settings")):
                origin_geo = data_config["origin_data_config"]["geometry_name"].GetString()
                origin_var = data_config["origin_data_config"]["name"].GetString()
                dest_geo = data_config["geometry_name"].GetString()
                dest_var = data_config["name"].GetString()
                if(self.HasMapper((origin_geo,dest_geo))):
                    mapper = self.GetMapper((origin_geo,dest_geo))
                elif(self.HasMapper((dest_geo,origin_geo))):
                    mapper = self.GetMapper((dest_geo,origin_geo))

                flags = data_config["mapper_settings"]
                mapper.Map(origin_var, dest_var, flags)
            pass
        else:
            # import data from remote solver
            pass


    ## ImportMesh :  used to import mesh from other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param self            The object pointer.
    #  @param mesh_config      python dictionary : configuration of the mesh to be imported.
    #                                             mesh will be imported into this dictionary.
    #  @param from_solver     python object : The solver from which mesh is to be imported.
    def ImportMesh(self, mesh_config, from_solver):
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
        if(to_solver): # IMPORTANT : exchanging data between python cosim solvers should be avoided here.
            # put data on to python cosim solver.
            print('print from kratos io')
            pass
        else:
            # export the given data to the remote solver
            pass

    ## ExportMesh :  used to export mesh to other solvers
    #                Follow EXAMPLE implementation below.
    #
    #  @param self            The object pointer.
    #  @param mesh_conig      python dictionary : configuration of the mesh to be exported.
    #                                             also contains the mesh data to export.
    #  @param to_solver       python object : The solver to which mesh is to be exported.
    def ExportMesh(self, mesh_conig, to_solver):
        """
        if(to_client not None): # IMPORTANT : exchanging mesh between python cosim solvers should be avoided here.
            # put mesh on to python cosim solver.
        else:
            # export the given mesh to the remote solver
        """
        raise NotImplementedError(tools.bcolors.FAIL + "From BaseIO : The method ExportMesh is not implemented in the IO class!" + tools.bcolors.ENDC)


    def HasMapper(self, mapper_tuple):
        return True

    def GetMapper(self, mapper_tuple):
        return self.mappers[(mapper_tuple)]