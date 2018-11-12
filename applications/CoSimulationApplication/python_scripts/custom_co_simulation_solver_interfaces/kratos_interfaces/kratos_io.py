try:
    import KratosMultiphysics
    # Check that applications were imported in the main script
    KratosMultiphysics.CheckRegisteredApplications("CoSimulationApplication")
    import KratosMultiphysics.CoSimulationApplication as CoSimApp
    import KratosMultiphysics.MappingApplication as MappingApp
except ModuleNotFoundError:
    print(tools.bcolors.FAIL + 'KRATOS is not available ! Please ensure that Kratos is available for usage !'+ tools.bcolors.ENDC)
    exit()

from base_co_simulation_classes.co_simulation_base_io import CoSimulationBaseIO
import co_simulation_tools as tools

# Other imports
import os

def Create(model, custom_settings):
    return KratosIo(model, custom_settings)

class KratosIo(CoSimulationBaseIO):
    def __init__(self, model, custom_settings):
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
        self.mapper_flags = {
            "add_values" : MappingApp.Mapper.ADD_VALUES,
            "swap_sign" : MappingApp.Mapper.SWAP_SIGN,
            "conservative" : MappingApp.Mapper.CONSERVATIVE,
        }

        super(KratosIo, self).__init__(model, self.settings)
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
                origin_geo_name = data_config["origin_data_config"]["geometry_name"].GetString()
                origin_var_name = data_config["origin_data_config"]["name"].GetString()
                dest_geo_name = data_config["geometry_name"].GetString()
                dest_var_name = data_config["name"].GetString()
                if(self.HasMapper((origin_geo_name,dest_geo_name))):
                    mapper = self.GetMapper((origin_geo_name,dest_geo_name))
                elif(self.HasMapper((dest_geo_name,origin_geo_name))):
                    mapper = self.GetMapper((dest_geo_name,origin_geo_name))
                else: # Create the mapper
                    mapper_settings = KratosMultiphysics.Parameters("""{
                                                                        "mapper_type" : ""
                                                                    }""")
                    mapper_settings["mapper_type"].SetString(data_config["mapper_settings"]["type"].GetString())
                    origin_geo = from_solver.model[origin_geo_name]
                    dest_geo = self.model[dest_geo_name]
                    print("origin_geo :: ",origin_geo)
                    print("dest_geo :: ",dest_geo)

                    mapper = self.SetupMapper(origin_geo, dest_geo, mapper_settings)

                flags = data_config["mapper_settings"]["flags"]
                origin_var = KratosMultiphysics.KratosGlobals.GetVariable(origin_var_name)
                dest_var = KratosMultiphysics.KratosGlobals.GetVariable(dest_var_name)

                set_flags = KratosMultiphysics.Flags()
                if data_config["mapper_settings"].Has("flags"):
                    num_flags = flags.size()
                    for i in range(num_flags):
                        flag_name = flags[i].GetString()
                        set_flags |= self.mapper_flags[flag_name]

                mapper.Map(origin_var, dest_var, set_flags)
            else:
                pass ## We can copy one to one instead of mapping here.
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
            print('############################## print from kratos io')
            adsfsdf
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
        return mapper_tuple in self.mappers.keys()

    def GetMapper(self, mapper_tuple):
        return self.mappers[(mapper_tuple)]

    def SetupMapper(self, geo_origin, geo_destination, mapper_settings):
        mapper = MappingApp.MapperFactory.CreateMapper(geo_origin, geo_destination, mapper_settings)
        self.mappers[(geo_origin.Name, geo_destination.Name)] = mapper
        return mapper