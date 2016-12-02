# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping
# The following check is needed in case of an mpi-parallel compilation that is run serial
try:
    import KratosMultiphysics.mpi as KratosMPI
except:
    pass

class NonMatchingGridMapper:
    def __init__(self, model_part_origin, model_part_destination, input_settings):
        self.model_part_origin = model_part_origin
        self.model_part_destination = model_part_destination
        self.settings = input_settings

        self.CheckAndValidateJson() # Check json and set default parameters if necessary

        self.ReadAndCheckInterfaceModelParts()

        self.search_radius_map = self.settings["mapper_settings"]["search_radius"].GetDouble()
        self.search_radius_inverse_map = self.settings["mapper_settings"]["search_radius"].GetDouble()
        self.search_iterations = self.settings["mapper_settings"]["search_iterations"].GetInt()

        if (self.search_radius_map < 0):
            self.search_radius_map = KratosMapping.ComputeSearchRadius(self.interface_model_part_origin)
            self.search_radius_inverse_map = KratosMapping.ComputeSearchRadius(self.interface_model_part_destination)

        self.type_of_mapper = self.settings["mapper_settings"]["mapper_type"].GetString()

        # Initialization of Mapper, add new mappers here
        if (self.type_of_mapper == "NearestNeighbor"):
            self.InitializeNearestNeighborMapper()
        elif (self.type_of_mapper == "NearestElement"):
            self.InitializeNearestElementMapper()
        elif (self.type_of_mapper == "Barycentric"):
            self.InitializeBarycentricMapper()
        elif (self.type_of_mapper == "RBF"):
            self.InitializeRBFMapper()
        elif (self.type_of_mapper == "IterativeMortar"):
            self.InitializeIterativeMortarMapper()
        elif (self.type_of_mapper == "Mortar"):
            self.InitializeMortarMapper()
        elif (self.type_of_mapper == "IGA"):
            self.InitializeIGAMapper()
        else:
            raise ValueError("Specified mapper not valid / implemented")

        # TODO Interface PreProcess


    ##### Mapping Functions ###################################################
    # origin_variable:      Variable Name on Origin Interface
    # destination_variable: Variable Name on Destination Interface
    # add_values:           Values that are being mapped overwrite (add_values=False)
                          # or add to the destination values they are being mapped to (add_values=True)
    # sign_positive:        sign of the values that are being mapped. Changes the sign if set to "False"

    # This Function maps from Origin (interface_model_part_origin) to Destination (interface_model_part_destination)
    def Map(self, origin_variable, destination_variable, add_values = False, sign_positive = True):
        if (self.type_of_mapper == "Mortar" and self.two_mortar_mappers == True and self.mortar_reverse_mapping == True):
            #self.mapper_inverse_map.InverseMap(origin_variable, destination_variable, add_values, sign_positive)
            print("Map; inverse_map.InverseMap")
        else: # usual Case
            self.mapper_map.Map(origin_variable, destination_variable, add_values, sign_positive)
            print("Map; map.Map")

    # This Function maps from Destination (interface_model_part_destination) to Origin (interface_model_part_origin)
    def InverseMap(self, origin_variable, destination_variable, add_values = False, sign_positive = True):
        if (self.type_of_mapper == "Mortar" and (self.two_mortar_mappers == False or
                (self.two_mortar_mappers == True and self.mortar_reverse_mapping == True))):
            #self.mapper_map.InverseMap(origin_variable, destination_variable, add_values, sign_positive)
            print("InverseMap; map.InverseMap")
        else:
            self.mapper_inverse_map.Map(origin_variable, destination_variable, add_values, sign_positive)
            print("InverseMap; inverse_map.Map")


    ##### Auxilliary Functions ################################################
    def CheckAndValidateJson(self):
        if not (self.settings["mapper_settings"].Has("mapper_type")):
            raise ValueError("No mapper defined in json; \"mapper_type\" ")

        if not (self.settings["mapper_settings"].Has("interface_submodel_part_origin")):
            raise ValueError("No origin_submodel_part defined in json; \"interface_submodel_part_origin\"")

        if not (self.settings["mapper_settings"].Has("interface_submodel_part_destination")):
            raise ValueError("No destination_submodel_part defined in json; \"interface_submodel_part_destination\"")

        self.default_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type"                           : "",
            "interface_submodel_part_origin"        : "",
            "interface_submodel_part_destination"   : "",
            "search_radius"                         : -1.0,
            "search_iterations"                     : 10
        }""")

        # Overwrite the default settings with user-provided parameters
        self.settings["mapper_settings"].RecursivelyValidateAndAssignDefaults(self.default_parameters)

    def ReadAndCheckInterfaceModelParts(self):
        name_interface_submodel_part = self.settings["mapper_settings"]["interface_submodel_part_origin"].GetString()
        self.interface_model_part_origin = self.model_part_origin.GetSubModelPart(name_interface_submodel_part)

        name_interface_submodel_part = self.settings["mapper_settings"]["interface_submodel_part_destination"].GetString()
        self.interface_model_part_destination = self.model_part_destination.GetSubModelPart(name_interface_submodel_part)

        if (KratosMapping.ComputeNumberOfNodes(self.interface_model_part_origin) < 1 and
            KratosMapping.ComputeNumberOfConditions(self.interface_model_part_origin) < 1):
            raise ValueError("Neither nodes nor conditions found in the origin model part")
        if (KratosMapping.ComputeNumberOfNodes(self.interface_model_part_destination) < 1 and
            KratosMapping.ComputeNumberOfConditions(self.interface_model_part_destination) < 1):
            raise ValueError("Neither nodes nor conditions found in the destination model part")

        # TODO is the following necessary? => comes form the nonconformant mapper
        domain_size_origin = self.interface_model_part_origin.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        domain_size_destination = self.interface_model_part_destination.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        if (domain_size_origin != domain_size_destination):
            raise ValueError("Domain sizes from two model parts are not compatible")


    ##### Mapper Initialize Functions #########################################
    def InitializeNearestNeighborMapper(self):
        print("Python, NearestNeighborMapper initialized")
        self.mapper_map = KratosMapping.NearestNeighborMapper(self.interface_model_part_origin, self.interface_model_part_destination,
                                                self.search_radius_map, self.search_iterations)
        self.mapper_inverse_map = KratosMapping.NearestNeighborMapper(self.interface_model_part_destination, self.interface_model_part_origin,
                                                        self.search_radius_inverse_map, self.search_iterations)
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    def InitializeNearestElementMapper(self):
        print("Python, NearestElementMapper initialized")
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    def InitializeBarycentricMapper(self):
        raise ValueError("BarycentricMapper not implemented")
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    def InitializeRBFMapper(self):
        raise ValueError("RBFMapper not implemented")
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    def InitializeIterativeMortarMapper(self):
        print("Python, IterativeMortarMapper initialized")

        # Initialize and set mapper - specific parameters
        default_parameters_iterative_mortar = KratosMultiphysics.Parameters( """
        {
            "convergence_tolerance"     : 1e-6,
            "convergence_iterations"    : 100
        }  """ )

        self.default_parameters.AddValue("iterative_mortar", default_parameters_iterative_mortar)

        # Overwrite the default settings with user-provided parameters
        self.settings["mapper_settings"].RecursivelyValidateAndAssignDefaults(self.default_parameters);
        # print(self.settings.PrettyPrintJsonString())
        convergence_tolerance = self.settings["mapper_settings"]["iterative_mortar"]["convergence_tolerance"].GetDouble()
        convergence_iterations = self.settings["mapper_settings"]["iterative_mortar"]["convergence_iterations"].GetInt()

        self.mapper_map = KratosMapping.IterativeMortarMapper(self.interface_model_part_origin, self.interface_model_part_destination,
                                                self.search_radius_map, self.search_iterations,
                                                convergence_tolerance, convergence_iterations)
        self.mapper_inverse_map = KratosMapping.IterativeMortarMapper(self.interface_model_part_destination, self.interface_model_part_origin,
                                                        self.search_radius_inverse_map, self.search_iterations,
                                                        convergence_tolerance, convergence_iterations)
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    def InitializeMortarMapper(self):
        print("Python, MortarMapper initialized")

        # Initialize and set mapper - specific parameters
        default_parameters_mortar = KratosMultiphysics.Parameters( """
        {   "two_mappers"       : false,
            "reverse_mapping"   : false,
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            }
        }  """ )

        self.default_parameters.AddValue("mortar", default_parameters_mortar)

        # Overwrite the default settings with user-provided parameters
        self.settings["mapper_settings"].RecursivelyValidateAndAssignDefaults(self.default_parameters);
        # print(self.settings.PrettyPrintJsonString())

        if (self.settings["mapper_settings"]["mortar"]["two_mappers"].GetBool()):
            self.two_mortar_mappers = True
            self.mortar_reverse_mapping = self.settings["mapper_settings"]["mortar"]["reverse_mapping"].GetBool()
        else:
            self.two_mortar_mappers = False
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


    def InitializeIGAMapper(self):
        raise ValueError("IGAMapper not implemented")
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
