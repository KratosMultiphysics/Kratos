import KratosMultiphysics
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.MeshingApplication as KratosMesh

from KratosMultiphysics.GeodataProcessingApplication.geo_processor import GeoProcessor

class GeoModel( GeoProcessor ):

    def __init__( self ):
        super(GeoModel, self).__init__()

        # self.HasModelPart = False
        self.HasCfdModelPart = False

        # we initialize json file
        self._parameter_initialization()


    def GenerateCfdModelPart(self):

        # current_model = KratosMultiphysics.Model()
        current_model = self.ModelPart.GetModel()
        model_part_name = "CFD_model_part"

        if current_model.HasModelPart(model_part_name):
            # clear the existing model part (to be sure)
            self.CfdModelPart = current_model.GetModelPart(model_part_name)
            self.CfdModelPart.Elements.clear()
            self.CfdModelPart.Conditions.clear()
            self.CfdModelPart.Nodes.clear()
        else:
            self.CfdModelPart = current_model.CreateModelPart(model_part_name)

        # self.CfdModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        # self.CfdModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

        # we set the DENSITY and DYNAMIC_VISCOSITY values
        prop = self.CfdModelPart.GetProperties()[0]
        prop.SetValue(KratosMultiphysics.DENSITY, 1)
        prop.SetValue(KratosMultiphysics.DYNAMIC_VISCOSITY, 0.002)

        self.HasCfdModelPart = True


    def GetGeoCfdModelPart(self):
        
        if self.HasCfdModelPart:
            return self.CfdModelPart
        else:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "No CFD model part can be returned")


    def FillPartsFluid(self, elem_sub_model_name):
        # we fill element into Parts_Fluid sub model part
        # we perform a check if the nodes and the elements are already into the main model part
        # and, if not, we create it and we added it into sub model part

        # if not self.HasCfdModelPart:
        #     KratosMultiphysics.Logger.PrintWarning("GeoModel", "CFD model part has to be set, first.")
        #     return
        
        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "A model part has to be set, first.")
            return

        # if not self.CfdModelPart.HasSubModelPart("Parts_Fluid"):
        if not self.ModelPart.HasSubModelPart("Parts_Fluid"):
            # self.CfdModelPart.CreateSubModelPart("Parts_Fluid")
            self.ModelPart.CreateSubModelPart("Parts_Fluid")
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Parts_Fluid\" created right now!")

        # KratosGeo.FillCfdModelpartUtilities(self.CfdModelPart).FillPartsFluid(self.ModelPart, elem_sub_model_name)
        KratosGeo.FillCfdModelpartUtilities(self.ModelPart).FillPartsFluid(elem_sub_model_name)


    def FillInlet(self, cond_sub_model_name):
        # we fill conditions into Inlet sub model part
        # we perform a check if the nodes and the condiions are already into the main model part
        # and, if not, we create it and we added it into sub model part

        if not self.HasCfdModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "CFD model part has to be set, first.")
            return

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "A model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("Inlet"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Inlet\" created right now!")
            self.CfdModelPart.CreateSubModelPart("Inlet")
        
        KratosGeo.FillCfdModelpartUtilities(self.CfdModelPart).FillInlet(self.ModelPart, cond_sub_model_name)


    def FillOutlet(self, cond_sub_model_name):
        # we fill conditions into Outlet sub model part
        # we perform a check if the nodes and the condiions are already into the main model part
        # and, if not, we create it and we added it into sub model part
        
        if not self.HasCfdModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "CFD model part has to be set, first.")
            return

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "A model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("Outlet"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Outlet\" created right now!")
            self.CfdModelPart.CreateSubModelPart("Outlet")
        
        KratosGeo.FillCfdModelpartUtilities(self.CfdModelPart).FillOutlet(self.ModelPart, cond_sub_model_name)


    def FillSlip(self, cond_sub_model_name):
        # we fill conditions into Slip sub model part
        # we perform a check if the nodes and the condiions are already into the main model part
        # and, if not, we create it and we added it into sub model part
        
        # if not self.HasCfdModelPart:
        #     KratosMultiphysics.Logger.PrintWarning("GeoModel", "CFD model part has to be set, first.")
        #     return

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "A model part has to be set, first.")
            return

        # if not self.CfdModelPart.HasSubModelPart("Slip"):
        if not self.ModelPart.HasSubModelPart("Slip"):
            # self.CfdModelPart.CreateSubModelPart("Slip")
            self.ModelPart.CreateSubModelPart("Slip")
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Slip\" created right now!")
        
        # KratosGeo.FillCfdModelpartUtilities(self.CfdModelPart).FillSlip(self.ModelPart, cond_sub_model_name)
        KratosGeo.FillCfdModelpartUtilities(self.ModelPart).FillSlip(cond_sub_model_name)


    def FillNoslip(self, cond_sub_model_name):
        # we fill conditions into NoSlip sub model part
        # we perform a check if the nodes and the condiions are already into the main model part
        # and, if not, we create it and we added it into sub model part
        
        # if not self.HasCfdModelPart:
        #     KratosMultiphysics.Logger.PrintWarning("GeoModel", "CFD model part has to be set, first.")
        #     return

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "A model part has to be set, first.")
            return

        # if not self.CfdModelPart.HasSubModelPart("NoSlip"):
        if not self.ModelPart.HasSubModelPart("NoSlip"):
            # self.CfdModelPart.CreateSubModelPart("NoSlip")
            self.ModelPart.CreateSubModelPart("NoSlip")
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"NoSlip\" created right now!")
        
        # KratosGeo.FillCfdModelpartUtilities(self.CfdModelPart).FillNoslip(self.ModelPart, cond_sub_model_name)
        KratosGeo.FillCfdModelpartUtilities(self.ModelPart).FillNoslip(cond_sub_model_name)


    # UGLY!!!
    def SetName(self, file_name):
        "update file name in json file"

        self.cfd_param["problem_data"]["problem_name"].SetString(file_name)
        self.cfd_param["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(file_name)
        # self.cfd_param["solver_settings"]["model_import_settings"]["input_filename"].SetString(file_name)


    def NoSlip(self, sub_model_name):
        self.cfd_param["solver_settings"]["skin_parts"].Append(sub_model_name)
        
        no_slip_param = KratosMultiphysics.Parameters("""{
            "python_module" : "apply_noslip_process",
            "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
            "process_name"  : "ApplyNoSlipProcess",
            "Parameters"    : {
                "model_part_name" : "FluidModelPart."""+ sub_model_name +""""
            }
        }""")
        self.cfd_param["processes"]["boundary_conditions_process_list"].Append(no_slip_param)


    def Slip (self, sub_model_name):
        self.cfd_param["solver_settings"]["skin_parts"].Append(sub_model_name)

        slip_param = KratosMultiphysics.Parameters("""{
            "python_module" : "apply_slip_process",
            "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
            "process_name"  : "ApplySlipProcess",
            "Parameters"    : {
                "model_part_name" : "FluidModelPart."""+ sub_model_name +""""
            }
        }""")
        self.cfd_param["processes"]["boundary_conditions_process_list"].Append(slip_param)

    
    def Inlet_Outlet(self, n_sectors=12, dir_in=1):
        """ function to set the boundary conditions "inlet" and "outlet"

        Note:
             with 12 directions:
                1   = East
                4   = North
                7   = West
                10  = South
            

        Args:
            n_sectors: number of all sectors
            dir_in: incoming wind direction

        Returns:
            update the ModelPart with the boundary conditions
        """

        import math

        if (dir_in > n_sectors or dir_in <= 0):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Unable to analyze the current direction. Out of range!")
            return

        sect_in = dir_in - int(n_sectors/4)
        
        sect_out = dir_in + int(n_sectors/2)
        sect_out -= int(n_sectors/4)
        if sect_out > n_sectors:
            sect_out -= n_sectors

        print("inlet direzione: ", dir_in)
        # fill Inlet BC
        for n_sect in range(sect_in, sect_in + int(n_sectors/2)):
            if (n_sect <= 0):
                n_sect += n_sectors
            elif (n_sect > n_sectors):
                n_sect -= n_sectors
            print("\tLateralSector_{}".format(n_sect))
            sub_model_name = "LateralSector_{}".format(n_sect)
            
            theta = (360.0 / n_sectors) * (dir_in - 1)  # angle with the X axis
            dir_x = math.cos(math.radians(theta))
            dir_y = math.sin(math.radians(theta))
            dir_z = 0
            direction = [-dir_x, -dir_y, dir_z] # the minus to consider the incoming direction of the wind

            self.Inlet(sub_model_name, direction)
        print()

        print("outlet direzione: ", dir_in + int(n_sectors/2))
        # fill Outlet BC
        for n_sect in range(sect_out, sect_out + int(n_sectors/2)):
            if (n_sect <= 0):
                n_sect += n_sectors
            elif (n_sect > n_sectors):
                n_sect -= n_sectors
            print("\tLateralSector_{}".format(n_sect))
            sub_model_name = "LateralSector_{}".format(n_sect)
            self.Outlet(sub_model_name)
        print()


    def Inlet(self, sub_model_name, direction="automatic_inwards_normal"):
        self.cfd_param["solver_settings"]["skin_parts"].Append(sub_model_name)

        # TODO:
        #   with "modulus": 6m/s          (change with logarithmic profile)
        #   with "interval": [0,"End"]    (check if it is correct)
        inlet_param = KratosMultiphysics.Parameters("""{
            "python_module" : "apply_inlet_process",
            "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
            "process_name"  : "ApplyInletProcess",
            "Parameters"    : {
                "model_part_name" : "FluidModelPart."""+ sub_model_name +"""",
                "variable_name"   : "VELOCITY",
                "interval"        : [0,"End"],
                "modulus"         : "6",
                "direction"       : """ + str(direction) + """
            }
        }""")
        self.cfd_param["processes"]["boundary_conditions_process_list"].Append(inlet_param)


    def Outlet(self, sub_model_name):
        self.cfd_param["solver_settings"]["skin_parts"].Append(sub_model_name)

        output_param = KratosMultiphysics.Parameters("""{
            "python_module" : "apply_outlet_process",
            "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
            "process_name"  : "ApplyOutletProcess",
            "Parameters"    : {
                "model_part_name"    : "FluidModelPart."""+ sub_model_name +"""",
                "variable_name"      : "PRESSURE",
                "constrained"        : true,
                "value"              : 0.0,
                "hydrostatic_outlet" : false,
                "h_top"              : 0.0
            }
        }""")
        self.cfd_param["processes"]["boundary_conditions_process_list"].Append(output_param)


    def PartsFluid(self, sub_model_name):
        self.cfd_param["solver_settings"]["volume_model_part_name"].SetString(sub_model_name)
        self.cfd_param["processes"]["gravity"][0]["Parameters"]["model_part_name"].SetString("FluidModelPart." + sub_model_name)


### --- auxiliary functions --- ### -------------------------------------

    def _parameter_initialization(self):
        "Parameters initialization with default values"

        self.cfd_param = KratosMultiphysics.Parameters("""{
            "problem_data"     : {
                "problem_name"  : "",
                "parallel_type" : "OpenMP",
                "echo_level"    : 0,
                "start_time"    : 0.0,
                "end_time"      : 1.0
            },
            "output_processes" : {
                "gid_output" : [{
                    "python_module" : "gid_output_process",
                    "kratos_module" : "KratosMultiphysics",
                    "process_name"  : "GiDOutputProcess",
                    "help"          : "This process writes postprocessing files for GiD",
                    "Parameters"    : {
                        "model_part_name"        : "FluidModelPart.fluid_computational_model_part",
                        "output_name"            : "",
                        "postprocess_parameters" : {
                            "result_file_configuration" : {
                                "gidpost_flags"               : {
                                    "GiDPostMode"           : "GiD_PostBinary",
                                    "WriteDeformedMeshFlag" : "WriteDeformed",
                                    "WriteConditionsFlag"   : "WriteConditions",
                                    "MultiFileFlag"         : "SingleFile"
                                },
                                "file_label"                  : "time",
                                "output_control_type"         : "step",
                                "output_interval"             : 1,
                                "body_output"                 : true,
                                "node_output"                 : false,
                                "skin_output"                 : false,
                                "plane_output"                : [],
                                "nodal_results"               : ["VELOCITY","PRESSURE"],
                                "gauss_point_results"         : [],
                                "nodal_nonhistorical_results" : []
                            },
                            "point_data_configuration"  : []
                        }
                    }
                }]
            },
            "solver_settings"  : {
                "model_part_name"             : "FluidModelPart",
                "domain_size"                 : 3,
                "solver_type"                 : "Monolithic",
                "model_import_settings"       : {
                    "input_type"     : "mdpa",
                    "input_filename" : ""
                },
                "echo_level"                  : 0,
                "compute_reactions"           : false,
                "maximum_iterations"          : 10,
                "relative_velocity_tolerance" : 0.001,
                "absolute_velocity_tolerance" : 1e-5,
                "relative_pressure_tolerance" : 0.001,
                "absolute_pressure_tolerance" : 1e-5,
                "volume_model_part_name"      : "",
                "skin_parts"                  : [],
                "no_skin_parts"               : [],
                "time_scheme"                 : "bossak",
                "time_stepping"               : {
                    "automatic_time_step" : false,
                    "time_step"           : 0.1
                },
                "formulation"                 : {
                    "element_type"             : "vms",
                    "use_orthogonal_subscales" : false,
                    "dynamic_tau"              : 1.0
                },
                "reform_dofs_at_each_step"    : false
            },
            "processes"        : {
                "initial_conditions_process_list"  : [],
                "boundary_conditions_process_list" : [],
                "gravity"                          : [{
                    "python_module" : "assign_vector_by_direction_process",
                    "kratos_module" : "KratosMultiphysics",
                    "process_name"  : "AssignVectorByDirectionProcess",
                    "Parameters"    : {
                        "model_part_name" : "",
                        "variable_name"   : "BODY_FORCE",
                        "modulus"         : 0.0,
                        "constrained"     : false,
                        "direction"       : [0.0,-1.0,0.0]
                    }
                }],
                "auxiliar_process_list"            : []
            }
        }""")
