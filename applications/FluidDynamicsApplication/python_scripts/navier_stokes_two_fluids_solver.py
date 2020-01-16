from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#Import basic python Libraries
import numpy as np

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtilities

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
have_conv_diff = KratosUtilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication")
if have_conv_diff:
    import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver
from KratosMultiphysics.FluidDynamicsApplication.read_distance_from_file import DistanceImportUtility

import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

def CreateSolver(model, custom_settings):
    return NavierStokesTwoFluidsSolver(model, custom_settings)

class NavierStokesTwoFluidsSolver(FluidSolver):

    @classmethod
    def GetDefaultSettings(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "two_fluids_solver_from_defaults",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "distance_reading_settings"    : {
                "import_mode"         : "from_mdpa",
                "distance_file_name"  : "no_distance_file"
            },
            "maximum_iterations": 7,
            "echo_level": 0,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-2,
                "maximum_delta_time"  : 1.0,
                "time_step"           : 0.0
            },
            "periodic": "periodic",
            "move_mesh_flag": false,
            "formulation": {
                "dynamic_tau": 1.0
            },
            "bfecc_convection" : false,
            "bfecc_number_substeps" : 10
        }""")

        default_settings.AddMissingParameters(super(NavierStokesTwoFluidsSolver, cls).GetDefaultSettings())
        return default_settings

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(NavierStokesTwoFluidsSolver,self).__init__(model,custom_settings)

        self.element_name = "TwoFluidNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.element_has_nodal_properties = True

        self.min_buffer_size = 3

        self._bfecc_convection = self.settings["bfecc_convection"].GetBool()

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidsSolver", "Construction of NavierStokesTwoFluidsSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)              # Distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.DISTANCE_AUX)                   # Auxiliary distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)     # Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.DISTANCE_GRADIENT_AUX)     # Auxiliary Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CURVATURE)                      # Store curvature as a nodal variable
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.AREA_VARIABLE_AUX)              # Auxiliary area_variable for parallel distance calculator    
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.NORMAL_VECTOR)                  # Auxiliary normal vector at interface
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.TANGENT_VECTOR)                 # Auxiliary tangent vector at contact line
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CONTACT_VECTOR)                 # Auxiliary contact vector     

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidsSolver", "Fluid solver variables added correctly.")

    def PrepareModelPart(self):
        # Initialize the level-set function
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Setting the nodal distance
            self._set_distance_function()

        # Call the base solver PrepareModelPart()
        super(NavierStokesTwoFluidsSolver, self).PrepareModelPart()

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        ## Construct the linear solver
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.computing_model_part, self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.computing_model_part, 10, 10)
        (self.neighbour_search).Execute()

        self.accelerationLimitationUtility = KratosCFD.AccelerationLimitationUtilities( self.computing_model_part, 5.0 )

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Set the time discretization utility to compute the BDF coefficients
        time_order = self.settings["time_order"].GetInt()
        if time_order == 2:
            self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
        else:
            raise Exception("Only \"time_order\" equal to 2 is supported. Provided \"time_order\": " + str(time_order))

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.level_set_convection_process = self._set_level_set_convection_process()

        #self.mass_conservation_correction = self._set_mass_conservation_correction()
        #(self.mass_conservation_correction).Initialize();

        self.parallel_distance_process = self._set_parallel_distance_process()
        #(self.parallel_distance_process).CalculateDistances(
        #            self.main_model_part, 
        #            KratosMultiphysics.DISTANCE, 
        #            KratosCFD.AREA_VARIABLE_AUX, 
        #            2000, 
        #            0.01,
        #            (self.parallel_distance_process).CALCULATE_EXACT_DISTANCES_TO_PLANE)

        self.variational_distance_process = self._set_variational_distance_process()
        #(self.variational_distance_process).Execute()

        #self.lumped_eikonal_distance_calculation = self._set_lumped_eikonal_distance_calculation()
        #(self.lumped_eikonal_distance_calculation).Execute()

        self.surface_smoothing_process = self._set_surface_smoothing_process()
        #(self.surface_smoothing_process).Execute()
        
        #for node in self.main_model_part.Nodes:
        #    smooth_distance = node.GetSolutionStepValue(KratosCFD.DISTANCE_AUX)
        #    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, smooth_distance)

        self.distance_gradient_process = self._set_distance_gradient_process()
        #(self.distance_gradient_process).Execute()

        self.curvature_calculation_process = self._set_curvature_calculation_process()
        #(self.curvature_calculation_process).Execute()

        self.variational_non_eikonal_distance = self._set_variational_non_eikonal_distance()
        #(self.distance_gradient_process).Execute()
        #(self.curvature_calculation_process).Execute()
        #(self.variational_non_eikonal_distance).Execute()
        #for node in self.main_model_part.Nodes:
        #    smooth_distance = node.GetSolutionStepValue(KratosCFD.DISTANCE_AUX)
        #    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, smooth_distance)

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],   # Domain size (2,3)
                                                                                        self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]+1) # DOFs (3,4)

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        (self.solver).Initialize() # Initialize the solver. Otherwise the constitutive law is not initializated.
        (self.solver).Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["formulation"]["dynamic_tau"].GetDouble())

        #for node in (self.main_model_part.GetSubModelPart("NoSlip3D_No_Slip_Auto1")).Nodes:
        #    node.SetValue(KratosMultiphysics.IS_STRUCTURE, 1.0)
            #NodeId = node.Id
            #KratosMultiphysics.Logger.PrintInfo("Wall", NodeId)
        #Set IS_STRUCTURE to define contact line.

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidsSolver", "Solver initialization finished.")

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Recompute the BDF2 coefficients
            (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

            TimeStep = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

            # Correct the distance function according to volume conservation
            # Check the JASON properties file for duplicate processes!
            #if (TimeStep % 1 == 0):
            #    KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidsSolver", "About to impose mass conservation!")
            #    (self.mass_conservation_correction).ExecuteInTimeStep();

            # Perform the level-set convection according to the previous step velocity
            if self._bfecc_convection:
                KratosMultiphysics.Logger.PrintInfo("LevelSetSolver", "BFECCconvect will be called")
                (self.level_set_convection_process).BFECCconvect(
                    self.main_model_part,
                    KratosMultiphysics.DISTANCE,
                    KratosMultiphysics.VELOCITY,
                    self.settings["bfecc_number_substeps"].GetInt())
            else:
                (self.level_set_convection_process).Execute()

            # Recompute the distance field according to the new level-set position
            #if (TimeStep % 1 == 0):
            #    (self.variational_distance_process).Execute()

            # Recompute the distance field according to the new level-set position
            #if (TimeStep % 20 == 0):
            #    (self.parallel_distance_process).CalculateInterfacePreservingDistances( #CalculateDistances(
            #        self.main_model_part, 
            #        KratosMultiphysics.DISTANCE, 
            #        KratosCFD.AREA_VARIABLE_AUX, 
            #        500, 
            #        0.3)#,
                    #(self.parallel_distance_process).CALCULATE_EXACT_DISTANCES_TO_PLANE)

            # Reinitialize distance according to time dependent Eikonal equation
            #if (TimeStep % 1 == 0):
            #    (self.lumped_eikonal_distance_calculation).Execute()

            # Reinitialize distance using a new variational process needs Compute the DISTANCE_GRADIENT and CURVATURE on nodes
            if (TimeStep % 1 == 0):
                #(self.distance_gradient_process).Execute()
                #(self.curvature_calculation_process).Execute()
                (self.variational_non_eikonal_distance).Execute()
                for node in self.main_model_part.Nodes:
                    smooth_distance = node.GetSolutionStepValue(KratosCFD.DISTANCE_AUX)
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, smooth_distance)

            # Smoothing the surface to filter oscillatory surface
            #(self.surface_smoothing_process).Execute()

            #for node in self.main_model_part.Nodes:
            #    smooth_distance = node.GetSolutionStepValue(KratosCFD.DISTANCE_AUX)
            #    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, smooth_distance)

            # Compute the DISTANCE_GRADIENT on nodes
            (self.distance_gradient_process).Execute()

            # Compute CURVATURE on nodes
            (self.curvature_calculation_process).Execute()

            #for node in self.main_model_part.Nodes:
            #    node.SetSolutionStepValue(KratosCFD.CURVATURE, 1000.0)

            # Update the DENSITY and DYNAMIC_VISCOSITY values according to the new level-set
            self._SetNodalProperties()

            # Initialize the solver current step
            (self.solver).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            (self.solver).FinalizeSolutionStep()
            (self.accelerationLimitationUtility).Execute()

    # TODO: Remove this method as soon as the subproperties are available
    def _SetPhysicalProperties(self):
        import os
        warn_msg  = '\nThe materials import mechanism used in the two fluids solver is DEPRECATED!\n'
        warn_msg += 'It will be removed to use the base fluid_solver.py one as soon as the subproperties are available.\n'
        KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)

        # Check if the fluid properties are provided using a .json file
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            with open(materials_filename,'r') as materials_file:
                materials = KratosMultiphysics.Parameters(materials_file.read())

            # Create and read an auxiliary materials file for each one of the fields
            for i_material in materials["properties"]:
                aux_materials = KratosMultiphysics.Parameters()
                aux_materials.AddEmptyArray("properties")
                aux_materials["properties"].Append(i_material)
                prop_id = i_material["properties_id"].GetInt()

                aux_materials_filename = materials_filename + "_" + str(prop_id) + ".json"
                with open(aux_materials_filename,'w') as aux_materials_file:
                    aux_materials_file.write(aux_materials.WriteJsonString())
                    aux_materials_file.close()

                aux_material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
                aux_material_settings["Parameters"]["materials_filename"].SetString(aux_materials_filename)
                KratosMultiphysics.ReadMaterialsUtility(aux_material_settings, self.model)

                os.remove(aux_materials_filename)

            materials_imported = True
        else:
            materials_imported = False

        # If the element uses nodal material properties, transfer them to the nodes
        if self.element_has_nodal_properties:
            self._SetNodalProperties()

        return materials_imported

    def _SetNodalProperties(self):
        # Get fluid 1 and 2 properties
        properties_1 = self.main_model_part.Properties[1]
        properties_2 = self.main_model_part.Properties[2]

        rho_1 = properties_1.GetValue(KratosMultiphysics.DENSITY)
        rho_2 = properties_2.GetValue(KratosMultiphysics.DENSITY)
        mu_1 = properties_1.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        mu_2 = properties_2.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)

        # Check fluid 1 and 2 properties
        if rho_1 <= 0.0:
            raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho_1, properties_1.Id))
        if rho_2 <= 0.0:
            raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho_2, properties_2.Id))
        if mu_1 <= 0.0:
            raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(mu_1, properties_1.Id))
        if mu_2 <= 0.0:
            raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(mu_2, properties_2.Id))

        # Transfer density and (dynamic) viscostity to the nodes
        for node in self.main_model_part.Nodes:
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) <= 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho_1)
                node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, mu_1)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho_2)
                node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, mu_2)

    def _set_distance_function(self):
        ## Set the nodal distance function
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            DistanceUtility = DistanceImportUtility(self.main_model_part, self.settings["distance_reading_settings"])
            DistanceUtility.ImportDistance()
        elif (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_mdpa"):
            KratosMultiphysics.Logger.PrintInfo("Navier Stokes Embedded Solver","Distance function taken from the .mdpa input file.")

    def _set_level_set_convection_process(self):
        # Construct the level set convection process
        if self._bfecc_convection:
            if have_conv_diff:
                if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                    locator = KratosMultiphysics.BinBasedFastPointLocator2D(self.main_model_part)
                    locator.UpdateSearchDatabase()
                    level_set_convection_process = KratosConvDiff.BFECCConvection2D(locator)
                else:
                    locator = KratosMultiphysics.BinBasedFastPointLocator3D(self.main_model_part)
                    locator.UpdateSearchDatabase()
                    level_set_convection_process = KratosConvDiff.BFECCConvection3D(locator)

            else:
                raise Exception("The BFECC level set convection requires the Kratos ConvectionDiffusionApplication compilation.")
        else:
            if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess2D(
                    KratosMultiphysics.DISTANCE,
                    self.main_model_part,
                    self.linear_solver)
            else:
                level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess3D(
                    KratosMultiphysics.DISTANCE,
                    self.main_model_part,
                    self.linear_solver)

        return level_set_convection_process

    def _set_variational_distance_process(self):
        # Construct the variational distance calculation process
        maximum_iterations = 2 #TODO: Make this user-definable
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                self.main_model_part,
                self.linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE)
        else:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                self.main_model_part,
                self.linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE)

        return variational_distance_process

    def _set_parallel_distance_process(self):
        # Construct the parallel distance process
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                locator = KratosMultiphysics.BinBasedFastPointLocator2D(self.main_model_part)
                locator.UpdateSearchDatabase()
                parallel_distance_process = KratosMultiphysics.ParallelDistanceCalculator2D()
        else:
                locator = KratosMultiphysics.BinBasedFastPointLocator3D(self.main_model_part)
                locator.UpdateSearchDatabase()
                parallel_distance_process = KratosMultiphysics.ParallelDistanceCalculator3D()
        
        return parallel_distance_process

    #def _set_lumped_eikonal_distance_calculation(self):
        # Construct the process for redistancing
    #    lumped_eikonal_distance_calculation = KratosCFD.LumpedEikonalDistanceCalculation(
    #            self.main_model_part,
    #            10,
    #            1.0e-8,
    #            0.1)
    #    
    #    return lumped_eikonal_distance_calculation

    def _set_variational_non_eikonal_distance(self):
        #Distance re-initialization by solving a diffusion problem
        variational_non_eikonal_distance = KratosCFD.VariationalNonEikonalDistance(
                self.main_model_part, 
                self.linear_solver)
    
        return variational_non_eikonal_distance

    def _set_distance_gradient_process(self):
        #Calculate DISTANCE_GRADIENT at nodes using ComputeNodalGradientProcess
        distance_gradient_process = KratosMultiphysics.ComputeNodalGradientProcess(
                self.main_model_part, 
                KratosMultiphysics.DISTANCE,  #KratosCFD.DISTANCE_AUX,
                KratosMultiphysics.DISTANCE_GRADIENT, #KratosCFD.DISTANCE_GRADIENT_AUX,
                KratosMultiphysics.NODAL_AREA)

        return distance_gradient_process

    def _set_curvature_calculation_process(self):
        #Calculate curvature as divergence of normalized DISTANCE_GRADIENT at nodes using ComputeNodalNormalDivergenceProcess
        curvature_calculation_process = KratosMultiphysics.ComputeNodalNormalDivergenceProcess(
                self.main_model_part, 
                KratosMultiphysics.DISTANCE_GRADIENT, #KratosCFD.DISTANCE_GRADIENT_AUX,
                KratosCFD.CURVATURE,
                KratosMultiphysics.NODAL_AREA)

        return curvature_calculation_process

    def _set_surface_smoothing_process(self):
        #Smoothing the surface (zero DISTANCE) by solving a diffusion problem
        surface_smoothing_process = KratosCFD.SurfaceSmoothingProcess(
                self.main_model_part, 
                self.linear_solver)

        return surface_smoothing_process

    def _set_mass_conservation_correction(self):
        #Restoring the distance function to compensate for the mass loss, proportional to curvature 
        #(could be velocity instead) but this is associated with the way diffusive surface smoothing process works!
        surface_smoothing_process = KratosCFD.MassConservationCorrection(
                self.main_model_part, 
                True,
                "mass_conservation.log")

        return surface_smoothing_process

    
