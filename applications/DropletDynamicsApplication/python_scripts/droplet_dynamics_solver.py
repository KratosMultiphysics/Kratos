# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import auxiliary_solver_utilities
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.DropletDynamicsApplication as KratosDroplet
import numpy as np

# Import base class file
#from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver
#from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_two_fluids_solver import NavierStokesTwoFluidsSolver

from pathlib import Path

def CreateSolver(model, custom_settings):
    return DropletDynamicsSolver(model, custom_settings)


class DropletDynamicsSolver(PythonSolver):  # Before, it was derived from NavierStokesTwoFluidsSolver

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "two_fluids",
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
            "maximum_iterations": 7,
            "echo_level": 0,
            "time_order": 2,
            "time_scheme": "bdf2",
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": false,
            "consider_periodic_conditions": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "assign_neighbour_elements_to_conditions": true,
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
            "acceleration_limitation": true,
            "formulation": {
                "dynamic_tau": 1.0
            },
            "levelset_convection_settings": {
                "max_CFL" : 1.0,
                "max_substeps" : 0,
                "eulerian_error_compensation" : false,
                "element_type" : "levelset_convection_supg",
                "element_settings" : {
                    "dynamic_tau" : 0.0,
                    "cross_wind_stabilization_factor" : 0.7
                }
            },
            "contact_angle_settings": {
                "theta_advancing" : 130,
                "theta_receding" : 130
            },                                               
            "distance_reinitialization": "variational",
            "parallel_redistance_max_layers" : 25,
            "distance_smoothing": false,
            "distance_smoothing_coefficient": 1.0,
            "distance_modification_settings": {
                "model_part_name": "",
                "distance_threshold": 1e-5,
                "continuous_distance": true,
                "check_at_each_time_step": true,
                "avoid_almost_empty_elements": false,
                "deactivate_full_negative_elements": false
            }
        }""")

        default_settings.AddMissingParameters(super(DropletDynamicsSolver, cls).GetDefaultParameters())
        return default_settings


    def __init__(self, model, custom_settings):
        """Initializing the solver."""
        # TODO: DO SOMETHING IN HERE TO REMOVE THE "time_order" FROM THE DEFAULT SETTINGS BUT KEEPING THE BACKWARDS COMPATIBILITY

        if custom_settings.Has("levelset_convection_settings"):
            if custom_settings["levelset_convection_settings"].Has("levelset_splitting"):
                custom_settings["levelset_convection_settings"].RemoveValue("levelset_splitting")
                KratosMultiphysics.Logger.PrintWarning("NavierStokesTwoFluidsSolver", "\'levelset_splitting\' has been temporarily deactivated. Using the standard levelset convection with no splitting.")

        #TODO: Remove this after the retrocompatibility period
        if custom_settings.Has("bfecc_convection"):
            KratosMultiphysics.Logger.PrintWarning("NavierStokesTwoFluidsSolver", "the semi-Lagrangian \'bfecc_convection\' is no longer supported. Using the standard Eulerian levelset convection.")
            custom_settings.RemoveValue("bfecc_convection")
            if custom_settings.Has("bfecc_number_substeps"):
                custom_settings.RemoveValue("bfecc_number_substeps")

        # FluidSolver.__init__(self,model,custom_settings)
        super(DropletDynamicsSolver,self).__init__(model, custom_settings)

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        self.element_name = "DropletDynamics"
        self.condition_name = "TwoFluidNavierStokesWallCondition"
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = True

        self.min_buffer_size = 3

        # Set the levelset characteristic variables and add them to the convection settings
        # These are required to be set as some of the auxiliary processes admit user-defined variables
        self._levelset_variable = KratosMultiphysics.DISTANCE
        self._levelset_gradient_variable = KratosMultiphysics.DISTANCE_GRADIENT
        self._levelset_convection_variable = KratosMultiphysics.VELOCITY
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_variable_name").SetString("DISTANCE")
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_gradient_variable_name").SetString("DISTANCE_GRADIENT")
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_convection_variable_name").SetString("VELOCITY")

        dynamic_tau = self.settings["formulation"]["dynamic_tau"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, dynamic_tau)

        surface_tension = True
        #if (self.settings["formulation"].Has("surface_tension")):
        #    surface_tension = self.settings["formulation"]["surface_tension"].GetBool()
        self.main_model_part.ProcessInfo.SetValue(KratosCFD.SURFACE_TENSION, surface_tension)

        self.momentum_correction = True
        #if self.settings["formulation"].Has("momentum_correction"):
        #    self.momentum_correction = self.settings["formulation"]["momentum_correction"].GetBool()
        self.main_model_part.ProcessInfo.SetValue(KratosCFD.MOMENTUM_CORRECTION, self.momentum_correction)

        self._reinitialization_type = self.settings["distance_reinitialization"].GetString()

        self._distance_smoothing = self.settings["distance_smoothing"].GetBool()
        smoothing_coefficient = self.settings["distance_smoothing_coefficient"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosCFD.SMOOTHING_COEFFICIENT, smoothing_coefficient)

        self._apply_acceleration_limitation = self.settings["acceleration_limitation"].GetBool()

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        #if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
        #    self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesTwoFluidsSolver finished.")


    def AddDofs(self):
        dofs_and_reactions_to_add = []
        dofs_and_reactions_to_add.append(["VELOCITY_X", "REACTION_X"])
        dofs_and_reactions_to_add.append(["VELOCITY_Y", "REACTION_Y"])
        dofs_and_reactions_to_add.append(["VELOCITY_Z", "REACTION_Z"])
        dofs_and_reactions_to_add.append(["PRESSURE", "REACTION_WATER_PRESSURE"])
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver DOFs added correctly.")


    def GetDofsList(self):
        """This function creates and returns a list with the DOFs defined in the conditions and elements specifications
        Note that this requires the main_model_part to be already set, that is to say to have already performed the element substitution (see PrepareModelPart).
        """
        return KratosMultiphysics.SpecificationsUtilities.GetDofsListFromSpecifications(self.main_model_part)


    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])


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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)     # Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.AUX_DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.DISTANCE_AUX)                   # Auxiliary distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.DISTANCE_AUX2)                  # Auxiliary distance function nodal values       
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.DISTANCE_GRADIENT_AUX)          # Auxiliary Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.CONVECTIVE_VELOCITY)            # Store conctive velocity for level-set process
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CURVATURE)                      # Store curvature as a nodal variable
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.AREA_VARIABLE_AUX)              # Auxiliary area_variable for parallel distance calculator
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.NORMAL_VECTOR)                  # Auxiliary normal vector at interface
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.TANGENT_VECTOR)                 # Auxiliary tangent vector at contact line
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.CONTACT_VECTOR)                 # Auxiliary contact vector
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.CONTACT_ANGLE)                  # Contact angle (may not be needed at nodes)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.CONTACT_VECTOR_MICRO)           # Auxiliary contact vector at micro-scale
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.CONTACT_ANGLE_MICRO)            # Contact angle (micro-scale)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.CONTACT_VELOCITY)               # Contact line tangential velocity (normal to the contact-line)
        #self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.VELOCITY_STAR)                  # Last known velocity
        #self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.PRESSURE_STAR)                  # Last known pressure
        # self.main_model_part.AddNodalSolutionStepVariable(KratosDroplet.PRESSURE_GRADIENT_AUX)          # Pressure gradient on positive and negative sides
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_VELOCITY)

    def PrepareModelPart(self):
        # Restarting the simulation is OFF (needs a careful implementation)
        # Set fluid properties from materials json file
        materials_imported = self._SetPhysicalProperties()
        if not materials_imported:
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, "Material properties have not been imported. Check \'material_import_settings\' in your ProjectParameters.json.")
        # Replace default elements and conditions
        self._ReplaceElementsAndConditions()
        # Set and fill buffer
        self._SetAndFillBuffer()

        # Executes the check and prepare model process. Always executed as it also assigns neighbors which are not saved in a restart
        self._ExecuteCheckAndPrepare()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Model reading finished.")


    def ExportModelPart(self):
        # Writing the model part
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Model export finished.")

   
    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    
    # Temporary name: "fluid_computational_model_part"
    def GetComputingModelPart(self):
        if not self.main_model_part.HasSubModelPart("fluid_computational_model_part"):
            raise Exception("The ComputingModelPart was not created yet!")
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")


    def Initialize(self):
        computing_model_part = self.GetComputingModelPart()
        # Calculate boundary normals
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(
            computing_model_part,
            computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Finding nodal and elemental neighbors
        data_communicator = computing_model_part.GetCommunicator().GetDataCommunicator()
        neighbour_search = KratosMultiphysics.FindGlobalNodalNeighboursProcess(
            data_communicator,
            computing_model_part)
        neighbour_search.Execute()

        elemental_neighbour_search = KratosMultiphysics.GenericFindElementalNeighboursProcess(
            computing_model_part)
        elemental_neighbour_search.Execute()

        # Set and initialize the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        # Set nodal properties after setting distance(level-set).
        self._SetNodalProperties()

        # Initialize the distance correction process
        self._GetDistanceModificationProcess().ExecuteInitialize()
        self._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        #Here the initial water volume of the system is calculated without considering inlet and outlet flow rate
        self.initial_system_volume=KratosCFD.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(self.GetComputingModelPart())

        # Instantiate the level set convection process
        # Note that is is required to do this in here in order to validate the defaults and set the corresponding distance gradient flag
        # Note that the nodal gradient of the distance is required either for the eulerian BFECC limiter or by the algebraic element antidiffusivity
        self._GetLevelSetConvectionProcess()

        # Just to be sure that no conflict would occur if the element is derived from the Navier Stokes Two Fluids.
        self.mass_source = False

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

        # The external interfacial force (per unit area) should be set in case of the presence of electromagnetic forces, etc.
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosDroplet.EXT_INT_FORCE, self.main_model_part.Elements)


    def AdvanceInTime(self, current_time):
        dt = self.settings["time_stepping"]["time_step"].GetDouble() # The automatic dt calculation is deactivated.
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    
    def InitializeSolutionStep(self):
        # Momentum correction is on by default!
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosCFD.DISTANCE_CORRECTION, 0.0, self.main_model_part.Nodes)

        # Recompute the BDF2 coefficients
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

        # Perform the level-set convection according to the previous step velocity
        self._PerformLevelSetConvection()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Level-set convection is performed.")

        # Clear any existing intersection points from previous steps
        KratosDroplet.IntersectionPointsUtility.ClearIntersectionPoints()
    
        # Collect intersection points from all elements
        for element in self.main_model_part.Elements:
            KratosDroplet.IntersectionPointsUtility.CollectElementIntersectionPoints(element)
        
        # # Run diagnostic to check how many elements are split by the level-set
        # KratosDroplet.IntersectionPointsUtility.DiagnosticOutput(self.main_model_part)
    
        # Get all intersection points
        points = KratosDroplet.IntersectionPointsUtility.GetIntersectionPoints()
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, f"Collected {len(points)} intersection points.")
    
        # Save intersection points to file
        KratosDroplet.IntersectionPointsUtility.SaveIntersectionPointsToFile("intersection_points.txt")
    
        # Fit curves to the points on an element-by-element basis
        KratosDroplet.IntersectionPointsUtility.ProcessIntersectionPointsAndFitCurves("element_curves.txt")
        KratosDroplet.IntersectionPointsUtility.ProcessIntersectionPointsAndFitCurvesparabola("element_curves_parabola.txt")
    

        # filtering noises is necessary for curvature calculation
        # distance gradient is used as a boundary condition for smoothing process
        self._GetDistanceGradientProcess().Execute()
        self._GetDistanceSmoothingProcess().Execute()
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Smoothing process is finished.")

        # distance gradient is called again to comply with the smoothed/modified DISTANCE
        self._GetDistanceGradientProcess().Execute()

        for node in self.main_model_part.Nodes:
            gx = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT_X)
            gy = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT_Y)
            gz = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT_Z)
            g = (gx**2+gy**2+gz**2)**0.5
            gx /= g
            gy /= g
            gz /= g
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT_X,gx)
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT_Y,gy)
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT_Z,gz)

        # self.ExtrapolateBoundaryValues(KratosMultiphysics.DISTANCE_GRADIENT_Y)
        # self.ExtrapolateBoundaryValues(KratosMultiphysics.DISTANCE_GRADIENT_X)

        # curvature is calculated using nodal distance gradient
        self._GetDistanceCurvatureProcess().Execute()

        ##########
        # Contact angle calculation
        # self._GetContactAngleEvaluatorProcess().Execute()
        # Store current level-set to check for wetting/dewetting used in contact_angle_evaluator
        for node in self.main_model_part.Nodes:
            old_distance = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            node.SetValue(KratosDroplet.DISTANCE_AUX, old_distance)
        # print("Contact Angle Evaluator: Finished")
        ##########

        # it is needed to store level-set consistent nodal PRESSURE_GRADIENT for stabilization purpose
        self._GetConsistentNodalPressureGradientProcess().Execute()

        # TODO: Performing mass conservation check and correction process

        # Perform distance correction to prevent ill-conditioned cuts
        self._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        # Update the DENSITY and DYNAMIC_VISCOSITY values according to the new level-set
        self._SetNodalProperties()

        # Initialize the solver current step
        self._GetSolutionStrategy().InitializeSolutionStep()

        # We set this value at every time step as other processes/solvers also use them
        dynamic_tau = self.settings["formulation"]["dynamic_tau"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, dynamic_tau)


    def Predict(self):
        self._GetSolutionStrategy().Predict()


    def SolveSolutionStep(self):
        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        if not is_converged:
            msg  = "Droplet dynamics solver did not converge for step " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
            msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, msg)
        return is_converged


    def FinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Mass and momentum conservation equations are solved.")

        # Recompute the distance field according to the new level-set position
        if self._reinitialization_type != "none":
            step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
            if step == 2 or step%12==0:
                self._GetDistanceReinitializationProcess().Execute()
                KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Redistancing process is finished.")

        # Prepare distance correction for next step
        self._GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()

        # Finalize the solver current step
        self._GetSolutionStrategy().FinalizeSolutionStep()
        # Limit the obtained acceleration for the next step
        # This limitation should be called on the second solution step onwards (e.g. STEP=3 for BDF2)
        # We intentionally avoid correcting the acceleration in the first resolution step as this might cause problems with zero initial conditions
        if self._apply_acceleration_limitation and self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.min_buffer_size:
            self._GetAccelerationLimitationUtility().Execute()


    def Check(self):
        self._GetSolutionStrategy().Check()


    def Clear(self):
        self._GetSolutionStrategy().Clear()


    def _SetPhysicalProperties(self):   # Might be removed from the NavierStokesTwoFluidsSolver!!! 
                                        # It is just a duplicate of the same routine defined in NavierStokesTwoFluidsSolver
        warn_msg  = '\nThe materials import mechanism used in the two fluids solver is DEPRECATED!\n'
        warn_msg += 'It will be removed to use the base fluid_solver.py one as soon as the subproperties are available.\n'
        KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)

        # Check if the fluid properties are provided using a .json file
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            data_comm = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator()

            def GetAuxMaterialsFileName(mat_file_name, prop_id):
                p_mat_file_name = Path(mat_file_name)
                new_stem = "{}_p{}".format(p_mat_file_name.stem, prop_id)
                return str(p_mat_file_name.with_name(new_stem).with_suffix(p_mat_file_name.suffix))

            with open(materials_filename,'r') as materials_file:
                materials = KratosMultiphysics.Parameters(materials_file.read())

            if data_comm.Rank() == 0:
                # Create and read an auxiliary materials file for each one of the fields (only on one rank)
                for i_material in materials["properties"]:
                    aux_materials = KratosMultiphysics.Parameters()
                    aux_materials.AddEmptyArray("properties")
                    aux_materials["properties"].Append(i_material)

                    aux_materials_filename = GetAuxMaterialsFileName(materials_filename, i_material["properties_id"].GetInt())
                    with open(aux_materials_filename,'w') as aux_materials_file:
                        aux_materials_file.write(aux_materials.WriteJsonString())

            data_comm.Barrier()

            # read the files on all ranks
            for i_material in materials["properties"]:
                aux_materials_filename = GetAuxMaterialsFileName(materials_filename, i_material["properties_id"].GetInt())
                aux_material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
                aux_material_settings["Parameters"]["materials_filename"].SetString(aux_materials_filename)
                KratosMultiphysics.ReadMaterialsUtility(aux_material_settings, self.model)

            data_comm.Barrier()

            materials_imported = True
        else:
            materials_imported = False

         # If the element uses nodal material properties, transfer them to the nodes
        if self.element_has_nodal_properties:
            self._SetNodalProperties()

        return materials_imported


    def _SetNodalProperties(self):
        # Keep it for now, this function might need more parameters for EHD
        # If the element uses nodal material properties, transfer them to the nodes
        if self.element_has_nodal_properties:
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
                if node.GetSolutionStepValue(self._levelset_variable) <= 0.0:
                    node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho_1)
                    node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, mu_1)
                else:
                    node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho_2)
                    node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, mu_2)

    
    # This routine is the duplicate of the same one defined in FluidSolver 
    def _ReplaceElementsAndConditions(self):
        ## Get number of nodes and domain size
        elem_num_nodes = self._GetElementNumNodes()
        cond_num_nodes = self._GetConditionNumNodes()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## If there are no elements and/or conditions, default to triangles/tetra meshes to avoid breaking the ReplaceElementsAndConditionsProcess
        ## This only affects the input name (if there are no elements or conditions to replace, nothing is replaced).
        if elem_num_nodes == 0:
            elem_num_nodes = domain_size + 1
        if cond_num_nodes == 0:
            cond_num_nodes = domain_size

        ## Complete the element name
        if (self.element_name is not None):
            new_elem_name = self.element_name + str(int(domain_size)) + "D" + str(int(elem_num_nodes)) + "N"
        else:
            raise Exception("There is no element name. Define the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if (self.condition_name is not None):
            new_cond_name = self.condition_name + str(int(domain_size)) + "D" + str(int(cond_num_nodes)) + "N"
        else:
            raise Exception("There is no condition name. Define the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        #self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""{}""")
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()


    # This routine is the duplicate of the same one defined in FluidSolver
    def _SetAndFillBuffer(self):
        init_dt = self.settings["time_stepping"]["time_step"].GetDouble() # The automatic dt calculation is deactivated.
        auxiliary_solver_utilities.SetAndFillBuffer(self.main_model_part, self.min_buffer_size, init_dt)


    # This routine is the duplicate of the same one defined in FluidSolver
    def _ExecuteCheckAndPrepare(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])
        prepare_model_part_settings.AddValue("assign_neighbour_elements_to_conditions",self.settings["assign_neighbour_elements_to_conditions"])

        # CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()
        if prepare_model_part_settings["volume_model_part_name"].GetString() == "":
            raise Exception("Please define the \"volume_model_part_name\" (string) argument.")

        volume_model_part_name = prepare_model_part_settings["volume_model_part_name"].GetString()
        skin_name_list = prepare_model_part_settings["skin_parts"]

        if self.main_model_part.Name == volume_model_part_name:
            self.volume_model_part = self.main_model_part
        else:
            self.volume_model_part = self.main_model_part.GetSubModelPart(volume_model_part_name)

        skin_parts = []
        for i in range(skin_name_list.size()):
            skin_parts.append(self.main_model_part.GetSubModelPart(skin_name_list[i].GetString()))

        # Temporary name: "fluid_computational_model_part"
        if self.main_model_part.HasSubModelPart("fluid_computational_model_part"):
            fluid_computational_model_part = self.main_model_part.GetSubModelPart("fluid_computational_model_part")
        else:
            fluid_computational_model_part = self.main_model_part.CreateSubModelPart("fluid_computational_model_part")
            fluid_computational_model_part.ProcessInfo = self.main_model_part.ProcessInfo

            for node in self.volume_model_part.Nodes:
                fluid_computational_model_part.AddNode(node,0)
            for elem in self.volume_model_part.Elements:
                fluid_computational_model_part.AddElement(elem,0)

            list_of_ids = set()
            for part in skin_parts:
                for cond in part.Conditions:
                    list_of_ids.add(cond.Id)

            fluid_computational_model_part.AddConditions(list(list_of_ids))

        # Orientation of the elements: only for trangles and tetrahedrons (simplex elements)
        geometry = self.main_model_part.Elements.__iter__().__next__().GetGeometry()
        is_simplex = geometry.LocalSpaceDimension() + 1 == geometry.PointsNumber()
        if not is_simplex:
            msg = "Geoemetry is not simplex. Orientation check is only available"
            msg += " for simplex geometries and hence it will be skipped."
            KratosMultiphysics.Logger.PrintWarning(type(self).__name__, msg)
            return 0

        tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse()
        # By default the neighboring elements are assigned
        flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS

        KratosMultiphysics.TetrahedralMeshOrientationCheck(fluid_computational_model_part,throw_errors, flags).Execute()

    
    def _GetElementNumNodes(self):
        if self.main_model_part.NumberOfElements() != 0:
            element_num_nodes = len(self.main_model_part.Elements.__iter__().__next__().GetNodes())
        else:
            element_num_nodes = 0

        element_num_nodes = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(element_num_nodes)
        return element_num_nodes


    def _GetConditionNumNodes(self):
        if self.main_model_part.NumberOfConditions() != 0:
            condition_num_nodes = len(self.main_model_part.Conditions.__iter__().__next__().GetNodes())
        else:
            condition_num_nodes = 0

        condition_num_nodes = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(condition_num_nodes)
        return condition_num_nodes
    
    
    def _GetScheme(self):
        if not hasattr(self, '_scheme'):
            self._scheme = self._CreateScheme()
        return self._scheme
    
    def _CreateScheme(self):
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # Here, the element incorporates the time integration scheme
        # It is required to perform the nodal update once the current time step is solved
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(
            domain_size,
            domain_size + 1)
        # Tthe BDF time discretization utility is required to update the BDF coefficients
        if (self.settings["time_scheme"].GetString() == "bdf2"):
            time_order = 2
            self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
        else:
            err_msg = "Requested time integration scheme \"" + self.settings["time_scheme"].GetString()+ "\" is not available.\n"
            raise Exception(err_msg)
        return scheme


    def _GetConvergenceCriterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._CreateConvergenceCriterion()
        return self._convergence_criterion

    def _CreateConvergenceCriterion(self):
        convergence_criterion = KratosMultiphysics.MixedGenericCriteria(
                [(KratosMultiphysics.VELOCITY, self.settings["relative_velocity_tolerance"].GetDouble(), self.settings["absolute_velocity_tolerance"].GetDouble()),
                (KratosMultiphysics.PRESSURE, self.settings["relative_pressure_tolerance"].GetDouble(), self.settings["absolute_pressure_tolerance"].GetDouble())])
        convergence_criterion.SetEchoLevel(self.settings["echo_level"].GetInt())
        return convergence_criterion


    def _GetLinearSolver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._CreateLinearSolver()
        return self._linear_solver

    def _CreateLinearSolver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        return linear_solver_factory.ConstructSolver(linear_solver_configuration)


    def _GetBuilderAndSolver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._CreateBuilderAndSolver()
        return self._builder_and_solver

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        if self.settings["consider_periodic_conditions"].GetBool():
            builder_and_solver = KratosCFD.ResidualBasedBlockBuilderAndSolverPeriodic(
                linear_solver,
                KratosCFD.PATCH_INDEX)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        return builder_and_solver


    def _GetSolutionStrategy(self):
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._CreateSolutionStrategy()
        return self._solution_strategy

    def _CreateSolutionStrategy(self):
        # Only the nonlinear (Newton-Raphson) strategy is available.
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
            computing_model_part,
            time_scheme,
            convergence_criterion,
            builder_and_solver,
            self.settings["maximum_iterations"].GetInt(),
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["move_mesh_flag"].GetBool())


    def _GetDistanceModificationProcess(self):
        if not hasattr(self, '_distance_modification_process'):
            self._distance_modification_process = self.__CreateDistanceModificationProcess()
        return self._distance_modification_process

    def __CreateDistanceModificationProcess(self):
        # Set suitable distance correction settings for free-surface problems
        # Note that the distance modification process is applied to the computing model part
        distance_modification_settings = self.settings["distance_modification_settings"]
        distance_modification_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["distance_modification_settings"])
        distance_modification_settings["model_part_name"].SetString(self.GetComputingModelPart().FullName())

        # Check user provided settings
        if not distance_modification_settings["continuous_distance"].GetBool():
            distance_modification_settings["continuous_distance"].SetBool(True)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'continuous_distance\' is \'False\'. Setting to \'True\'.")
        if not distance_modification_settings["check_at_each_time_step"].GetBool():
            distance_modification_settings["check_at_each_time_step"].SetBool(True)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'check_at_each_time_step\' is \'False\'. Setting to \'True\'.")
        if distance_modification_settings["avoid_almost_empty_elements"].GetBool():
            distance_modification_settings["avoid_almost_empty_elements"].SetBool(False)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'avoid_almost_empty_elements\' is \'True\'. Setting to \'False\' to avoid modifying the distance sign.")
        if distance_modification_settings["deactivate_full_negative_elements"].GetBool():
            distance_modification_settings["deactivate_full_negative_elements"].SetBool(False)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'deactivate_full_negative_elements\' is \'True\'. Setting to \'False\' to avoid deactivating the negative volume (e.g. water).")

        # Create and return the distance correction process
        return KratosCFD.DistanceModificationProcess(
            self.model,
            distance_modification_settings)


    def _GetAccelerationLimitationUtility(self):
        if not hasattr(self, '_acceleration_limitation_utility'):
            self._acceleration_limitation_utility = self.__CreateAccelerationLimitationUtility()
        return self._acceleration_limitation_utility
    
    def __CreateAccelerationLimitationUtility(self):
        maximum_multiple_of_g_acceleration_allowed = 5.0
        acceleration_limitation_utility = KratosCFD.AccelerationLimitationUtilities(
            self.GetComputingModelPart(),
            maximum_multiple_of_g_acceleration_allowed)

        return acceleration_limitation_utility


    def _PerformLevelSetConvection(self):
        # Solve the levelset convection problem
        self._GetLevelSetConvectionProcess().Execute()
    
    def _GetLevelSetConvectionProcess(self):
        if not hasattr(self, '_level_set_convection_process'):
            self._level_set_convection_process = self._CreateLevelSetConvectionProcess()
        return self._level_set_convection_process

    def _CreateLevelSetConvectionProcess(self):
        # Construct the level set convection process
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        computing_model_part = self.GetComputingModelPart()
        linear_solver = self._GetLevelsetLinearSolver()
        levelset_convection_settings = self.settings["levelset_convection_settings"]
        if domain_size == 2:
            level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess2D(
                computing_model_part,
                linear_solver,
                levelset_convection_settings)
        else:
            level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess3D(
                computing_model_part,
                linear_solver,
                levelset_convection_settings)

        return level_set_convection_process

    def _GetLevelsetLinearSolver(self):
        # A linear solver configured specifically for the level-set convection process
        if not hasattr(self, '_levelset_linear_solver'):
            self._levelset_linear_solver = self._CreateLinearSolver() # TODO: add customized configuration
        return self._levelset_linear_solver


    def _GetDistanceReinitializationProcess(self):
        if not hasattr(self, '_distance_reinitialization_process'):
            self._distance_reinitialization_process = self._CreateDistanceReinitializationProcess()
        return self._distance_reinitialization_process

    def _CreateDistanceReinitializationProcess(self):
        # Construct the variational distance calculation process
        if (self._reinitialization_type == "variational"):
            maximum_iterations = 2 #TODO: Make this user-definable
            linear_solver = self._GetRedistancingLinearSolver()
            computing_model_part = self.GetComputingModelPart()
            if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                distance_reinitialization_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                    computing_model_part,
                    linear_solver,
                    maximum_iterations,
                    KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE)
            else:
                distance_reinitialization_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                    computing_model_part,
                    linear_solver,
                    maximum_iterations,
                    KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE)

        elif (self._reinitialization_type == "parallel"):
            #TODO: move all this to solver settings
            layers = self.settings["parallel_redistance_max_layers"].GetInt()
            parallel_distance_settings = KratosMultiphysics.Parameters("""{
                "max_levels" : 25,
                "max_distance" : 1.0,
                "calculate_exact_distances_to_plane" : true
            }""")
            parallel_distance_settings["max_levels"].SetInt(layers)
            if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                distance_reinitialization_process = KratosMultiphysics.ParallelDistanceCalculationProcess2D(
                    self.main_model_part,
                    parallel_distance_settings)
            else:
                distance_reinitialization_process = KratosMultiphysics.ParallelDistanceCalculationProcess3D(
                    self.main_model_part,
                    parallel_distance_settings)
        elif (self._reinitialization_type == "none"):
                KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Redistancing is turned off.")
        else:
            raise Exception("Please use a valid distance reinitialization type or set it as \'none\'. Valid types are: \'variational\' and \'parallel\'.")

        return distance_reinitialization_process

    def _GetRedistancingLinearSolver(self):
        # A linear solver configured specifically for distance re-initialization process
        if not hasattr(self, '_redistancing_linear_solver'):
            self._redistancing_linear_solver = self._CreateLinearSolver() # TODO: add customized configuration
        return self._redistancing_linear_solver


    def _GetDistanceSmoothingProcess(self):
        if not hasattr(self, '_distance_smoothing_process'):
            self._distance_smoothing_process = self._CreateDistanceSmoothingProcess()
        return self._distance_smoothing_process

    def _CreateDistanceSmoothingProcess(self):
        # construct the distance smoothing process
        linear_solver = self._GetSmoothingLinearSolver()
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            distance_smoothing_process = KratosCFD.DistanceSmoothingProcess2D(
            self.main_model_part,
            linear_solver)
        else:
            distance_smoothing_process = KratosCFD.DistanceSmoothingProcess3D(
            self.main_model_part,
            linear_solver)

        return distance_smoothing_process

    def _GetSmoothingLinearSolver(self):
        # A linear solver configured specifically for the distance smoothing process
        if not hasattr(self, '_smoothing_linear_solver'):
            self._smoothing_linear_solver = self._CreateLinearSolver() # TODO: add customized configuration
        return self._smoothing_linear_solver


    def _GetDistanceGradientProcess(self):
        if not hasattr(self, '_distance_gradient_process'):
            self._distance_gradient_process = self._CreateDistanceGradientProcess()
        return self._distance_gradient_process

    def _CreateDistanceGradientProcess(self):
        distance_gradient_process = KratosMultiphysics.ComputeNodalGradientProcess(
                self.main_model_part,
                self._levelset_variable,
                self._levelset_gradient_variable,
                KratosMultiphysics.NODAL_AREA)

        return distance_gradient_process


    def _GetDistanceCurvatureProcess(self):
        if not hasattr(self, '_distance_curvature_process'):
            self._distance_curvature_process = self._CreateDistanceCurvatureProcess()
        return self._distance_curvature_process

    def _CreateDistanceCurvatureProcess(self):
        distance_curvature_process = KratosMultiphysics.ComputeNonHistoricalNodalNormalDivergenceProcess(
                self.main_model_part,
                self._levelset_gradient_variable,
                KratosCFD.CURVATURE,
                KratosMultiphysics.NODAL_AREA)

        return distance_curvature_process

    
    def _GetConsistentNodalPressureGradientProcess(self):
        if not hasattr(self, '_consistent_nodal_pressure_gradient_process'):
            self._consistent_nodal_pressure_gradient_process = self._CreateConsistentNodalPressureGradientProcess()
        return self._consistent_nodal_pressure_gradient_process

    def _CreateConsistentNodalPressureGradientProcess(self):
        consistent_nodal_pressure_gradient_process = KratosCFD.CalulateLevelsetConsistentNodalGradientProcess(
                self.main_model_part)

        return consistent_nodal_pressure_gradient_process
    
    def _GetContactAngleEvaluatorProcess(self):
        if not hasattr(self, '_distance_curvature_process'):
            self._distance_curvature_process = self._CreateContactAngleEvaluatorProcess()
        return self._distance_curvature_process

    def _CreateContactAngleEvaluatorProcess(self):
        contact_angle_settings = self.settings["contact_angle_settings"]
        contact_angle_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["contact_angle_settings"])
        contact_angle_evaluator = KratosDroplet.ContactAngleEvaluatorProcess(self.main_model_part, contact_angle_settings)

        return contact_angle_evaluator
    



    def FitPlaneAndExtrapolate(self, neighbor_nodes, variable):
        """
        Fit a linear plane (z = ax + by + c) using least squares and return the coefficients.

        :param neighbor_nodes: List of neighbor nodes.
        :param variable: The Kratos variable to be extrapolated.
        :return: (a, b, c) coefficients of the plane equation.
        """
        X, Y, Z = [], [], []

        for node in neighbor_nodes:
            X.append(node.X)
            Y.append(node.Y)
            Z.append(node.GetSolutionStepValue(variable))  # Get the value of the variable

        # Solve the least squares problem: [X Y 1] * [a b c] = Z
        A = np.vstack([X, Y, np.ones(len(X))]).T
        coeffs, _, _, _ = np.linalg.lstsq(A, Z, rcond=None)  # Solve for [a, b, c]

        return coeffs  # Returns (a, b, c)
    
    def ExtrapolateBoundaryValues(self, variable):
        """
        Extrapolate values for a given variable on boundary nodes using polynomial fitting.
    
        :param model_part: The Kratos ModelPart containing nodes and elements.
        :param variable: The Kratos variable to extrapolate (e.g., KM.DISTANCE, KM.TEMPERATURE).
        """
        for node in self.main_model_part.Nodes:
            if node.Is(KratosMultiphysics.BOUNDARY):  # Only process boundary nodes
                neighbor_nodes = set()

                # Find elements that contain this node
                for elem in self.main_model_part.Elements:
                    if node.Id in [n.Id for n in elem.GetNodes()]:  # Check if the node is in the element
                        for neighbor in elem.GetNodes():
                            if not neighbor.Is(KratosMultiphysics.BOUNDARY):  # Ensure we get non-boundary neighbors
                                neighbor_nodes.add(neighbor)

                # Ensure at least 3 neighbors for polynomial fitting
                if len(neighbor_nodes) >= 3:
                    a, b, c = self.FitPlaneAndExtrapolate(neighbor_nodes, variable)
                    estimated_value = a * node.X + b * node.Y + c
                    node.SetSolutionStepValue(variable, estimated_value)



