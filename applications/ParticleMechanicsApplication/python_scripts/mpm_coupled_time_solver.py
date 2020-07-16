from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

from KratosMultiphysics.ParticleMechanicsApplication.mpm_solver import MPMSolver

# Other imports
from KratosMultiphysics import auxiliary_solver_utilities

def CreateSolver(model, custom_settings):
    return MPMCoupledTimeSolver(model, custom_settings)

class MPMCoupledTimeSolver(MPMSolver):

    ### Solver constructor
    def __init__(self, model, custom_settings):
        print('\n\n\n ====================== USING COUPLED TIME SOLVER ====================== \n\n\n')
        #self._validate_settings_in_baseclass=True # To be removed eventually
        super(MPMCoupledTimeSolver, self).__init__(model, custom_settings)

        # Add model part containers
        self._AddModelPartContainers()

        # Default settings
        self.min_buffer_size = 2
        self.tolerance = 1e-6

        self.time_step_1 = self.settings["time_stepping"]["time_step"].GetDouble()
        self.time_step_ratio = self.settings["time_stepping"]["timestep_ratio"].GetInt()
        self.time_step_2 = self.time_step_1/self.time_step_ratio
        self.time_step_index_j = 0;

        self.gamma_1 = 10.0
        self.gamma_2 = 10.0

        if (self.settings["time_stepping"]["time_integration_method_1"].GetString() == "implicit"):
            self.gamma_1 = 0.5
        elif (self.settings["time_stepping"]["time_integration_method_1"].GetString() == "explicit"):
            self.gamma_1 = 1.0
        else:
            raise Exception("Check the time_integration_method_1 in project parameters file")
        if (self.settings["time_stepping"]["time_integration_method_2"].GetString() == "implicit"):
            self.gamma_2 = 0.5
        elif (self.settings["time_stepping"]["time_integration_method_2"].GetString() == "explicit"):
            self.gamma_2 = 1.0
        else:
            raise Exception("Check the time_integration_method_2 in project parameters file")

        # initialize mpm coupling utility class
        self.coupling_utility = KratosParticle.MPMTemporalCouplingUtility(self.model_sub_domain_1,
                                                                          self.model_sub_domain_2,
                                                                          self.time_step_ratio,
                                                                          self.time_step_2,
                                                                          self.gamma_1, self.gamma_2)

        # TODO read from parameters
        self.interface_criteria_normal = [1,0,0]
        self.interface_criteria_origin = [5,0,0]

        KratosMultiphysics.Logger.PrintInfo("::[MPMCoupledTimeSolver]:: ", "Construction finished.")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""
        {
            "analysis_type"   : "linear",
            "scheme_type"   : "newmark",
            "stress_update" : "usl",
            "newmark_beta"  : 0.25,
            "consistent_mass_matrix"  : false,
            "time_stepping"            : {
            "time_step"                          : 1,
            "timestep_ratio"                     : 1,
            "time_integration_method_1"          : "implicit",
            "time_integration_method_2"          : "implicit"
            }
        }""")
        this_defaults.AddMissingParameters(super(MPMCoupledTimeSolver, cls).GetDefaultSettings())
        return this_defaults

    ### Solver public functions

    def GetComputingSubModelPart(self, index):
        if index == 1:
            return self.model_sub_domain_1
        elif index == 2:
            return self.model_sub_domain_2


    def Initialize(self):
        # The particle solution strategy is created here if it does not already exist.
        particle_solution_strategy_1 = self._GetSolutionStrategy(1)
        particle_solution_strategy_1.SetEchoLevel(self.settings["echo_level"].GetInt())
        if self.gamma_1 == 1.0:
            self.grid_model_part.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT, True)
        else:
            self.grid_model_part.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT, False)

        particle_solution_strategy_2 = self._GetSolutionStrategy(2)
        particle_solution_strategy_2.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Generate material points
        self._GenerateMaterialPoint()

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","Solver is initialized correctly.")


    def AdvanceInTime(self, current_time):
        new_time = current_time + self.time_step_1

        self.grid_model_part.ProcessInfo[KratosMultiphysics.STEP] += self.time_step_ratio
        self.grid_model_part.CloneTimeStep(new_time)

        self.model_sub_domain_1.ProcessInfo[KratosMultiphysics.STEP] += self.time_step_ratio
        self.model_sub_domain_1.CloneTimeStep(new_time)

        return new_time


    def _is_time_step_solved(self, index):
        if index == 1:
            if is_integer(self.new_time/self.time_step_1):
                return True
            else:
                return False
        if index == 2:
            if is_integer(self.new_time/self.time_step_2):
                return True
            else:
                return False


    def InitializeSolutionStep(self):
        if self.gamma_1 == 1.0:
            self.model_sub_domain_1.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT, True)
        else:
            self.model_sub_domain_1.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT, False)
        self._SearchElement()
        #print('Initializing sd1')
        #print("Subdomain 1 time = ", self.model_sub_domain_1.ProcessInfo[KratosMultiphysics.TIME])
        self._GetSolutionStrategy(1).Initialize()
        self._GetSolutionStrategy(1).InitializeSolutionStep()
        self.coupling_utility.SetActiveNodes(self.model_sub_domain_1)
        self.coupling_utility.InitializeSubDomain1Coupling()


    def Predict(self):
        #print('Predicting sd1')
        self._GetSolutionStrategy(1).Predict()


    def SolveSolutionStep(self):
        #print('solving sd1')
        is_converged = self._GetSolutionStrategy(1).SolveSolutionStep()
        if self.gamma_1 == 0.5:
            K1 = self._GetSolutionStrategy(1).GetSystemMatrix()
            self.coupling_utility.StoreFreeVelocitiesSubDomain1(K1)
        else:
            self.coupling_utility.StoreFreeVelocitiesSubDomain1Explicit()

        # store interface velocities in coupling class vector
        if self.gamma_2 == 1.0:
            self.model_sub_domain_2.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT, True)
        else:
            self.model_sub_domain_2.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT, False)
        for j in range(1,self.time_step_ratio+1):
            #print('Advance SD2 time')
            time = self.model_sub_domain_2.ProcessInfo[KratosMultiphysics.TIME]
            time += self.time_step_2
            self.model_sub_domain_2.CloneTimeStep(time)
            self.model_sub_domain_2.ProcessInfo[KratosMultiphysics.STEP] += 1
            print("Subdomain 2 timestep", j, "of",self.time_step_ratio,
                  "(SD1 time = ",self.model_sub_domain_1.ProcessInfo[KratosMultiphysics.TIME], ", SD2 time = ",time,")")

            #print('Initializing sd2')
            self._GetSolutionStrategy(2).Initialize()
            self._GetSolutionStrategy(2).InitializeSolutionStep()
            self.coupling_utility.SetActiveNodes(self.model_sub_domain_2)
            #print('Predicting sd2')
            self._GetSolutionStrategy(2).Predict()
            #print('solving sd2')
            is_converged = self._GetSolutionStrategy(2).SolveSolutionStep()

            if self.gamma_2 == 0.5:
                K2 = self._GetSolutionStrategy(2).GetSystemMatrix()
                self.coupling_utility.CalculateCorrectiveLagrangianMultipliers(K2)
            else:
                self.coupling_utility.CalculateCorrectiveLagrangianMultipliersExplicit()
            print('finalizing sd2')
            self._GetSolutionStrategy(2).FinalizeSolutionStep()
            self._GetSolutionStrategy(2).Clear()
        if self.gamma_1 == 1.0:
            self.model_sub_domain_1.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT, True)
        else:
            self.model_sub_domain_1.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT, False)
        print('correcting sd1')
        self.coupling_utility.CorrectSubDomain1()
        return is_converged


    def FinalizeSolutionStep(self):
        print('finalizing sd1')
        self._GetSolutionStrategy(1).FinalizeSolutionStep()
        self._GetSolutionStrategy(1).Clear()
        #print('finalize sd1')


    def Check(self):
        self._GetSolutionStrategy(1).Check()
        self._GetSolutionStrategy(2).Check()


    def Clear(self):
        self._GetSolutionStrategy(1).Clear()
        self._GetSolutionStrategy(2).Clear()


    def _GetLinearSolver(self,index):
        if index == 1:
            if not hasattr(self, '_linear_solver_1'):
                self._linear_solver_1 = self._CreateLinearSolver()
                print("::[MPMCoupledTimeSolver]::    Created linear solver 1")
            return self._linear_solver_1
        if index == 2:
            if not hasattr(self, '_linear_solver_2'):
                self._linear_solver_2 = self._CreateLinearSolver()
                print("::[MPMCoupledTimeSolver]::    Created linear solver 2")
            return self._linear_solver_2

    def _GetConvergenceCriteria(self, index):
        if index == 1:
            if not hasattr(self, '_convergence_criterion_1'):
                self._convergence_criterion_1 = self._CreateConvergenceCriteria(index)
                print("::[MPMCoupledTimeSolver]::    Created convergence criteria 1")
            return self._convergence_criterion
        if index == 2:
            if not hasattr(self, '_convergence_criterion_2'):
                self._convergence_criterion_2 = self._CreateConvergenceCriteria(index)
                print("::[MPMCoupledTimeSolver]::    Created convergence criteria 2")
            return self._convergence_criterion

    def _GetBuilderAndSolver(self, index):
        if index == 1:
            if not hasattr(self, '_builder_and_solver_1'):
                self._builder_and_solver_1 = self._CreateBuilderAndSolver()
                print("::[MPMCoupledTimeSolver]::    Created builder and solver 1")
            return self._builder_and_solver
        if index == 2:
            if not hasattr(self, '_builder_and_solver_2'):
                self._builder_and_solver_2 = self._CreateBuilderAndSolver()
                print("::[MPMCoupledTimeSolver]::    Created builder and solver 2")
            return self._builder_and_solver


    def _GetSolutionScheme(self, index):
        if index == 1:
            if not hasattr(self, '_solution_scheme_1'):
                self._solution_scheme_1 = self._CreateSolutionScheme(index)
                print("::[MPMCoupledTimeSolver]::    Created scheme 1")
            return self._solution_scheme_1
        if index == 2:
            if not hasattr(self, '_solution_scheme_2'):
                self._solution_scheme_2 = self._CreateSolutionScheme(index)
                print("::[MPMCoupledTimeSolver]::    Created scheme 2")
            return self._solution_scheme_2


    def _GetSolutionStrategy(self, index):
        if index == 1:
            if not hasattr(self, '_solution_strategy_1'):
                self._solution_strategy_1 = self._CreateSolutionStrategy(index)
                print("::[MPMCoupledTimeSolver]::    Created strategy 1")
            return self._solution_strategy_1
        if index == 2:
            if not hasattr(self, '_solution_strategy_2'):
                self._solution_strategy_2 = self._CreateSolutionStrategy(index)
                print("::[MPMCoupledTimeSolver]::    Created strategy 2")
            return self._solution_strategy_2

    ### Solver protected functions

    def _GenerateMaterialPoint(self):
        pressure_dofs          = self.settings["pressure_dofs"].GetBool()
        axis_symmetric_flag    = self.settings["axis_symmetric_flag"].GetBool()

        # Assigning extra information to the main model part
        self.material_point_model_part.SetNodes(self.grid_model_part.GetNodes())
        self.material_point_model_part.ProcessInfo = self.grid_model_part.ProcessInfo
        self.material_point_model_part.SetBufferSize(self.grid_model_part.GetBufferSize())

        self.model_sub_domain_1.SetNodes(self.grid_model_part.GetNodes())
        self.model_sub_domain_2.SetNodes(self.grid_model_part.GetNodes())

        # Transfer temporal interface into sub domains
        self.model_sub_domain_1.CreateSubModelPart("temporal_interface")
        sub_domain_1_temporal_interface = self.model_sub_domain_1.GetSubModelPart("temporal_interface")
        self.model_sub_domain_2.CreateSubModelPart("temporal_interface")
        sub_domain_2_temporal_interface = self.model_sub_domain_2.GetSubModelPart("temporal_interface")
        grid_temporal_interface = self.grid_model_part.GetSubModelPart("temporal_interface")
        id_offset = 100 # TODO this should be the max node ID in the grid modelpart
        id_counter = 0
        #print('before interface add',self.model_sub_domain_2)
        for interface_node in grid_temporal_interface.Nodes:
            sub_domain_1_temporal_interface.AddNodes([interface_node.Id]) #add into sub domain 1 sub model part
            sub_domain_2_temporal_interface.AddNodes([interface_node.Id]) # TODO update
            #self.model_sub_domain_2.RemoveNode(interface_node.Id)
            #id_counter += 1
            #self.model_sub_domain_2.CreateNewNode(id_offset+id_counter,interface_node.X,interface_node.Y,interface_node.Z)
            #sub_domain_2_temporal_interface.AddNodes([id_offset+id_counter])
        #sub_domain_2_temporal_interface.Properties = self.initial_mesh_model_part.Properties
        #print('after interface add',self.model_sub_domain_2)
        #self._AddVariablesToModelPart(sub_domain_2_temporal_interface)
        #sub_domain_2_temporal_interface.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        print(sub_domain_1_temporal_interface)
        print(sub_domain_2_temporal_interface)


        self.model_sub_domain_1.SetBufferSize(self.grid_model_part.GetBufferSize())
        self.model_sub_domain_2.SetBufferSize(self.grid_model_part.GetBufferSize())

        print('----------- SETTING GRID NODES INTO SUB DOMAIN MODEL PARTS ----------------')

        # Generate MP Element and Condition
        KratosParticle.GenerateMaterialPointElement(self.grid_model_part, self.initial_mesh_model_part, self.material_point_model_part, pressure_dofs)
        KratosParticle.GenerateMaterialPointCondition(self.grid_model_part, self.initial_mesh_model_part, self.material_point_model_part)

        # TODO generalize interface definition
        for element in self.material_point_model_part.Elements:
            mp_coord = element.CalculateOnIntegrationPoints(KratosParticle.MP_COORD, self.material_point_model_part.ProcessInfo)
            for i in range(3):
                if self.interface_criteria_normal[i] == 1:
                    if mp_coord[0][i] < self.interface_criteria_origin[i]:
                        print('element ',element.Id,' with x = ',mp_coord[0][i],' added to subdomain 1')
                        self.model_sub_domain_1.AddElement(element,0)
                    else:
                        print('element ',element.Id,' with x = ',mp_coord[0][i],' added to subdomain 2')
                        self.model_sub_domain_2.AddElement(element,0)

        self.model_sub_domain_1.CreateSubModelPart("active_nodes")
        self.model_sub_domain_2.CreateSubModelPart("active_nodes")
        sub_domain_1_active_nodes = self.model_sub_domain_1.GetSubModelPart("active_nodes")
        sub_domain_2_active_nodes = self.model_sub_domain_2.GetSubModelPart("active_nodes")
        for node in self.grid_model_part.Nodes:
            if self.interface_criteria_normal[0] == 1:
                if node.X == self.interface_criteria_origin[0]:
                    sub_domain_1_active_nodes.AddNodes([node.Id])
                    sub_domain_2_active_nodes.AddNodes([node.Id])
                elif node.X < self.interface_criteria_origin[0]:
                    sub_domain_1_active_nodes.AddNodes([node.Id])
                else:
                    sub_domain_2_active_nodes.AddNodes([node.Id])
            elif self.interface_criteria_normal[1] == 1:
                if node.Y == self.interface_criteria_origin[1]:
                    sub_domain_1_active_nodes.AddNodes([node.Id])
                    sub_domain_2_active_nodes.AddNodes([node.Id])
                elif node.Y < self.interface_criteria_origin[1]:
                    sub_domain_1_active_nodes.AddNodes([node.Id])
                else:
                    sub_domain_2_active_nodes.AddNodes([node.Id])
            elif self.interface_criteria_normal[2] == 1:
                if node.Z == self.interface_criteria_origin[2]:
                    sub_domain_1_active_nodes.AddNodes([node.Id])
                    sub_domain_2_active_nodes.AddNodes([node.Id])
                elif node.Z < self.interface_criteria_origin[2]:
                    sub_domain_1_active_nodes.AddNodes([node.Id])
                else:
                    sub_domain_2_active_nodes.AddNodes([node.Id])
            else:
                raise Exception("Only simple interface definitions allowed so far.")
        print("\n\nsub domain 1 active nodes")
        for node in sub_domain_1_active_nodes.Nodes:
            print("node x = ",node.X)
        print("\n\nsub domain 2 active nodes")
        for node in sub_domain_2_active_nodes.Nodes:
            print("node x = ",node.X)


    def _SearchElement(self):
        for element in self.model_sub_domain_1.Elements:
            mp_coord = element.CalculateOnIntegrationPoints(KratosParticle.MP_COORD, self.material_point_model_part.ProcessInfo)
            for i in range(3):
                if self.interface_criteria_normal[i] == 1:
                    if not mp_coord[0][i] < self.interface_criteria_origin[i]:
                        self.model_sub_domain_2.AddElement(element)
                        self.model_sub_domain_1.RemoveElement(element)

        for element in self.model_sub_domain_2.Elements:
            mp_coord = element.CalculateOnIntegrationPoints(KratosParticle.MP_COORD, self.material_point_model_part.ProcessInfo)
            for i in range(3):
                if self.interface_criteria_normal[i] == 1:
                    if not mp_coord[0][i] >= self.interface_criteria_origin[i]:
                        self.model_sub_domain_1.AddElement(element)
                        self.model_sub_domain_2.RemoveElement(element)

        searching_alg_type = self.settings["element_search_settings"]["search_algorithm_type"].GetString()
        max_number_of_search_results = self.settings["element_search_settings"]["max_number_of_results"].GetInt()
        searching_tolerance          = self.settings["element_search_settings"]["searching_tolerance"].GetDouble()
        if (searching_alg_type == "bin_based"):
            KratosParticle.SearchElement(self.grid_model_part, self.model_sub_domain_1, max_number_of_search_results, searching_tolerance)
            KratosParticle.SearchElement(self.grid_model_part, self.model_sub_domain_2, max_number_of_search_results, searching_tolerance)
        else:
            err_msg  = "The requested searching algorithm \"" + searching_alg_type
            err_msg += "\" is not available for ParticleMechanicsApplication!\n"
            err_msg += "Available options are: \"bin_based\""
            raise Exception(err_msg)


    def _AddModelPartContainers(self):
        domain_size = self._GetDomainSize()
        if domain_size not in [2,3]:
            err_msg  = "The input \"domain_size\" is wrong!"
            err_msg += "Available options are: \"2\" or \"3\""
            raise Exception(err_msg)

        ## In MPM three model parts are needed
        # Material model part definition
        material_point_model_part_name = self.settings["model_part_name"].GetString()
        if not self.model.HasModelPart(material_point_model_part_name):
            self.material_point_model_part = self.model.CreateModelPart(material_point_model_part_name) # Equivalent to model_part3 in the old format
            self.material_point_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        # Initial material model part definition
        initial_mesh_model_part_name = "Initial_" + material_point_model_part_name
        if not self.model.HasModelPart(initial_mesh_model_part_name):
            self.initial_mesh_model_part = self.model.CreateModelPart(initial_mesh_model_part_name) #Equivalent to model_part2 in the old format
            self.initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        # Grid model part definition
        if not self.model.HasModelPart("Background_Grid"):
            self.grid_model_part = self.model.CreateModelPart("Background_Grid") #Equivalent to model_part1 in the old format
            self.grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        # model sub domain 1
        if not self.model.HasModelPart("model_sub_domain_1"):
            self.model_sub_domain_1 = self.model.CreateModelPart("model_sub_domain_1") #Equivalent to model_part1 in the old format
            self.model_sub_domain_1.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        # model sub domain 2
        if not self.model.HasModelPart("model_sub_domain_2"):
            self.model_sub_domain_2 = self.model.CreateModelPart("model_sub_domain_2") #Equivalent to model_part1 in the old format
            self.model_sub_domain_2.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)


    def _AddVariablesToModelPart(self, model_part):
        # Add displacements and reaction
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # Add specific variables for the problem conditions
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

        # MPM specific nodal variables
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MOMENTUM)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_INERTIA)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MIXED_TIME_LAGRANGE)

        # MPM dynamic variables
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)

        # Adding explicit variables
        self.grid_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        self.grid_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)

        # Add variables that the user defined in the ProjectParameters
        auxiliary_solver_utilities.AddVariables(model_part, self.settings["auxiliary_variables_list"])

        # Add variables for specific cases
        if self.settings["pressure_dofs"].GetBool():
            # add specific variables for the problem (pressure dofs)
            model_part.AddNodalSolutionStepVariable(KratosParticle.PRESSURE_REACTION)
            model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MPRESSURE)


    def _ModelPartReading(self):
        # reading the model part of the background grid
        if(self.settings["grid_model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.grid_model_part, self.settings["grid_model_import_settings"])
            self._AddTemporalInterfaceSubModelPartToGrid()
        else:
            raise Exception("Other input options are not implemented yet.")

        # reading the model part of the material point
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.initial_mesh_model_part, self.settings["model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")


    def _AddTemporalInterfaceSubModelPartToGrid(self):
        print('----------- adding temporal interface to grid  ----------------')
        self.grid_model_part.CreateSubModelPart("temporal_interface")
        temporal_interface = self.grid_model_part.GetSubModelPart("temporal_interface")
        for node in self.grid_model_part.Nodes:
            if self.interface_criteria_normal[0] == 1:
                if abs(node.X - self.interface_criteria_origin[0]) < self.tolerance:
                    temporal_interface.AddNodes([node.Id])
                    print('Added interface node ', node.Id, ' at x = ',node.X)
            elif self.interface_criteria_normal[1] == 1:
                if abs(node.Y - self.interface_criteria_origin[1]) < self.tolerance:
                    temporal_interface.AddNodes([node.Id])
            elif self.interface_criteria_normal[2] == 1:
                if abs(node.Z - self.interface_criteria_origin[2]) < self.tolerance:
                    temporal_interface.AddNodes([node.Id])
            else:
                raise Exception("Only simple interface definitions allowed so far.")


    def _AddDofsToModelPart(self, model_part):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, model_part)

        if self.settings["pressure_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosParticle.PRESSURE_REACTION, model_part)

        # Add dofs that the user defined in the ProjectParameters
        auxiliary_solver_utilities.AddDofs(model_part, self.settings["auxiliary_dofs_list"], self.settings["auxiliary_reaction_list"])

    def _GetDomainSize(self):
        if not hasattr(self, '_domain_size'):
            self._domain_size = self.settings["domain_size"].GetInt()
        return self._domain_size

    def _GetConvergenceCriteriaSettings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])

        return conv_params

    def _CreateConvergenceCriteria(self):
        convergence_criterion_parameters = self._GetConvergenceCriteriaSettings()
        if (convergence_criterion_parameters["convergence_criterion"].GetString() == "residual_criterion"):
            R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
            R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()
            convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            convergence_criterion.SetEchoLevel(convergence_criterion_parameters["echo_level"].GetInt())
        else:
            err_msg  = "The requested convergence criteria \"" + convergence_criterion_parameters["convergence_criterion"].GetString()
            err_msg += "\" is not supported for ParticleMechanicsApplication!\n"
            err_msg += "Available options are: \"residual_criterion\""
            raise Exception(err_msg)

        return convergence_criterion

    def _CreateSolutionScheme(self, index):
        grid_model_part = self.GetGridModelPart()
        domain_size = self._GetDomainSize()
        block_size  = domain_size

        is_implicit = True
        if index == 1:
            if (self.settings["time_stepping"]["time_integration_method_1"].GetString() == "explicit"):
                is_implicit = False
        if index == 2:
            if (self.settings["time_stepping"]["time_integration_method_2"].GetString() == "explicit"):
                is_implicit = False
        stress_update_option = 1
        grid_model_part.ProcessInfo.SetValue(KratosParticle.EXPLICIT_STRESS_UPDATE_OPTION, stress_update_option)
        grid_model_part.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT_CENTRAL_DIFFERENCE, False)
        # Setting the time integration schemes
        scheme_type = self.settings["scheme_type"].GetString()
        if(scheme_type == "newmark"):
            damp_factor_m = 0.0
            newmark_beta = self.settings["newmark_beta"].GetDouble()
        elif scheme_type == "forward_euler":
            stress_update_option = 1
            grid_model_part.ProcessInfo.SetValue(KratosParticle.EXPLICIT_STRESS_UPDATE_OPTION, stress_update_option)
            grid_model_part.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT_CENTRAL_DIFFERENCE, False)
        else:
            err_msg = "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"forward_euler\""
            raise Exception(err_msg)

        is_dynamic = self._IsDynamic()

        if (is_implicit):
            return KratosParticle.MPMResidualBasedBossakScheme( grid_model_part,
                                                                domain_size,
                                                                block_size,
                                                                damp_factor_m,
                                                                newmark_beta,
                                                                is_dynamic)
        else:
            stress_update_option = 1
            grid_model_part.ProcessInfo.SetValue(KratosParticle.EXPLICIT_STRESS_UPDATE_OPTION, stress_update_option)
            grid_model_part.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT_CENTRAL_DIFFERENCE, False)
            if index==1:
                self.model_sub_domain_1.ProcessInfo.SetValue(KratosParticle.EXPLICIT_STRESS_UPDATE_OPTION, stress_update_option)
            elif index == 2:
                self.model_sub_domain_2.ProcessInfo.SetValue(KratosParticle.EXPLICIT_STRESS_UPDATE_OPTION, stress_update_option)
            return KratosParticle.MPMExplicitScheme( grid_model_part)


    def _IsDynamic(self):
        return True

    def _CreateSolutionStrategy(self, index):
        analysis_type = self.settings["analysis_type"].GetString()
        is_consistent_mass_matrix = self.settings["consistent_mass_matrix"].GetBool()
        if is_consistent_mass_matrix:
            self.grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX, False)
            self.model_sub_domain_1.ProcessInfo.SetValue(KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX, False)
            self.model_sub_domain_2.ProcessInfo.SetValue(KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX, False)
        else:
            self.grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX, True)
            self.model_sub_domain_1.ProcessInfo.SetValue(KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX, True)
            self.model_sub_domain_2.ProcessInfo.SetValue(KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX, True)
        if analysis_type == "non_linear":
            solution_strategy = self._CreateNewtonRaphsonStrategy(index)
        elif analysis_type == 'linear':
            if index == 1:
                self.model_sub_domain_1.ProcessInfo.SetValue(KratosParticle.IGNORE_GEOMETRIC_STIFFNESS, True)
            if index == 2:
                self.model_sub_domain_2.ProcessInfo.SetValue(KratosParticle.IGNORE_GEOMETRIC_STIFFNESS, True)
            solution_strategy = self._CreateLinearStrategy(index)
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return solution_strategy

    def _CreateNewtonRaphsonStrategy(self, index):
        computing_model_part = self.GetComputingSubModelPart(index)
        solution_scheme = self._GetSolutionScheme(index)
        linear_solver = self._GetLinearSolver(index)
        convergence_criterion = self._GetConvergenceCriteria()
        builder_and_solver = self._GetBuilderAndSolver()
        reform_dofs_at_each_step = False ## hard-coded, but can be changed upon implementation
        return KratosParticle.MPMResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                                        solution_scheme,
                                                                        linear_solver,
                                                                        convergence_criterion,
                                                                        builder_and_solver,
                                                                        self.settings["max_iteration"].GetInt(),
                                                                        self.settings["compute_reactions"].GetBool(),
                                                                        reform_dofs_at_each_step,
                                                                        self.settings["move_mesh_flag"].GetBool())

    def _CreateLinearStrategy(self, index):
        is_implicit = True
        if index == 1:
            if (self.settings["time_stepping"]["time_integration_method_1"].GetString() == "explicit"):
                is_implicit = False
        if index == 2:
            if (self.settings["time_stepping"]["time_integration_method_2"].GetString() == "explicit"):
                is_implicit = False
        if is_implicit:
            computing_model_part = self.GetComputingSubModelPart(index)
            solution_scheme = self._GetSolutionScheme(index)
            linear_solver = self._GetLinearSolver(index)
            reform_dofs_at_each_step = False ## hard-coded, but can be changed upon implementation
            calc_norm_dx_flag = False ## hard-coded, but can be changed upon implementation
            return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                                  solution_scheme,
                                                                  linear_solver,
                                                                  self.settings["compute_reactions"].GetBool(),
                                                                  reform_dofs_at_each_step,
                                                                  calc_norm_dx_flag,
                                                                  self.settings["move_mesh_flag"].GetBool())
        else:
            computing_model_part = self.GetComputingSubModelPart(index)
            solution_scheme = self._GetSolutionScheme(index)
            reform_dofs_at_each_step = False ## hard-coded, but can be changed upon implementation
            move_mesh_flag = self.settings["move_mesh_flag"].GetBool()
            move_mesh_flag = False ## hard-coded
            return KratosParticle.MPMExplicitStrategy(computing_model_part,
                                                          solution_scheme,
                                                          self.settings["compute_reactions"].GetBool(),
                                                          reform_dofs_at_each_step,
                                                          move_mesh_flag)

    def _SetBufferSize(self):
        current_buffer_size = self.grid_model_part.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.grid_model_part.SetBufferSize(self.min_buffer_size)
        else:
            self.grid_model_part.SetBufferSize(current_buffer_size)

        current_buffer_size = self.initial_mesh_model_part.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.initial_mesh_model_part.SetBufferSize(self.min_buffer_size)
        else:
            self.initial_mesh_model_part.SetBufferSize(current_buffer_size)

    ### Solver private functions

    def __ExecuteCheckAndPrepare(self):
        # Specific active node and element check for particle MPM solver
        for node in self.grid_model_part.Nodes:
            if (node.Is(KratosMultiphysics.ACTIVE)):
                KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","WARNING: This grid node have been set active: ", node.Id)

        # Setting active initial elements
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, self.initial_mesh_model_part.Elements)

        # Read material property
        materials_imported = self.__ImportConstitutiveLaws()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","Constitutive law was successfully imported.")
        else:
            KratosMultiphysics.Logger.PrintWarning("::[MPMSolver]:: ","Constitutive law was not imported.")

        # Clone property of model_part2 to model_part3
        self.material_point_model_part.Properties = self.initial_mesh_model_part.Properties

    def __ImportConstitutiveLaws(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported
