# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.CoSimulationApplication as CoSim

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver

# Other imports
from KratosMultiphysics import auxiliary_solver_utilities
from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory


def Create(settings, models, solver_name):
    return FetiDynamicCoupledSolver(settings, models, solver_name)

class FetiDynamicCoupledSolver(CoSimulationCoupledSolver):
    def __init__(self, settings, models, solver_name):
        super().__init__(settings, models, solver_name)

        #get solvers and add lagrange multiplier solution variable
        self.solver_wrappers_vector = []
        for solver_name, solver in self.solver_wrappers.items():
            self.solver_wrappers_vector.append(solver)
            structure = solver.model.GetModelPart("Structure")
            structure.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)

        self.is_initialized = False


    def SolveSolutionStep(self):
        if self.is_initialized == False:
            self.__InitializeFetiMethod()

        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.InitializeCouplingIteration()

        solver_index = 0
        for solver_name, solver in self.solver_wrappers.items():
            solver.SolveSolutionStep()

            system_matrix = solver.GetSolverStrategy().GetSystemMatrix()
            self.feti_coupling.SetEffectiveStiffnessMatrices(system_matrix,solver_index)

            solver_index += 1

        self.feti_coupling.EquilibrateDomains()

        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.FinalizeCouplingIteration()

        return True

    def __InitializeFetiMethod(self):
        print('\n\n\n ========= __InitializeFetiMethod \n\n\n')

        # get mapper parameters
        self.mapper_parameters = self.data_transfer_operators_dict["mapper"].settings["mapper_settings"]
        mapper_type = self.mapper_parameters["mapper_type"].GetString()


        # get mapper origin and destination modelparts
        origin_modelpart_name = self.mapper_parameters["modeler_parameters"]["origin_interface_sub_model_part_name"].GetString()
        self.model_part_origin_interface = self.solver_wrappers_vector[0].model.GetModelPart(origin_modelpart_name)

        destination_modelpart_name = self.mapper_parameters["modeler_parameters"]["destination_interface_sub_model_part_name"].GetString()
        self.model_part_destination_interface = self.solver_wrappers_vector[1].model.GetModelPart(destination_modelpart_name)

        # manually create mapper
        mapper_create_fct = KratosMapping.MapperFactory.CreateMapper
        self.mapper = mapper_create_fct(self.model_part_origin_interface, self.model_part_destination_interface, self.mapper_parameters.Clone())


        # get interface modelparts created by the mapper modeler
        self.modelpart_interface_origin_from_mapper = self.mapper.GetInterfaceModelPart(0)
        self.modelpart_interface_destination_from_mapper = self.mapper.GetInterfaceModelPart(1)


        # Get time integration parameters
        origin_newmark_beta = self.settings["origin_newmark_beta"].GetDouble()
        origin_newmark_gamma = self.settings["origin_newmark_gamma"].GetDouble()
        destination_newmark_beta = self.settings["destination_newmark_beta"].GetDouble()
        destination_newmark_gamma = self.settings["destination_newmark_gamma"].GetDouble()


        # create the solver
        linear_solver = self._CreateLinearSolver()


        # Create feti class instance
        self.feti_coupling = CoSim.FetiDynamicCouplingUtilities(
            self.modelpart_interface_origin_from_mapper,
            self.modelpart_interface_destination_from_mapper,
            origin_newmark_beta, origin_newmark_gamma,
            destination_newmark_beta, destination_newmark_gamma)

        # The origin and destination interfaces from the mapper submitted above are both
        # stored within the origin modelpart. Now we submit the 'original' origin and destination
        # interface model parts stored on the origin and destination models to get access to the
        # origin and destination models.
        self.feti_coupling.SetOriginAndDestinationDomainsWithInterfaceModelParts(
            self.model_part_origin_interface,
            self.model_part_destination_interface)

        self.feti_coupling.SetLinearSolver(linear_solver)


        # set the mapper
        if mapper_type == "coupling_geometry":
            p_mapping_matrix = self.mapper.pGetDenseMappingMatrix()
            self.feti_coupling.SetMappingMatrix(p_mapping_matrix)

        self.is_initialized = True

    def _CreateLinearSolver(self):
        return linear_solver_factory.CreateFastestAvailableDirectLinearSolver()
        #linear_solver_configuration = self.settings["linear_solver_settings"]
        #if linear_solver_configuration.Has("solver_type"): # user specified a linear solver
        #    return linear_solver_factory.ConstructSolver(linear_solver_configuration)
        #else:
        #    KratosMultiphysics.Logger.PrintInfo('::[MPMSolver]:: No linear solver was specified, using fastest available solver')
        #    return linear_solver_factory.CreateFastestAvailableDirectLinearSolver()

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "origin_newmark_beta" : -1.0,
            "origin_newmark_gamma" : -1.0,
            "destination_newmark_beta" : -1.0,
            "destination_newmark_gamma" : -1.0
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultSettings())

        return this_defaults
