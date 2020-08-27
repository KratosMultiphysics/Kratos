# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver
#from KratosMultiphysics.CoSimulationApplication.custom_utilities import FetiDynamicCouplingUtilities


def Create(settings, models, solver_name):
    return FetiDynamicCoupledSolver(settings, models, solver_name)

class FetiDynamicCoupledSolver(CoSimulationCoupledSolver):
    def __init__(self, settings, models, solver_name):
        super().__init__(settings, models, solver_name)

        self.is_initialized = False


    def SolveSolutionStep(self):
        if self.is_initialized == False:
            self.__InitializeFetiMethod()

        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.InitializeCouplingIteration()

        solver_index = 0
        for solver_name, solver in self.solver_wrappers.items():
            self._SynchronizeInputData(solver_name)
            solver.SolveSolutionStep()

            system_matrix = solver.get_mechanical_solution_strategy().GetSystemMatrix()
            self.feti_coupling_solver.SetEffectiveStiffnessMatrices(system_matrix,solver_index)

            self._SynchronizeOutputData(solver_name)

            solver_index += 1

        self.feti_coupling_solver.EquilibrateDomains()

        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.FinalizeCouplingIteration()

        return True

    def __InitializeFetiMethod(self):
        print('\n\n\n ========= __InitializeFetiMethod \n\n\n')

        # get mapper parameters
        self.mapper_parameters = self.data_transfer_operators_dict["mapper"].settings["mapper_settings"]

        #get solvers
        self.solver_wrappers_vector = []
        for solver_name, solver in self.solver_wrappers.items():
            self.solver_wrappers_vector.append(solver)

        # get mapper origin and destination modelparts
        origin_modelpart_name = self.mapper_parameters["modeler_parameters"]["origin_interface_sub_model_part_name"].GetString()
        self.model_part_origin = self.solver_wrappers_vector[0].model.GetModelPart(origin_modelpart_name)
        print(self.model_part_origin)

        destination_modelpart_name = self.mapper_parameters["modeler_parameters"]["destination_interface_sub_model_part_name"].GetString()
        self.model_part_destination = self.solver_wrappers_vector[1].model.GetModelPart(destination_modelpart_name)
        print(self.model_part_destination)

        # manually create mapper
        mapper_create_fct = KratosMapping.MapperFactory.CreateMapper
        self.mapper = mapper_create_fct(self.model_part_origin, self.model_part_destination, self.mapper_parameters.Clone())







        # get interface modelparts created by the mapper modeler
        self.modelpart_interface_origin = self.mapper.GetInterfaceModelPart(0)
        self.modelpart_interface_destination = self.mapper.GetInterfaceModelPart(1)
        print(self.modelpart_interface_origin)
        print(self.modelpart_interface_destination)

        # Get time integration parameters
        origin_newmark_beta = self.settings["origin_newmark_beta"].GetDouble()
        origin_newmark_gamma = self.settings["origin_newmark_gamma"].GetDouble()
        destination_newmark_beta = self.settings["destination_newmark_beta"].GetDouble()
        destination_newmark_gamma = self.settings["destination_newmark_gamma"].GetDouble()

        # Create feti class instance
        self.feti_coupling = FetiDynamicCouplingUtilities(
            self.modelpart_interface_origin,
            self.modelpart_interface_destination,
            origin_newmark_beta, origin_newmark_gamma,
            destination_newmark_beta, destination_newmark_gamma)


        # set the mapper

        self.feti_coupling.SetMappingMatrix(self.mapper.pGetMappingMatrix())

        self.is_initialized = True

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
