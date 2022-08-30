from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import utilities
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory         # Import the FSI convergence accelerator factory
from KratosMultiphysics.FSIApplication import fsi_coupling_interface

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Import applications
import KratosMultiphysics.TrilinosApplication as KratosTrilinos
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Import base class file
from KratosMultiphysics.FSIApplication import partitioned_fsi_base_solver

def CreateSolver(model, project_parameters):
    return TrilinosPartitionedFSIBaseSolver(model, project_parameters)

class TrilinosPartitionedFSIBaseSolver(partitioned_fsi_base_solver.PartitionedFSIBaseSolver):
    def __init__(self, model, project_parameters):
        super(TrilinosPartitionedFSIBaseSolver, self).__init__(model, project_parameters)

    @classmethod
    def GetDefaultParameters(cls):
        """This function returns the default-settings used by this class
        """
        this_defaults = KratosMultiphysics.Parameters("""{
            "parallel_type": "MPI"
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def _GetEpetraCommunicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = self._CreateEpetraCommunicator()
        return self._epetra_communicator

    def _CreateEpetraCommunicator(self):
        return KratosTrilinos.CreateCommunicator()

    def _CreateConvergenceAccelerator(self):
        """ Create the MPI parallel convergence accelerator and assign it to the structure FSI coupling interface"""
        # Create the MPI parallel convergence accelerator
        convergence_accelerator = convergence_accelerator_factory.CreateTrilinosConvergenceAccelerator(
            self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
            self._GetEpetraCommunicator(),
            self.settings["coupling_settings"]["coupling_strategy_settings"])
        # Assign the new convergence accelerator to the structure FSI coupling interface
        self._GetFSICouplingInterfaceStructure().SetConvergenceAccelerator(convergence_accelerator)
        KratosMultiphysics.Logger.PrintInfo('TrilinosPartitionedFSIBaseSolver', 'Coupling strategy construction finished.')
        return convergence_accelerator

    def _CreateFSICouplingInterfaceStructure(self):
        """Create the structure FSI coupling interface

        Create the structure FSI coupling interface with a default MPI convergence accelerator
        The final convergence accelerator will be assigned to the structure FSI coupling interface in the _CreateConvergenceAccelerator
        This is required since the MPI parallel convergence accelerators require the residual minimization model part to be instantiated
        """
        # Set auxiliary settings
        aux_settings = KratosMultiphysics.Parameters(
        """{
            "model_part_name": "FSICouplingInterfaceStructure",
            "parent_model_part_name": "",
            "input_variable_list": ["SURFACE_LOAD"],
            "output_variable_list": ["DISPLACEMENT"]
        }""")
        aux_settings["parent_model_part_name"].SetString(self.interfaces_dict['structure'])

        # Construct the FSI coupling interface
        fsi_coupling_interface_structure = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            aux_settings,
            KratosTrilinos.TrilinosConvergenceAccelerator())

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Structure FSI coupling interface created')

        return fsi_coupling_interface_structure

    def _CreatePartitionedFSIUtilities(self):
        if self._GetDomainSize() == 2:
            return KratosTrilinos.TrilinosPartitionedFSIUtilitiesArray2D(self._GetEpetraCommunicator())
        elif self._GetDomainSize() == 3:
            return KratosTrilinos.TrilinosPartitionedFSIUtilitiesArray3D(self._GetEpetraCommunicator())
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self._GetDomainSize()))

    @classmethod
    def _CreateStructureToFluidInterfaceMapper(self, structure_interface, fluid_interface):
        mapper_params = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
        }""")
        structure_to_fluid_interface_mapper = KratosMapping.MapperFactory.CreateMPIMapper(
            structure_interface,
            fluid_interface,
            mapper_params)

        return structure_to_fluid_interface_mapper
