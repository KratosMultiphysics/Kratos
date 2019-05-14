from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI

class FSICouplingInterface():

    def _ValidateSettings(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name": "",
            "father_model_part_name": "",
            "input_variable_name": "",
            "output_variable_name": ""
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        if settings["model_part_name"].GetString() == "":
            raise Exception("Provided \'model_part_name\' is empty.")
        if settings["father_model_part_name"].GetString() == "":
            raise Exception("Provided \'father_model_part_name\' is empty.")
        if settings["input_variable_name"].GetString() == "":
            raise Exception("Provided \'input_variable_name\' is empty.")
        if settings["output_variable_name"].GetString() == "":
            raise Exception("Provided \'output_variable_name\' is empty.")

        return settings

    def __init__(self, model, settings, convergence_accelerator = None):

        # Validate settings
        settings = self._ValidateSettings(settings)

        # Set required member variables
        self.model = model
        self.convergence_accelerator = convergence_accelerator
        self.model_part_name = settings["model_part_name"].GetString()
        self.father_model_part_name = settings["father_model_part_name"].GetString()
        self.input_variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["input_variable_name"].GetString())
        self.output_variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["output_variable_name"].GetString())

    def GetInterfaceModelPart(self):
        if not hasattr(self, '_fsi_interface_model_part'):
            self._fsi_interface_model_part = self._create_fsi_interface_model_part()
        return self._fsi_interface_model_part

    def GetFatherModelPart(self):
        return self.model.GetModelPart(self.father_model_part_name)

    def Update(self):
        # Set the interface (residual size)
        residual_size = self._get_partitioned_fsi_utilities().GetInterfaceResidualSize(self.GetInterfaceModelPart())

        # Set and fill the iteration value vector with the current output variable values
        # Note that these are the ones of the previous non-linear iteration
        iteration_value = KratosMultiphysics.Vector(residual_size)
        self._get_partitioned_fsi_utilities().InitializeInterfaceVector(
            self.GetInterfaceModelPart(),
            self.output_variable,
            iteration_value)

        # Get the output variable from the father model part
        # Note that these is to obtain the current non-linear iteration unrelaxed values
        self.GetValuesFromFatherModelPart(self.output_variable)

        # Compute the interface residual using the output variable
        output_variable_residual = KratosMultiphysics.Vector(residual_size)
        self._get_partitioned_fsi_utilities().ComputeInterfaceResidualVector(
            self.GetInterfaceModelPart(),
            self.output_variable,
            KratosMultiphysics.VECTOR_PROJECTED,
            KratosMultiphysics.FSI_INTERFACE_RESIDUAL,
            output_variable_residual,
            "nodal",
            KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM)

        # Compute the convergence accelerator correction
        self._get_convergence_accelerator().UpdateSolution(output_variable_residual, iteration_value)

        # Apply the corrected solution to the FSI coupling interface
        self._get_partitioned_fsi_utilities().UpdateInterfaceValues(
            self.GetInterfaceModelPart(),
            self.output_variable,
            iteration_value)

        # And save it for the next residual computation
        self._get_partitioned_fsi_utilities().UpdateInterfaceValues(
            self.GetInterfaceModelPart(),
            KratosMultiphysics.VECTOR_PROJECTED,
            iteration_value)

        # Return the current interface residual norm
        return self.GetInterfaceModelPart().ProcessInfo[KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM]

    def UpdatePosition(self):
        for node in self.GetInterfaceModelPart().Nodes:
            aux_disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            node.X = node.X0 + aux_disp[0]
            node.Y = node.Y0 + aux_disp[1]
            node.Z = node.Z0 + aux_disp[2]

    def GetValuesFromFatherModelPart(self, variable):
        buffer_step = 0
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
            variable,
            self.GetFatherModelPart(),
            self.GetInterfaceModelPart(),
            buffer_step)

    def TransferValuesToFatherModelPart(self, variable):
        buffer_step = 0
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
            variable,
            self.GetInterfaceModelPart(),
            self.GetFatherModelPart(),
            buffer_step)

    def _create_fsi_interface_model_part(self):
        # Add the FSI coupling interface to the model
        self._fsi_interface_model_part = self.model.CreateModelPart(self.model_part_name)

        # Add the required variables to the FSI coupling interface model part
        self._fsi_interface_model_part.AddNodalSolutionStepVariable(self.input_variable)
        self._fsi_interface_model_part.AddNodalSolutionStepVariable(self.output_variable)
        self._fsi_interface_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_PROJECTED)
        self._fsi_interface_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FSI_INTERFACE_RESIDUAL)

        # Set the FSI coupling interface entities (nodes and elements)
        self._get_partitioned_fsi_utilities().CreateCouplingElementBasedSkin(
            self.GetFatherModelPart(),
            self._fsi_interface_model_part)

        return self._fsi_interface_model_part

    def _get_convergence_accelerator(self):
        return self.convergence_accelerator

    def _get_partitioned_fsi_utilities(self):
        if not hasattr(self, '_partitioned_fsi_utilities'):
            self._partitioned_fsi_utilities = self._create_partitioned_fsi_utilities()
        return self._partitioned_fsi_utilities

    def _create_partitioned_fsi_utilities(self):
        domain_size = self.GetFatherModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size == 2:
            return KratosFSI.PartitionedFSIUtilitiesArray2D()
        elif domain_size == 3:
            return KratosFSI.PartitionedFSIUtilitiesArray3D()
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))
