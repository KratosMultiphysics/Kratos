from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI

class FSICouplingInterface():

    def __ValidateSettings(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name": "",
            "parent_model_part_name": "",
            "input_variable_list": [],
            "output_variable_list": []
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        if settings["model_part_name"].GetString() == "":
            raise Exception("Provided \'model_part_name\' is empty.")
        if settings["parent_model_part_name"].GetString() == "":
            raise Exception("Provided \'parent_model_part_name\' is empty.")
        if settings["input_variable_list"].size() == 0:
            raise Exception("Provided \'input_variable_list\' is empty.")
        if settings["output_variable_list"].size() == 0:
            raise Exception("Provided \'output_variable_list\' is empty.")

        return settings

    def __init__(self, model, settings, convergence_accelerator = None):
        # Validate settings
        settings = self.__ValidateSettings(settings)

        # Set required member variables
        self.model = model
        self.convergence_accelerator = convergence_accelerator
        self.model_part_name = settings["model_part_name"].GetString()
        self.parent_model_part_name = settings["parent_model_part_name"].GetString()
        self.input_variable_list = settings["input_variable_list"]
        self.output_variable_list = settings["output_variable_list"]

    def GetInterfaceModelPart(self):
        if not hasattr(self, '_fsi_interface_model_part'):
            self._fsi_interface_model_part = self._create_fsi_interface_model_part()
        return self._fsi_interface_model_part

    def GetFatherModelPart(self):
        return self.model.GetModelPart(self.parent_model_part_name)

    def Update(self):
        # Check that the output variable that defines the residual is unique and get it
        if (self.output_variable_list.size() == 1):
            output_variable_name = self.output_variable_list[0].GetString()
            output_variable = KratosMultiphysics.KratosGlobals.GetVariable(output_variable_name)
        else:
            err_msg = "Update() can only be performed with a unique output variable.\n"
            err_msg = "Number of variables in \'output_variable_list\' is " + int(self.output_variable_list.size())
            raise Exception(err_msg)

        # Get the output variable from the father model part
        # Note that this are the current non-linear iteration unrelaxed values (\tilde{u}^{k+1})
        self.GetValuesFromFatherModelPart(output_variable)

        # Save the previously existent RELAXED_DISPLACEMENT in OLD_RELAXED_DISPLACEMENT before doing the update
        for node in self.GetInterfaceModelPart().Nodes:
            relaxed_disp = node.GetSolutionStepValue(KratosMultiphysics.RELAXED_DISPLACEMENT)
            node.SetSolutionStepValue(KratosMultiphysics.OLD_RELAXED_DISPLACEMENT, relaxed_disp)

        # Get the interface residual size
        residual_size = self._get_partitioned_fsi_utilities().GetInterfaceResidualSize(self.GetInterfaceModelPart())

        # Set and fill the iteration value vector with the previous non-linear iteration relaxed values (u^{k})
        iteration_value_vector = KratosMultiphysics.Vector(residual_size)
        self._get_partitioned_fsi_utilities().InitializeInterfaceVector(
            self.GetInterfaceModelPart(),
            KratosMultiphysics.OLD_RELAXED_DISPLACEMENT,
            iteration_value_vector)

        # Compute the current non-linear iteration interface residual using the output variable
        # Note that the residual is computed as r^{k+1} = \tilde{u}^{k+1} - u^{k}
        output_variable_residual_vector = KratosMultiphysics.Vector(residual_size)
        self._get_partitioned_fsi_utilities().ComputeInterfaceResidualVector(
            self.GetInterfaceModelPart(),
            output_variable,
            KratosMultiphysics.OLD_RELAXED_DISPLACEMENT,
            KratosMultiphysics.FSI_INTERFACE_RESIDUAL,
            output_variable_residual_vector,
            "nodal",
            KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM)

        # Compute the convergence accelerator correction
        self._get_convergence_accelerator().UpdateSolution(output_variable_residual_vector, iteration_value_vector)

        # Apply the corrected solution to the FSI coupling interface nodes
        self._get_partitioned_fsi_utilities().UpdateInterfaceValues(
            self.GetInterfaceModelPart(),
            KratosMultiphysics.RELAXED_DISPLACEMENT,
            iteration_value_vector)

        # Return the current interface residual norm
        return self.GetInterfaceModelPart().ProcessInfo[KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM]

    def UpdatePosition(self):
        KratosMultiphysics.VariableUtils().UpdateCurrentPosition(self.GetInterfaceModelPart().Nodes)

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
        for variable_name in self.input_variable_list:
            input_variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name.GetString())
            self._fsi_interface_model_part.AddNodalSolutionStepVariable(input_variable)
        for variable_name in self.output_variable_list:
            output_variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name.GetString())
            self._fsi_interface_model_part.AddNodalSolutionStepVariable(output_variable)
        self._fsi_interface_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.RELAXED_DISPLACEMENT)
        self._fsi_interface_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.OLD_RELAXED_DISPLACEMENT)
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
