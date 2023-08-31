# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI
if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

class FSICouplingInterface:

    def __ValidateSettings(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name": "",
            "parent_model_part_name": "",
            "input_variable_list": [],
            "output_variable_list": [],
            "auxiliary_variable_list": []
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
        self.auxiliary_variable_list = settings["auxiliary_variable_list"]

        # Check the variables list sizes
        n_input_vars = self.input_variable_list.size()
        if n_input_vars != 1:
            err_msg = "\'input_variable_list\' size is {0}. Expected 1.".format(n_input_vars)
            raise Exception(err_msg)
        n_output_vars = self.output_variable_list.size()
        if (n_output_vars == 0) or (n_output_vars > 2):
            err_msg = "\'output_variable_list\' size is {0}. Expected 1 or 2.".format(n_output_vars)
            raise Exception(err_msg)

        # Check and get the output variable(s) data type
        aux_list = []
        for variable in self.output_variable_list.values():
            if KratosMultiphysics.KratosGlobals.Kernel.HasDoubleVariable(variable.GetString()):
                aux_list.append(True)
            elif KratosMultiphysics.KratosGlobals.Kernel.HasArrayVariable(variable.GetString()):
                aux_list.append(False)
            else:
                aux_list.append(None)
                break

        if aux_list.count(None) != 0:
            err_msg = "Non-expected variable type in \'output_variable_list\'. Only scalar and array variables are supported."
            raise Exception(err_msg)

        self.scalar_output = None
        if aux_list.count(True) == n_output_vars:
            self.scalar_output = True
        elif aux_list.count(False) == n_output_vars:
            self.scalar_output = False
        else:
            err_msg = "Scalar and array variable types are mixed in the \'output_variable_list\'."
            raise Exception(err_msg)

    def SetConvergenceAccelerator(self, convergence_accelerator):
        """Set the provided convergence accelerator to the current FSI coupling interface

        This function sets the convergence accelerator of the current FSI coupling interface
        This auxiliary method is understood to set the convergence accelerator in those situations
        in which its constructor requires the FSI coupling interface model part to be set (e.g. MPI)
        """
        self.convergence_accelerator = convergence_accelerator

    def GetInterfaceModelPart(self):
        if not hasattr(self, '_fsi_interface_model_part'):
            self._fsi_interface_model_part = self._create_fsi_interface_model_part()
        return self._fsi_interface_model_part

    def GetFatherModelPart(self):
        return self.model.GetModelPart(self.parent_model_part_name)

    def ComputeResidualVector(self):
        # Check that the output variable that defines the residual is unique and get it
        if (self.output_variable_list.size() == 1):
            output_variable_name = self.output_variable_list[0].GetString()
            output_variable = KratosMultiphysics.KratosGlobals.GetVariable(output_variable_name)
        else:
            err_msg = "ComputeResidualVector() can only be performed with a unique output variable.\n"
            err_msg += "Number of variables in \'output_variable_list\' is " + int(self.output_variable_list.size())
            raise Exception(err_msg)

        # Get the output variable from the father model part
        # Note that these are the current non-linear iteration unrelaxed values (\tilde{u}^{k+1})
        # These values will be used below to calculate the interface residual vector
        self.GetValuesFromFatherModelPart(output_variable)

        # Save the previously existent RELAXED_DISPLACEMENT in OLD_RELAXED_DISPLACEMENT before doing the update
        # These values will be used below to calculate the interface residual vector (u^{k})
        KratosMultiphysics.VariableUtils().CopyVariable(
            self._relaxed_variable,
            self._old_relaxed_variable,
            self.GetInterfaceModelPart().Nodes)

        # Compute the current non-linear iteration interface residual using the output variable
        # Note that the residual is computed as r^{k+1} = \tilde{u}^{k+1} - u^{k}
        self._output_variable_residual_vector = self._get_partitioned_fsi_utilities().SetUpInterfaceVector(self.GetInterfaceModelPart())
        self._get_partitioned_fsi_utilities().ComputeInterfaceResidualVector(
            self.GetInterfaceModelPart(),
            self._old_relaxed_variable,
            output_variable,
            self._residual_variable,
            self._output_variable_residual_vector,
            "nodal",
            KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM)

        # Return the current interface residual norm
        return self.GetInterfaceModelPart().ProcessInfo[KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM]

    def Update(self):
        # Set and fill the iteration value vector with the previous non-linear iteration relaxed values (u^{k})
        iteration_value_vector = self._get_partitioned_fsi_utilities().SetUpInterfaceVector(self.GetInterfaceModelPart())
        self._get_partitioned_fsi_utilities().InitializeInterfaceVector(
            self.GetInterfaceModelPart(),
            self._old_relaxed_variable,
            iteration_value_vector)

        # Compute the convergence accelerator correction
        self._get_convergence_accelerator().UpdateSolution(self._output_variable_residual_vector, iteration_value_vector)

        # Apply the corrected solution to the FSI coupling interface nodes
        self._get_partitioned_fsi_utilities().UpdateInterfaceValues(
            self.GetInterfaceModelPart(),
            self._relaxed_variable,
            iteration_value_vector)

    def UpdateDisplacement(self):
        # Set and fill the iteration value vector with the previous non-linear iteration relaxed values (u^{k})
        iteration_value_vector = self._get_partitioned_fsi_utilities().SetUpInterfaceVector(self.GetInterfaceModelPart())
        self._get_partitioned_fsi_utilities().InitializeInterfaceVector(
            self.GetInterfaceModelPart(),
            self._relaxed_variable,
            iteration_value_vector)

        # Check that the output variable that defines the residual is unique and get it (i.e. DISPLACEMENT)
        if (self.output_variable_list.size() == 1):
            output_variable_name = self.output_variable_list[0].GetString()
            output_variable = KratosMultiphysics.KratosGlobals.GetVariable(output_variable_name)
        else:
            err_msg = "ComputeResidualVector() can only be performed with a unique output variable.\n"
            err_msg += "Number of variables in \'output_variable_list\' is " + int(self.output_variable_list.size())
            raise Exception(err_msg)

        # Get the output variable from the father model part
        # Note that these are the current non-linear iteration unrelaxed values (\tilde{u}^{k+1})
        self.GetValuesFromFatherModelPart(output_variable)

        # Set and fill the displacement value vector with the current non-linear iteration values (\tilde{u}^{k+1})
        # Note that these are the current uncorrected displacement obtained with the previously corrected load
        iteration_value_u_vector = self._get_partitioned_fsi_utilities().SetUpInterfaceVector(self.GetInterfaceModelPart())
        self._get_partitioned_fsi_utilities().InitializeInterfaceVector(
            self.GetInterfaceModelPart(),
            output_variable,
            iteration_value_u_vector)

        # Set and fill the traction value vector with the current non-linear iteration values (f^{k+1})
        # Note that these are the values that were employed to obtain the current displacements
        iteration_value_f_vector = self._get_partitioned_fsi_utilities().SetUpInterfaceVector(self.GetInterfaceModelPart())
        self._get_partitioned_fsi_utilities().InitializeInterfaceVector(
            self.GetInterfaceModelPart(),
            self._traction_relaxed_variable,
            iteration_value_f_vector)

        # Compute the convergence accelerator correction
        self._get_convergence_accelerator().UpdateSolutionRight(
            iteration_value_f_vector,
            iteration_value_u_vector,
            iteration_value_vector)

        # Apply the corrected solution to the FSI coupling interface nodes
        self._get_partitioned_fsi_utilities().UpdateInterfaceValues(
            self.GetInterfaceModelPart(),
            self._relaxed_variable,
            iteration_value_vector)

    def UpdateTraction(self):
        # Set and fill the iteration value vector with the previous non-linear iteration relaxed values (f^{k})
        iteration_value_vector = self._get_partitioned_fsi_utilities().SetUpInterfaceVector(self.GetInterfaceModelPart())
        self._get_partitioned_fsi_utilities().InitializeInterfaceVector(
            self.GetInterfaceModelPart(),
            self._traction_relaxed_variable,
            iteration_value_vector)

        # Set and fill the displacement value vector with the current non-linear iteration values (\tilde{f}^{k+1})
        if self.input_variable_list.size() == 1:
            input_variable_name = self.input_variable_list[0].GetString()
            input_variable = KratosMultiphysics.KratosGlobals.GetVariable(input_variable_name)
        else:
            err_msg = "Input variable list has more than one variable. One is expected for the IBQN update."
            raise Exception(err_msg)

        iteration_value_f_vector = self._get_partitioned_fsi_utilities().SetUpInterfaceVector(self.GetInterfaceModelPart())
        self._get_partitioned_fsi_utilities().InitializeInterfaceVector(
            self.GetInterfaceModelPart(),
            input_variable,
            iteration_value_f_vector)

        # Set and fill the traction value vector with the current non-linear iteration values (u^{k+1})
        # Note that these are the values that were employed to obtain the current tractions
        iteration_value_u_vector = self._get_partitioned_fsi_utilities().SetUpInterfaceVector(self.GetInterfaceModelPart())
        self._get_partitioned_fsi_utilities().InitializeInterfaceVector(
            self.GetInterfaceModelPart(),
            self._relaxed_variable,
            iteration_value_u_vector)

        # Compute the convergence accelerator correction
        self._get_convergence_accelerator().UpdateSolutionLeft(
            iteration_value_u_vector,
            iteration_value_f_vector,
            iteration_value_vector)

        # Apply the corrected solution to the FSI coupling interface nodes
        self._get_partitioned_fsi_utilities().UpdateInterfaceValues(
            self.GetInterfaceModelPart(),
            self._traction_relaxed_variable,
            iteration_value_vector)

    def UpdatePosition(self, update_variable = KratosMultiphysics.DISPLACEMENT):
        KratosMultiphysics.VariableUtils().UpdateCurrentPosition(
            self.GetInterfaceModelPart().Nodes,
            update_variable)

    def GetValuesFromFatherModelPart(self, variable):
        buffer_step = 0
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
            variable,
            self.GetFatherModelPart(),
            self.GetInterfaceModelPart(),
            buffer_step)
        #TODO: Remove as soon as the CopyModelPartNodalVar synchronizes internally
        self.GetInterfaceModelPart().GetCommunicator().SynchronizeVariable(variable)

    def TransferValuesToFatherModelPart(self, variable):
        buffer_step = 0
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
            variable,
            self.GetInterfaceModelPart(),
            self.GetFatherModelPart(),
            buffer_step)
        #TODO: Remove as soon as the CopyModelPartNodalVar synchronizes internally
        self.GetFatherModelPart().GetCommunicator().SynchronizeVariable(variable)

    def _create_fsi_interface_model_part(self):
        # Add the FSI coupling interface to the model
        self._fsi_interface_model_part = self.model.CreateModelPart(self.model_part_name)

        # Set required ProcessInfo data
        domain_size = self.GetFatherModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self._fsi_interface_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = domain_size

        # Add the required variables to the FSI coupling interface model part
        for variable_name in self.input_variable_list.values():
            input_variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name.GetString())
            self._fsi_interface_model_part.AddNodalSolutionStepVariable(input_variable)
        for variable_name in self.output_variable_list.values():
            output_variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name.GetString())
            self._fsi_interface_model_part.AddNodalSolutionStepVariable(output_variable)
        for variable_name in self.auxiliary_variable_list.values():
            auxiliary_variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name.GetString())
            self._fsi_interface_model_part.AddNodalSolutionStepVariable(auxiliary_variable)

        # Add the variables required for the residual minimization
        if self._get_convergence_accelerator() is not None:
            # Get the required variables
            if self.scalar_output:
                self._relaxed_variable = KratosMultiphysics.RELAXED_SCALAR
                self._old_relaxed_variable = KratosMultiphysics.OLD_RELAXED_SCALAR
                self._residual_variable = KratosMultiphysics.SCALAR_INTERFACE_RESIDUAL
                if self._get_convergence_accelerator().IsBlockNewton():
                    self._traction_relaxed_variable = KratosMultiphysics.RELAXED_SCALAR_TRACTION
                    self._old_traction_relaxed_variable = KratosMultiphysics.OLD_RELAXED_SCALAR_TRACTION
                    self._traction_residual_variable = KratosMultiphysics.SCALAR_TRACTION_INTERFACE_RESIDUAL
            else:
                self._relaxed_variable = KratosMultiphysics.RELAXED_DISPLACEMENT
                self._old_relaxed_variable = KratosMultiphysics.OLD_RELAXED_DISPLACEMENT
                self._residual_variable = KratosMultiphysics.FSI_INTERFACE_RESIDUAL
                if self._get_convergence_accelerator().IsBlockNewton():
                    self._traction_relaxed_variable = KratosMultiphysics.RELAXED_TRACTION
                    self._old_traction_relaxed_variable = KratosMultiphysics.OLD_RELAXED_TRACTION
                    self._traction_residual_variable = KratosMultiphysics.TRACTION_INTERFACE_RESIDUAL

            # Add the required variables to the nodal database
            self._fsi_interface_model_part.AddNodalSolutionStepVariable(self._relaxed_variable)
            self._fsi_interface_model_part.AddNodalSolutionStepVariable(self._old_relaxed_variable)
            self._fsi_interface_model_part.AddNodalSolutionStepVariable(self._residual_variable)
            if self._get_convergence_accelerator().IsBlockNewton():
                self._fsi_interface_model_part.AddNodalSolutionStepVariable(self._traction_relaxed_variable)
                self._fsi_interface_model_part.AddNodalSolutionStepVariable(self._old_traction_relaxed_variable)
                self._fsi_interface_model_part.AddNodalSolutionStepVariable(self._traction_residual_variable)

        # If parallel, add the PARTITION_INDEX variable
        if self.GetFatherModelPart().IsDistributed():
            self._fsi_interface_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        # Set the FSI coupling interface entities (nodes and elements)
        self._get_partitioned_fsi_utilities().CreateCouplingSkin(
            self.GetFatherModelPart(),
            self._fsi_interface_model_part)

        # Create the communication plan for the current coupling interface
        # Note that we retrieve the fill communicator from the ParallelEnvironment so nothing would be done if non MPI
        fill_communicator = KratosMultiphysics.ParallelEnvironment.CreateFillCommunicatorFromGlobalParallelism(
            self._fsi_interface_model_part,
            self._fsi_interface_model_part.GetCommunicator().GetDataCommunicator())
        fill_communicator.Execute()

        return self._fsi_interface_model_part

    def _get_convergence_accelerator(self):
        return self.convergence_accelerator

    def _get_partitioned_fsi_utilities(self):
        if not hasattr(self, '_partitioned_fsi_utilities'):
            self._partitioned_fsi_utilities = self._create_partitioned_fsi_utilities()
        return self._partitioned_fsi_utilities

    def _create_partitioned_fsi_utilities(self):
        domain_size = self.GetFatherModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if not self.GetFatherModelPart().IsDistributed():
            if domain_size == 2:
                if self.scalar_output:
                    return KratosFSI.PartitionedFSIUtilitiesDouble2D()
                else:
                    return KratosFSI.PartitionedFSIUtilitiesArray2D()
            elif domain_size == 3:
                if self.scalar_output:
                    return KratosFSI.PartitionedFSIUtilitiesDouble3D()
                else:
                    return KratosFSI.PartitionedFSIUtilitiesArray3D()
            else:
                raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))
        else:
            self._epetra_communicator = KratosTrilinos.CreateCommunicator()
            if domain_size == 2:
                if self.scalar_output:
                    return KratosTrilinos.TrilinosPartitionedFSIUtilitiesDouble2D(self._epetra_communicator)
                else:
                    return KratosTrilinos.TrilinosPartitionedFSIUtilitiesArray2D(self._epetra_communicator)
            elif domain_size == 3:
                if self.scalar_output:
                    return KratosTrilinos.TrilinosPartitionedFSIUtilitiesDouble3D(self._epetra_communicator)
                else:
                    return KratosTrilinos.TrilinosPartitionedFSIUtilitiesArray3D(self._epetra_communicator)
            else:
                raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))
