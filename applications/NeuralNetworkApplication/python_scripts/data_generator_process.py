import KratosMultiphysics as KM
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
import KratosMultiphysics.HDF5Application as KratosHDF5


import numpy as np
import operator as op

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DataGeneratorProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class DataGeneratorProcess(KM.Process):
    """This class generates a dataset for Neural Network training and testing

    Public member variables:
    model -- the container of the different model parts.
    settings -- Kratos parameters containing process settings.
    """

    def __init__(self, model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        KM.Process.__init__(self) # calling the baseclass constructor

        default_settings = KM.Parameters("""{
            "help"                     : "This process generates the data series for the neural network training and testing",
            "model_part_name"          : "",
            "interval"                 : [0, 1e30],
            "write_output_file"        : true,
            "input_variables"          : [""],
            "input_sources"            : [""],
            "input_model_part"         : "",
            "output_variables"         : [""],
            "output_sources"            : [""],
            "output_format"            : "ascii",
            "perturbate_variables"     : [""],
            "random_distribution"      : [""],
            "random_parameters"        : [],
            "print_format"             : ".8f",
            "training_file_settings": {
                    "input_file_name"  : "training_in",
                    "output_file_name" : "training_out",
                    "file_format"      : "dat",
                    "output_path"      : "data/"
            }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        # Detect 'End' as a tag and replace it by a large number
        if(settings.Has('interval')):
            if(settings['interval'][1].IsString()):
                if(settings['interval'][1].GetString() == 'End'):
                    settings['interval'][1].SetDouble(1e30)
                else:
                    raise Exception('The second value of interval can be \'End\' or a number, interval currently:' +
                                    settings['interval'].PrettyPrintJsonString())

        self.write_output_file = settings["write_output_file"].GetBool()
        perturbate_names = settings["perturbate_variables"]
        perturbate_variable_names = [ perturbate_names[i].GetString() for i in range( perturbate_names.size() ) ]
        self.perturbate_variables = [ KM.KratosGlobals.GetVariable( var ) for var in perturbate_variable_names ]
        input_dist_names = settings["random_distribution"]
        self.random_distribution = [ input_dist_names[i].GetString() for i in range( input_dist_names.size() ) ]
        

        # getting the list of distribution parameters for each variable
        random_parameters_matrix = settings["random_parameters"].GetMatrix()
        self.random_parameters =[]
        for distribution in range(random_parameters_matrix.Size1()):
            distribution_parameters = []
            for parameter in range(random_parameters_matrix.Size2()):
                distribution_parameters.append(random_parameters_matrix[distribution,parameter])
            self.random_parameters.append(distribution_parameters)
        self.interval = settings["interval"].GetVector()
        self.format = settings["print_format"].GetString()
        self.output_format = settings["output_format"].GetString()

        # getting the ModelPart from the Model
        model_part_name = settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        self.model_part = model[model_part_name]

         # getting the input ModelPart from the Model
        input_model_part_name = settings["input_model_part"].GetString()
        if input_model_part_name == "":
            raise Exception('No "input_model_part" was specified!')
        self.input_model_part = model[input_model_part_name]

        # retrieving the input variables
        input_var_names = settings["input_variables"]
        variable_names = [ input_var_names[i].GetString() for i in range( input_var_names.size() ) ]
        input_sources_names = settings["input_sources"]
        self.input_sources = [ input_sources_names[i].GetString() for i in range( input_sources_names.size() ) ]
        self.input_variables = [ KM.KratosGlobals.GetVariable( var ) for var in variable_names ]
        if len(self.input_variables) == 0:
            raise Exception('No variables specified for input!')
        if not (len(self.input_variables) == len(self.input_sources)):
            raise Exception('The number of input variables and sources are different.')
        self.dict_input = dict(zip(self.input_variables, self.input_sources))
        # # validate types of variables
        # for var in self.input_variables:
        #     if type(var) == KM.DoubleVariable:
        #         continue
        #     elif type(var) == KM.Array1DVariable3:
        #         continue
        #     else:
        #         err_msg  = 'Type of variable "' + var.Name() + '" is not valid\n'
        #         err_msg += 'It can only be double, component or array3d!'
        #         raise Exception(err_msg)

        # retrieving the output variables
        output_var_names = settings["output_variables"]
        variable_names = [ output_var_names[i].GetString() for i in range( output_var_names.size() ) ]
        output_sources_names = settings["output_sources"]
        self.output_sources = [ output_sources_names[i].GetString() for i in range( output_sources_names.size() ) ]
        self.output_variables = [ KM.KratosGlobals.GetVariable( var ) for var in variable_names ]
        if len(self.output_variables) == 0:
            raise Exception('No variables specified for output!')
        if not (len(self.output_variables) == len(self.output_sources)):
            raise Exception('The number of output variables and sources are different.')
        self.dict_output = dict(zip(self.output_variables, self.output_sources))
        # validate types of variables
        # for var in self.output_variables:
        #     if type(var) == KM.DoubleVariable:
        #         continue
        #     elif type(var) == KM.Array1DVariable3:
        #         continue
        #     else:
        #         err_msg  = 'Type of variable "' + var.Name() + '" is not valid\n'
        #         err_msg += 'It can only be double, component or array3d!'
        #         raise Exception(err_msg)

        # Output file names handling
        if (self.write_output_file):

            model_part_name = settings["model_part_name"].GetString()
            if model_part_name == "":
                raise Exception('No "model_part_name" was specified!')
            self.model_part = model[model_part_name]

            training_input_file_name = settings["training_file_settings"]["input_file_name"].GetString(
            ) + "." + settings["training_file_settings"]["file_format"].GetString()
            training_output_file_name = settings["training_file_settings"]["output_file_name"].GetString(
            ) + "." + settings["training_file_settings"]["file_format"].GetString()

            training_input_file_handler_params = KM.Parameters()
            training_output_file_handler_params = KM.Parameters()

            if self.output_format == "ascii":

                training_input_file_handler_params.AddEmptyValue("file_name")
                training_input_file_handler_params["file_name"].SetString(training_input_file_name)
                training_input_file_handler_params.AddEmptyValue("output_path")
                training_input_file_handler_params["output_path"].SetString(settings["training_file_settings"]["output_path"].GetString())
                training_output_file_handler_params.AddEmptyValue("file_name")
                training_output_file_handler_params["file_name"].SetString(training_output_file_name)
                training_output_file_handler_params.AddEmptyValue("output_path")
                training_output_file_handler_params["output_path"].SetString(settings["training_file_settings"]["output_path"].GetString())

                # Training input
                file_header = ''
                self.training_input_file = TimeBasedAsciiFileWriterUtility(self.model_part,
                                                            training_input_file_handler_params, file_header).file

                # Training output
                self.training_output_file = TimeBasedAsciiFileWriterUtility(self.model_part,
                                                            training_output_file_handler_params, file_header).file
            elif self.output_format == "h5":
                # Training input
                training_input_file_name = settings["training_file_settings"]["output_path"].GetString() + training_input_file_name
                training_input_file_handler_params.AddEmptyValue("file_access_mode")
                training_input_file_handler_params["file_access_mode"].SetString("truncate")
                training_input_file_handler_params.AddEmptyValue("file_name")
                training_input_file_handler_params["file_name"].SetString(training_input_file_name)
                self.hdf5_file_training_input = KratosHDF5.HDF5FileSerial(training_input_file_handler_params)

                self.hdf5_input_parameters = KM.Parameters()
                self.hdf5_input_parameters.AddEmptyArray("list_of_variables")
                self.hdf5_input_parameters["list_of_variables"].SetStringArray(input_var_names.GetStringArray())
                self.hdf5_input_parameters.AddEmptyValue("prefix")
                self.hdf5_input_parameters["prefix"].SetString("/<time>/<model_part_name>/InputData")               

                # Training output
                training_output_file_name = settings["training_file_settings"]["output_path"].GetString() + training_output_file_name
                training_output_file_handler_params.AddEmptyValue("file_access_mode")
                training_output_file_handler_params["file_access_mode"].SetString("truncate")
                training_output_file_handler_params.AddEmptyValue("file_name")
                training_output_file_handler_params["file_name"].SetString(training_output_file_name)
                self.hdf5_file_training_output = KratosHDF5.HDF5FileSerial(training_output_file_handler_params)

                self.hdf5_output_parameters = KM.Parameters()
                self.hdf5_output_parameters.AddEmptyArray("list_of_variables")
                self.hdf5_output_parameters["list_of_variables"].SetStringArray(output_var_names.GetStringArray())
                self.hdf5_output_parameters.AddEmptyValue("prefix")
                self.hdf5_output_parameters["prefix"].SetString("/<time>/<model_part_name>/ResultsData")

            else:
                raise Exception("Format not supported. Supported formats are: ascii, h5.")


    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KM.TIME]
        current_step = self.model_part.ProcessInfo[KM.STEP]
        if((current_time >= self.interval[0]) and (current_time < self.interval[1])):
            input_value_list=[]
            # Perturbation of the load conditions
            # for var,dist,params in zip(self.input_variables,self.random_distribution,self.random_parameters):
            #     if self.perturbate :
            #         factor = getattr(np.random, dist)(*params)
            #         for condition in self.input_model_part.GetConditions():
            #             input_value = op.imul(condition.GetValue(var),(0.0+factor))
            #             for node in condition.GetNodes():
            #                 node.SetSolutionStepValue(var,input_value)
            #                 input_value_list.append(input_value+condition.GetValue(var))
            #     else:
            #         for condition in self.input_model_part.GetConditions():
            #             for node in condition.GetNodes():
            #                 input_value_list.append(condition.GetValue(var))
            for variable, source in self.dict_input.items():
                if variable in self.perturbate_variables:
                    index = self.perturbate_variables.index(variable)
                    dist = self.random_distribution[index]
                    rnd_parameters = self.random_parameters[index]
                    factor = getattr(np.random, dist)(*rnd_parameters)
                    # Process related variables (e.g. TIME, STEP)
                    if source == 'process':
                        raise Exception("Process variables cannot be perturbated.")
                    # Node properties (e.g. position)
                    elif source == 'node':
                        print("Warning: The node properties are being perturbated.")
                        for node in self.input_model_part.Nodes:
                            input_value = op.imul(getattr(node,variable.Name()))
                            node.SetValue(variable, input_value)
                            input_value_list.append(getattr(node,variable.Name()))
                    # Node step values (e.g. variables like displacement)
                    elif source == "solution_step":
                        for node in self.input_model_part.Nodes:
                            input_value = op.imul(node.GetSolutionStepValue(variable), (1.0+factor))
                            node.SetSolutionStepValue(variable, input_value)
                            input_value_list.append(node.GetSolutionStepValue(variable,0))
                    # Condition values
                    elif source == "condition":
                        for condition in self.input_model_part.GetConditions():
                            input_value = op.imul(condition.GetValue(variable), (1.0+factor))
                            condition.SetValue(variable, input_value)
                            input_value_list.append(condition.GetValue(variable))
                else:
                    # Process related variables (e.g. TIME, STEP)
                    if source == 'process':
                        input_value_list.append(self.input_model_part.ProcessInfo[variable])
                    # Node properties (e.g. position)
                    elif source == 'node':
                        for node in self.input_model_part.Nodes:
                            input_value_list.append(getattr(node,variable.Name()))
                    # Node step values (e.g. variables like displacement)
                    elif source == "solution_step":
                        for node in self.input_model_part.Nodes:
                            input_value_list.append(node.GetSolutionStepValue(variable,0))
                    # Condition values
                    elif source == "condition":
                        for condition in self.input_model_part.GetConditions():
                            input_value_list.append(condition.GetValue(variable))

            # Writing input file
            if (self.write_output_file):
                if self.output_format == "ascii":
                    self.training_input_file.write(' '.join(str(v) for v in input_value_list) + '\n')
                if self.output_format == "h5":
                    self.hdf5_input_parameters["prefix"].SetString("/"+str(current_step)+"/InputData") 
                    nodal_io = KratosHDF5.HDF5NodalSolutionStepDataIO(self.hdf5_input_parameters, self.hdf5_file_training_input)
                    nodal_io.WriteNodalResults(self.input_model_part,1)
                    


    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KM.TIME]
        current_step = self.model_part.ProcessInfo[KM.STEP]
        if((current_time >= self.interval[0]) and (current_time < self.interval[1])):
            output_value =[]

            for variable, source in self.dict_output.items():
                # Process related variables (e.g. TIME, STEP)
                if source == 'process':
                    output_value.append(self.model_part.ProcessInfo[variable])
                # Node properties (e.g. position)
                elif source == 'node':
                    for node in self.model_part.Nodes:
                        output_value.append(getattr(node,variable.Name()))
                # Node step values (e.g. variables like displacement)
                elif source == "solution_step":
                    for node in self.model_part.Nodes:
                        output_value.append(node.GetSolutionStepValue(variable,0))
                # Condition values
                elif source == "condition":
                    for condition in self.model_part.GetConditions():
                        output_value.append(condition.GetValue(variable,0))

            if (self.write_output_file):
                if self.output_format == "ascii":
                    self.training_output_file.write(' '.join(str(v) for v in output_value) + '\n')
                if self.output_format == "h5":
                    self.hdf5_output_parameters["prefix"].SetString("/"+str(current_step)+"/OutputData")
                    nodal_io = KratosHDF5.HDF5NodalSolutionStepDataIO(self.hdf5_output_parameters, self.hdf5_file_training_output)
                    nodal_io.WriteNodalResults(self.model_part,1)


