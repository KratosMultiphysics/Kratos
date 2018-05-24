from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Import time library
import time

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

class ParticleMPMGiDOutputProcess(KratosMultiphysics.Process):
    defaults = KratosMultiphysics.Parameters("""{
        "result_file_configuration": {
            "gidpost_flags": {
                "GiDPostMode": "GiD_PostBinary",
                "WriteDeformedMeshFlag": "WriteUndeformed",
                "WriteConditionsFlag": "WriteElementsOnly",
                "MultiFileFlag": "SingleFile"
            },
            "file_label": "time",
            "output_control_type": "step",
            "output_frequency": 1.0,
            "body_output": true,
            "node_output": false,
            "skin_output": false,
            "plane_output": [],
            "nodal_results": [],
            "nodal_nonhistorical_results": [],
            "nodal_flags_results": [],
            "gauss_point_results": [],
            "additional_list_files": []
        },
        "point_data_configuration": []
    }""")


    def __init__(self, model_part, file_name, param):

        KratosMultiphysics.Process.__init__(self)

        if param is None:
            param = self.defaults
        else:
            param.ValidateAndAssignDefaults(self.defaults)

        # Default
        self.param = param
        self.base_file_name = file_name
        self.model_part = model_part

        self.step_count = 0
        self.printed_step_count = 0
        self.next_output = 0.0


    # Public Functions
    def ExecuteInitialize(self):

        result_file_configuration = self.param["result_file_configuration"]
        result_file_configuration.ValidateAndAssignDefaults(self.defaults["result_file_configuration"])

        # Set up output frequency and format
        output_file_label = result_file_configuration["file_label"].GetString()
        if output_file_label == "time":
            self.output_label_is_time = True
        elif output_file_label == "step":
            self.output_label_is_time = False
        else:
            msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__,output_file_label,"file_label")
            raise Exception(msg)
        
        output_control_type = result_file_configuration["output_control_type"].GetString()
        if output_control_type == "time":
            self.output_control_is_time = True
        elif output_control_type == "step":
            self.output_control_is_time = False
        else:
            msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__,output_file_label,"file_label")
            raise Exception(msg)

        self.output_frequency = result_file_configuration["output_frequency"].GetDouble()

        # Set Variable list to print
        self.variable_name_list = result_file_configuration["gauss_point_results"]
        self.variable_list      = []
        for i in range(self.variable_name_list.size()):
            var_name = self.variable_name_list[i].GetString()
            variable = self._get_variable(var_name)
            self.variable_list.append(variable)

    def ExecuteBeforeSolutionLoop(self):
        # Initiate Output Mesh
        self.mesh_file = open(self.base_file_name + ".post.msh",'w')
        self.mesh_file.write("MESH \"")
        self.mesh_file.write("outmesh")
        self.mesh_file.write("\" dimension 3 ElemType Point Nnode 1\n")
        self.mesh_file.write("Coordinates\n")
        for mpm in self.model_part.Elements:
            coord = mpm.GetValue(KratosParticle.GAUSS_COORD)
            self.mesh_file.write("{} {} {} {}\n".format( mpm.Id, coord[0], coord[1], coord[2]))
        self.mesh_file.write("End Coordinates\n")
        self.mesh_file.write("Elements\n")
        for mpm in self.model_part.Elements:
            self.mesh_file.write("{} {}\n".format(mpm.Id, mpm.Id))
        self.mesh_file.write("End Elements\n")
        self.mesh_file.flush()

        # Initiate Output File
        self.result_file = open(self.base_file_name + ".post.res",'w')
        self.result_file.write("GiD Post Results File 1.0\n")

    def ExecuteInitializeSolutionStep(self):
        self.step_count += 1

    def ExecuteFinalizeSolutionStep(self): pass

    def ExecuteFinalize(self): pass

    def PrintOutput(self):
        # Print the output
        time = self._get_pretty_time(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        self.printed_step_count += 1
        self.model_part.ProcessInfo[KratosMultiphysics.PRINTED_STEP] = self.printed_step_count
        
        # Write results to the initiated result file
        self._write_mp_results(time)

        # Schedule next output
        if self.output_frequency > 0.0: # Note: if == 0, we'll just always print
            if self.output_control_is_time:
                while self._get_pretty_time(self.next_output) <= time:
                    self.next_output += self.output_frequency
            else:
                while self.next_output <= self.step_count:
                    self.next_output += self.output_frequency


    def IsOutputStep(self):
        if self.output_control_is_time:
            time = self._get_pretty_time(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return (time >= self._get_pretty_time(self.next_output))
        else:
            return ( self.step_count >= self.next_output )
    

    # Private Functions
    def _get_pretty_time(self,time):
        pretty_time = "{0:.12g}".format(time)
        pretty_time = float(pretty_time)
        return pretty_time

    def _get_attribute(self, my_string, function_pointer, attribute_type):
        """Return the python object named by the string argument.

        To be used with functions from KratosGlobals

        Examples:
        variable = self._get_attribute("DISPLACEMENT",
                                       KratosMultiphysics.ParticleMechanicsApplication.GetVariable,
                                       "Variable")
        """
        splitted = my_string.split(".")

        if len(splitted) == 0:
            raise Exception("Something wrong. Trying to split the string " + my_string)
        if len(splitted) > 3:
            raise Exception("Something wrong. String " + my_string + " has too many arguments")

        attribute_name = splitted[-1]

        if len(splitted) == 2 or len(splitted) == 3:
            warning_msg =  "Ignoring \"" +  my_string.rsplit(".",1)[0]
            warning_msg += "\" for " + attribute_type +" \"" + attribute_name + "\""
            KratosMultiphysics.Logger.PrintInfo("Warning in mpm gid output", warning_msg)

        return function_pointer(attribute_name) # This also checks if the application has been imported

    def _get_variable(self, my_string):
        """Return the python object of a Variable named by the string argument.

        Examples:
        recommended usage:
        variable = self._get_variable("MP_VELOCITY")
        deprecated:
        variable = self._get_variables("KratosMultiphysics.ParticleMechanicsApplication.MP_VELOCITY")
        """
        return self._get_attribute(my_string, KratosMultiphysics.KratosGlobals.GetVariable, "Variable")

    def _write_mp_results(self, step_label=None):
        clock_time = self._start_time_measure()

        for i in range(self.variable_name_list.size()):
            var_name = self.variable_name_list[i].GetString()
            variable = self.variable_list[i]

            is_scalar = self._is_scalar(variable)

            # Write in result file
            self.result_file.write("Result \"")
            self.result_file.write(var_name)
            
            if is_scalar:
                self.result_file.write('" "Kratos" {} Scalar OnNodes\n'.format(step_label))
            else:
                self.result_file.write('" "Kratos" {} Vector OnNodes\n'.format(step_label))
            
            self.result_file.write("Values\n")
            for mpm in self.model_part.Elements:
                print_variable = mpm.GetValue(variable)
                # Check whether variable is a scalar or vector
                if isinstance(print_variable, float) or isinstance(print_variable, int):
                    print_size = 1
                else:
                    print_size = print_variable.Size()

                # Write variable as formated
                if print_size == 1:
                    self.result_file.write("{} {}\n".format(mpm.Id, print_variable))
                elif print_size == 3:
                    self.result_file.write("{} {} {} {}\n".format(mpm.Id, print_variable[0], print_variable[1], print_variable[2]))
                elif print_size == 6:
                    self.result_file.write("{} {} {} {} {} {} {}\n".format(mpm.Id, print_variable[0], print_variable[1], print_variable[2], print_variable[3], print_variable[4], print_variable[5]))
                else:
                    KratosMultiphysics.Logger.PrintInfo("Warning in mpm gid output", "Printing format is not defined for variable: ", var_name, "with size: ", print_size)

            self.result_file.write("End Values\n")

        self._stop_time_measure(clock_time)

    def _start_time_measure(self):
        return time.time()

    def _stop_time_measure(self, time_ip):
        time_fp = time.time()
        KratosMultiphysics.Logger.PrintInfo("ParticleMPMGidOutputUtility", "[Spent time for output = ", time_fp - time_ip, "sec]")

    def _is_scalar(self,variable):
        is_scalar = False
        if (isinstance(variable,KratosMultiphysics.IntegerVariable) or isinstance(variable,KratosMultiphysics.DoubleVariable) or isinstance(variable,KratosMultiphysics.BoolVariable)):
            is_scalar = True
        elif (isinstance(variable,KratosMultiphysics.StringVariable)):
            raise Exception("String variable cant be printed.")
        return is_scalar