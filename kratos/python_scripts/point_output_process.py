from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# other imports
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PointOutputProcess(Model, settings["Parameters"])

class PointOutputProcess(KratosMultiphysics.Process):
    """This process writes results from a geometrical position (point) in the model to a file
    It first searches the entity containing the requested output location and then interpolates
    the requested variable(s)
    The output can be requested for elements, conditions and nodes. For nodes no geometrical
    interpolation is performed, the exact coordinates have to be specified.

    This process works in MPI as well as with restarts

    It can serve as a basis for other processes (e.g. MultiplePointsOutputProcess)
    Furthermore it can be used for testing in MPI where the node numbers can change
    """
    def __init__(self, model, params):

        default_settings = KratosMultiphysics.Parameters('''{
            "position"         : [],
            "model_part_name"  : "",
            "output_file_name" : "",
            "output_variables" : [],
            "entity_type"      : "element"
        }''')

        self.model = model

        self.params = params
        self.params.ValidateAndAssignDefaults(default_settings)

        # These quantites are lists such that they can be looped
        # => needed for mpi in case the point is in a different partition
        self.output_file = []
        self.entity = []
        self.area_coordinates = []
        self.output_variables = []

    def ExecuteInitialize(self):
        # getting the ModelPart from the Model
        model_part_name = self.params["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        self.model_part = self.model[model_part_name]

        # retrieving the position of the entity
        point_position = self.params["position"].GetVector()
        if point_position.Size() != 3:
            raise Exception('The position has to be provided with 3 coordinates!')
        point = KratosMultiphysics.Point(point_position[0],
                                         point_position[1],
                                         point_position[2])

        # retrieving the output variables
        output_var_names = self.params["output_variables"]
        variable_names = [ output_var_names[i].GetString() for i in range( output_var_names.size() ) ]
        output_vars = [ KratosMultiphysics.KratosGlobals.GetVariable( var ) for var in variable_names ]
        if len(output_vars) == 0:
            raise Exception('No variables specified for output!')
        self.output_variables.append(output_vars)
        # validate types of variables
        for var in self.output_variables[0]:
            self.__CheckVariableIsSolutionStepVariable(var)
            if type(var) == KratosMultiphysics.DoubleVariable:
                continue
            elif type(var) == KratosMultiphysics.Array1DVariable3:
                continue
            elif type(var) == KratosMultiphysics.Array1DComponentVariable:
                continue
            else:
                err_msg  = 'Type of variable "' + var.Name() + '" is not valid\n'
                err_msg += 'It can only be double, component or array3d!'
                raise Exception(err_msg)

        # retrieving the entity type
        entity_type = self.params["entity_type"].GetString()

        if entity_type == "node":
            tol = 1e-12
            found_id = KratosMultiphysics.PointLocator(self.model_part).FindNode(point, tol)
            if found_id > -1:
                self.entity.append(self.model_part.Nodes[found_id]) # note that this is a find!
                self.area_coordinates.append("dummy") # needed for looping later
        elif entity_type == "element":
            self.sf_values = KratosMultiphysics.Vector()
            found_id = KratosMultiphysics.PointLocator(self.model_part).FindElement(point, self.sf_values)
            if found_id > -1:
                self.entity.append(self.model_part.Elements[found_id]) # note that this is a find!
                self.area_coordinates.append(self.sf_values)
        elif entity_type == "condition":
            self.sf_values = KratosMultiphysics.Vector()
            found_id = KratosMultiphysics.PointLocator(self.model_part).FindCondition(point, self.sf_values)
            if found_id > -1:
                self.entity.append(self.model_part.Conditions[found_id]) # note that this is a find!
                self.area_coordinates.append(self.sf_values)
        else:
            err_msg  = 'Invalid "entity_type" specified, it can only be:\n'
            err_msg += '"node", "element", "condition"'
            raise Exception(err_msg)

        # Check if a point was found, and initalize output
        # NOTE: If the search was not successful (i.e. found_id = -1), we fail silently and
        # do nothing. This is BY DESIGN, as we are supposed to work on MPI too, and the point
        # in question might lie on a different partition.
        if found_id > -1:
            # setting up the output_file
            output_file_name = self.params["output_file_name"].GetString()
            if output_file_name == "":
                raise Exception('No "output_file_name" was specified!')
            if not output_file_name.endswith('.dat'):
                output_file_name += ".dat"

            if self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                # if the simulation is restarted we search for a  file
                # existing from the previous run to append the data
                restart_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
                existing_file_is_valid, out_file = AddToExistingOutputFile(output_file_name,
                                                                           restart_time)

                if existing_file_is_valid:
                    self.output_file.append(out_file)

                # if no valid file can be found we create a new one
                # and issue a warning
                else:
                    warn_msg  = "No data file was found after restarting,\n"
                    warn_msg += "writing to a new file"
                    KratosMultiphysics.Logger.PrintWarning("PointOutputProcess", warn_msg)

                    self.output_file.append(InitializeOutputFile(output_file_name,
                                                                 entity_type,
                                                                 found_id,
                                                                 point,
                                                                 self.output_variables[0]))

            else: # no restart, regular simulation
                self.output_file.append(InitializeOutputFile(output_file_name,
                                                             entity_type,
                                                             found_id,
                                                             point,
                                                             self.output_variables[0]))

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        # zip works with the shortes list, which is what we want here
        # i.e. if no entity was found then also no output_file will be
        # initialized which means that the loop body will never be executed
        for var_list,ent,coord,f in zip(self.output_variables, self.entity, self.area_coordinates, self.output_file):
            out = "{0:.12g}".format(time) # making the format more readable
            for var in var_list:
                value = Interpolate(var, ent, coord)

                if IsArrayVariable(var):
                    out += " " + " ".join( "{0:.12g}".format(v) for v in value )
                else:
                    out += " " + "{0:.12g}".format(value)

            out += "\n"
            f.write(out)

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        self.__CloseOutputFiles()

    def __del__(self):
         # in case "ExecuteFinalize" is not called
         # this can happen if a simulation is forcefully stopped on a cluster
        self.__CloseOutputFiles()

    def __CloseOutputFiles(self):
        '''Close output files.'''
        for f in self.output_file:
            f.close()

    def __CheckVariableIsSolutionStepVariable(self, var):
        # if the requested Variable is a component we check the source Variable
        if type(var) == KratosMultiphysics.Array1DComponentVariable:
            var = var.GetSourceVariable()

        if not self.model_part.HasNodalSolutionStepVariable(var):
            err_msg  = 'ModelPart "' + self.model_part.Name + '" does not have'
            err_msg += ' "' + var.Name() + '" as SolutionStepVariable!'
            raise Exception(err_msg)


def InitializeOutputFile(output_file_name, entity_type, entity_id, point, output_variables):
    output_file = open(output_file_name,"w")
    out  = '# Results for "' + entity_type + '" '
    out += 'with ID: ' + entity_id + " at position: '
    out += 'x: ' + "{0:.12g}".format(point.X) + '; '
    out += 'y: ' + "{0:.12g}".format(point.Y) + '; '
    out += 'z: ' + "{0:.12g}".format(point.Z) + '\n'
    out += '# time'
    for var in output_variables:
        # if this is a Variable< array_1d< double,3 > >
        if IsArrayVariable(var):
            out += " {0}_X {0}_Y {0}_Z".format(var.Name())
        else:
            out += " " + var.Name()

    out += "\n"
    output_file.write(out)

    return output_file

def AddToExistingOutputFile(output_file_name, restart_time):
    if not os.path.isfile(output_file_name):
        return False, None

    try: # We try to open the file and transfer the info
        with open(output_file_name,'r') as out_file:
            lines_existing_file = out_file.readlines()

        output_file = open(output_file_name,"w") # this overwrites the old file

        # search for time, return false if it was not found
        # copy corresponding lines to new file and open it
        is_found = False

        for line in lines_existing_file:
            output_file.write(line)
            if line.startswith(str(restart_time)):
                is_found = True
                break

        if not(is_found):
            warn_msg  = "No line was found in " + output_file_name + " after restarting containing indicated restart time, \n"
            warn_msg += "appending results after restart from time " + str(restart_time) + " not possible"
            KratosMultiphysics.Logger.PrintWarning("PointOutputProcess", warn_msg)

        return True, output_file
    except:
        return False, None


def Interpolate(variable, entity, sf_values):
    if type(entity) == KratosMultiphysics.Node:
        return entity.GetSolutionStepValue(variable)
    else: # entity is element or condition
        nodes = entity.GetNodes()
        # Initializing 'value' like this, i don't need to know its type
        # => this way it works both for scalar and array3 variables
        value = nodes[0].GetSolutionStepValue(variable) * sf_values[0]
        for n,c in zip(nodes[1:], sf_values[1:]):
            value = value + c * n.GetSolutionStepValue(variable)

        return value

def IsArrayVariable(var):
    return type(var) == KratosMultiphysics.Array1DVariable3
