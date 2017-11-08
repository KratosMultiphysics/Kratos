# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya, https://github.com/sunethwarna
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
from contextlib import contextmanager
import json
import shutil
import os
import sys
import ctypes
import io
import tempfile

# importing the Kratos Library
from KratosMultiphysics import *
import KratosMultiphysics

class OptimizationSetup:
    def __init__ ( self,  inputModelPart, setupSettings ):
        self.setup = None
        self.model_part = inputModelPart
        self.iteration_post_fix = None
        # default settings string in json format
        default_settings = Parameters("""
        {
            "path"              : "",
            "python_module"     : "",
            "parameters_file"   : "",
            "mdpa_input_file"   : ""
        }""")
        
        self.parameters = setupSettings
        # overwrite the default settings with user-provided parameters
        # self.parameters["design_variables"].RecursivelyValidateAndAssignDefaults(default_settings)
        
        self.python_module = "%s_%s" % (self.parameters["path"].GetString(), self.parameters["python_module"].GetString())
        print(self.python_module)
        shutil.copy('%s/%s.py' %  (self.parameters["path"].GetString(), self.parameters["python_module"].GetString()), '%s.py' % self.python_module)

        self.ChangePythonModule()

        self.module = __import__(self.python_module)

        with open("%s/%s" % (self.parameters["path"].GetString(),  setupSettings["parameters_file"].GetString()), "r") as file_input:
            self.setup_parameters = json.load(file_input)
        file_input.close()

    def ChangePythonModule( self ):
        with open('%s.py' % self.python_module, 'r') as file_input:
            lines = file_input.readlines()
        file_input.close()

        python_nesting_system = ''
        found_while = False
        found_nesting = False
        for line in lines:
            if (line[0:5]=='while' or line[0:3]=='for'):
                found_while = True
                continue
            python_nesting_system = ''
            if found_while:
                for character in line:
                    if (character == ' ' or character == '\t'):
                        python_nesting_system = '%s%s' % (python_nesting_system, character)
                    elif (character != '\n'):
                        found_nesting = True
                        break
            
            if found_nesting:
                break

        new_headers = []
        new_lines = []
        found_headers = False

        for line in lines:
            if (line[0:4] == 'from'):
                new_headers.append(line)

        for line in lines:
            if (line[0:4] != 'from'):
                if not found_headers:
                    new_lines.append('def main():\n')
                    found_headers = True
                new_lines.append('%s%s' % (python_nesting_system, line))
        
        file_output = open('%s.py' % self.python_module, 'w')
        for line in new_headers:
            file_output.write(line)        
        for line in new_lines:
            file_output.write(line)
        file_output.close()                

    def AddCustomProcess( self, customProcess, overwrite = False, existingProcessPythonModuleName=""):
        if self.setup_parameters.get("list_other_processes") is not None:
            has_input_primal_solution_process = False
            if overwrite:
                for i in range(0, len(self.setup_parameters["list_other_processes"])):
                    process = self.setup_parameters["list_other_processes"][i]
                    if process["python_module"]==existingProcessPythonModuleName:
                        self.setup_parameters["list_other_processes"][i] = customProcess
                        has_input_primal_solution_process = True
                        break
            
            if not has_input_primal_solution_process:
                self.setup_parameters["list_other_processes"].append(customProcess)
        else:
            self.setup_parameters["list_other_processes"] = [customProcess]
    
    def SetCurrentPath( self, Path):
        self.path = Path
    
    def SetIterationPostFix( self, iterationPostFix):
        self.iteration_post_fix = iterationPostFix

    def Execute( self ):
        # copy the current setup to new path
        if not os.path.isdir(self.path):
            shutil.copytree(self.parameters["path"].GetString(), self.path)

        has_executed_simulation = False
        if os.path.isfile("%s/log.kratos" % (self.path)):
            print (">- --- Found existing log file.")
            with open("%s/log.kratos" % (self.path), "r") as file_input:
                lines = file_input.readlines()
            file_input.close()
            if len(lines)>0:
                if lines[-1].strip() == "KRATOS TERMINATED CORRECTLY":
                    has_executed_simulation = True
                    print(">- --- Simulation already completed.")
        
        if not has_executed_simulation:
            file_out = open("%s/%s" % (self.path, self.parameters["parameters_file"].GetString()), "w")
            file_out.write(json.dumps(self.setup_parameters, indent=4))
            file_out.close()            
            _temp = os.getcwd()
            os.chdir(self.path)
            self.WriteModelPart()
            # python_module = self.parameters["python_module"].GetString()
            # _ = os.popen("runkratos %s.py 2>&1 > log.kratos" % python_module).read()
            self.module.main()
            os.chdir(_temp)

    def WriteModelPart( self ):
        model_part_file_name = self.parameters["mdpa_input_file"].GetString()
        with open("%s.mdpa" % model_part_file_name, "r") as file_input:
            lines = file_input.readlines()
        file_input.close()

        node_lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]

        for node in self.model_part.Nodes:
            node_id = node.Id
            old_line = node_lines[node_id] # assumes consecutive node numbering starting with 1
            components = old_line.split()
            
            if int(components[0]) != node_id:
                raise RuntimeError('Error parsing file ' + model_part_file_name)
            new_line = ''
            
            if self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
                new_line = '{:5d}'.format(node_id) + ' ' \
                    + '{:19.10f}'.format(node.X0) + ' ' \
                    + '{:19.10f}'.format(node.Y0) + ' ' \
                    + '{:19.10f}'.format(node.Z0) + '\n'
            else:
                new_line = '{:5d}'.format(node_id) + ' ' \
                            + '{:19.10f}'.format(node.X0) + ' ' \
                            + '{:19.10f}'.format(node.Y0) + ' ' \
                            + '0.0000000000\n'
            lines[lines.index(old_line)] = new_line

        with open(model_part_file_name + '.mdpa', 'w') as model_part_file:
            model_part_file.writelines(lines)
        model_part_file.close()
