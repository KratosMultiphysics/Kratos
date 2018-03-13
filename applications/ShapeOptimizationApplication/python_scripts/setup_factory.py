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
                if line.strip()[0] == '#':
                    continue
                if line.strip() == '':
                    continue
                for character in line:
                    if (character == ' ' or character == '\t'):
                        python_nesting_system = '%s%s' % (python_nesting_system, character)
                    elif (character != '\n'):
                        found_nesting = True
                        break
            
            if found_nesting:
                break

        new_headers = []
        for line in lines:
            if (line[0:4] == 'from'):
                new_headers.append(line)

        new_lines = []
        is_def_main_added = False
        for line in lines:
            if (line[0:4] != 'from'):
                if not is_def_main_added:
                    new_lines.append('def main():\n')
                    is_def_main_added = True
                
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
                if lines[-1].strip() == "Simulation terminated successfully.":
                    has_executed_simulation = True
                    print(">- --- Simulation already completed.")
        
        if not has_executed_simulation:
            file_out = open("%s/%s" % (self.path, self.parameters["parameters_file"].GetString()), "w")
            file_out.write(json.dumps(self.setup_parameters, indent=4))
            file_out.close()            
            
            _temp = os.getcwd()
            os.chdir("%s/%s" % (_temp, self.path))
            sys.path.append(os.getcwd())

            libc = ctypes.CDLL(None)
            c_stdout = ctypes.c_void_p.in_dll(libc, 'stdout')
            std_io_out = io.StringIO()
            with stdout_redirector( std_io_out, c_stdout, libc ):
                self.WriteModelPart()
                self.module.main()
                print("Simulation terminated successfully.")
            std_file_out = open("log.kratos", "w")
            std_file_out.write(std_io_out.getvalue()[2:-2])
            std_file_out.close()
            
            sys.path.remove(os.getcwd())
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

@contextmanager
def stdout_redirector( stream, c_stdout, libc ):
    # The original fd stdout points to. Usually 1 on POSIX systems.
    original_stdout_fd = sys.stdout.fileno()

    def _redirect_stdout(to_fd):
        """Redirect stdout to the given file descriptor."""
        # Flush the C-level buffer stdout
        libc.fflush(c_stdout)
        # Flush and close sys.stdout - also closes the file descriptor (fd)
        sys.stdout.close()
        # Make original_stdout_fd point to the same file as to_fd
        os.dup2(to_fd, original_stdout_fd)
        # Create a new sys.stdout that points to the redirected fd
        sys.stdout = io.TextIOWrapper(os.fdopen(original_stdout_fd, 'wb'))

    # Save a copy of the original stdout fd in saved_stdout_fd
    saved_stdout_fd = os.dup(original_stdout_fd)
    try:
        # Create a temporary file and redirect stdout to it
        tfile = tempfile.TemporaryFile(mode='w+b')
        _redirect_stdout(tfile.fileno())
        # Yield to caller, then redirect stdout back to the saved fd
        yield
        _redirect_stdout(saved_stdout_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        stream.write(str(tfile.read()).replace('\\n', '\n' ))
    finally:
        tfile.close()
        os.close(saved_stdout_fd)        
