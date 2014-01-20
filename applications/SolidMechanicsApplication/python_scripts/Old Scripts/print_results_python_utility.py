from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()


class PrintResultsUtility:
    #

    def __init__(self, model_part, problem_type, solver_type, problem_name, output_mode, files_mode):

        self.model_part = model_part

        # set problem type (Mechanical,Thermal,..)
        self.problem_type = problem_type

        # set solver type (Static,Dynamic,...)
        self.solver_type = solver_type

        # set problem name
        self.problem_name = problem_name

        # set gid options:
        self.gid_output_mode = output_mode         # ascii or binary
        self.gid_files_mode = files_mode          # single or multiple files

        # set gid print flags
        self.write_deformed = WriteDeformedMeshFlag.WriteUndeformed
        self.write_conditions = WriteConditionsFlag.WriteElementsOnly

        # set gid print options:
        self.write_particles = False

        # set print variables
        self.print_id = 0
        self.print_step = 0
        self.print_frequency = 1

    #
    def Initialize(self, initial_time, steps_number, initial_print_step, initial_print_id):

        # set initial time
        self.initial_time = initial_time

        # set initial print step
        self.print_step = initial_print_step

        # set initial print id
        self.print_id = initial_print_id

        # set step number
        self.total_steps_number = steps_number - 1

        # set gid input-output class
        self.gid_io = GidIO(self.problem_name, self.gid_output_mode, self.gid_files_mode, self.write_deformed, self.write_conditions)

    #
    def SetPrintOptions(self, write_particles, write_deformed, write_conditions, write_frequency):

        if(write_particles == "True"):
            self.write_particles = True
        else:
            self.write_particles = False

        if(write_deformed == "Undeformed"):
            self.write_deformed = WriteDeformedMeshFlag.WriteUndeformed
        else:
            self.write_deformed = WriteDeformedMeshFlag.WriteDeformed

        if(write_conditions == "True"):
            self.write_conditions = WriteConditionsFlag.WriteElementsOnly
        else:
            self.write_conditions = WriteConditionsFlag.WriteConditions

        # set print frequency
        self.print_frequency = write_frequency

    #
    def PrintInitialMesh(self):
        initial_id = 1
        self.gid_io.InitializeMesh(initial_id)
        self.gid_io.WriteMesh((self.model_part).GetMesh())
        self.gid_io.FinalizeMesh()

    #
    def PrintResults(self, current_time, current_step, list_files):

        step_printed = False

        if(current_step == self.total_steps_number or current_step == self.print_step):
            print("WRITING RESULTS: [STEP: ", self.print_id, "] [TIME: ", current_time - self.initial_time, "]")
            # initialize results
            self.gid_io.InitializeResults(self.print_id, (self.model_part).GetMesh())

            print(self.gid_files_mode)
            write_mesh = True
            if(self.gid_files_mode == "SingleFile"):
                write_mesh = False

            if(write_mesh):
                # initialize mesh
                self.gid_io.InitializeMesh(self.print_id)

                # write mesh nodes
                if(self.write_particles):
                    self.gid_io.WriteNodeMesh((self.model_part).GetMesh())

                # write mesh elements and conditions
                self.gid_io.WriteMesh((self.model_part).GetMesh())

                # finalize mesh
                self.gid_io.FinalizeMesh()

            # set total time
            total_time = current_time - self.initial_time

            print("Writing Results", self.problem_type)
            # print variables
            if(self.problem_type == "Mechanical" or self.problem_type == "ThermoMechanical"):
                self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part.Nodes, total_time, 0)
                # self.gid_io.WriteNodalResults(PRESSURE,self.model_part.Nodes,total_time,0)
                self.gid_io.WriteNodalResults(REACTION, self.model_part.Nodes, total_time, 0)
                print("line125")
                self.gid_io.PrintOnGaussPoints(CAUCHY_STRESS_TENSOR, self.model_part, total_time)
                print("line126")
                self.gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, total_time)
                print("line128")
                self.gid_io.PrintOnGaussPoints(VON_MISES_STRESS, self.model_part, total_time)
                # self.gid_io.PrintOnGaussPoints(PLASTIC_STRAIN,self.model_part,total_time)
                # self.gid_io.PrintOnGaussPoints(DELTA_PLASTIC_STRAIN,self.model_part,total_time)

                # self.gid_io.WriteNodalResults(NORMAL,self.model_part.Nodes,total_time,0)
                # self.gid_io.WriteNodalResults(MEAN_ERROR,self.model_part.Nodes,total_time,0)
                # self.gid_io.WriteNodalResults(OFFSET,self.model_part.Nodes,total_time,0)

                # self.gid_io.WriteNodalResults(FORCE_INTERNAL,self.model_part.Nodes,total_time,0)
                self.gid_io.WriteNodalResults(FORCE_EXTERNAL, self.model_part.Nodes, total_time, 0)
                # self.gid_io.WriteNodalResults(FORCE_CONTACT_NORMAL,self.model_part.Nodes,total_time,0)
                # self.gid_io.WriteNodalResults(FORCE_CONTACT_TANGENT,self.model_part.Nodes,total_time,0)

            if(self.solver_type == "DynamicSolver"):
                print("line143")
                self.gid_io.WriteNodalResults(VELOCITY, self.model_part.Nodes, total_time, 0)
                self.gid_io.WriteNodalResults(ACCELERATION, self.model_part.Nodes, total_time, 0)

            # flush gid writing
            self.gid_io.Flush()

            # finalize results
            self.gid_io.FinalizeResults()

            # print list files
            list_files.PrintListFiles(current_step)

            # update time variables
            if(self.print_step == 0):
                self.print_step = self.print_frequency
            else:
                self.print_step = self.print_step + self.print_frequency

            self.print_id = self.print_id + 1
            self.model_part.ProcessInfo[WRITE_ID] = self.print_id;
            step_printed = True;

            print(" -Run_GID_for_viewing_the_results_of_the_analysis-")

        return step_printed


    #
