# This script couples DEM and Nvidia Flex...

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import sys

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.NvidiaFlexApplication import *

import main_script

class SolutionDEM(main_script.Solution):

    def __init__(self):
        super(SolutionDEM, self).__init__()
        self.changing_gravity_option = True
        self.velocity_threshold_for_gravity_change = 0.02
        self.time_at_last_gravity_change = 0.0
        self.min_time_between_gravity_changes = 0.25
        self.max_time_between_gravity_changes = 3.0
        self.gravities_filename = "TimeAngle.csv"
        self.list_of_gravities = []
        self.gravity_iterator_position = 0
        self.stop_signal = False

    def Initialize(self):
        super(SolutionDEM, self).Initialize()
        self.creator_destructor.DestroyParticlesOutsideBoundingBox(self.spheres_model_part) #TODO: why this?
        if self.changing_gravity_option:
            self._ReadFileWithGravities()
            if not self.list_of_gravities:
                print("No gravities found in the list of file:" + self.gravities_filename)
            self.spheres_model_part.ProcessInfo[GRAVITY_X]= self.list_of_gravities[0][0]
            self.spheres_model_part.ProcessInfo[GRAVITY_Y]= self.list_of_gravities[0][1]
            self.spheres_model_part.ProcessInfo[GRAVITY_Z]= self.list_of_gravities[0][2]

    def SolverSolve(self):
        super(SolutionDEM, self).SolverSolve()
        if self.changing_gravity_option:
            if NvidiaFlexPreUtilities().CheckIfItsTimeToChangeGravity(self.spheres_model_part, self.time_at_last_gravity_change, self.velocity_threshold_for_gravity_change, self.min_time_between_gravity_changes, self.max_time_between_gravity_changes):
                #TODO: utility of max_time_between_gravity_changes??
                self._ChangeGravity()

    def _ReadFileWithGravities(self):
        print("Reading gravity values from columns 3 to 5, starting from 0, from file " + self.gravities_filename + "...")
        with open(self.gravities_filename) as gravities_file:
            for line in gravities_file:
                array_of_tokens = line.split(",")
                gravities_array = [9.81*float(array_of_tokens[3]), 9.81*float(array_of_tokens[4]), 9.81*float(array_of_tokens[5])]
                self.list_of_gravities.append(gravities_array)

    def _ChangeGravity(self):
        print("Changing gravity direction to gravity number " + str(self.gravity_iterator_position + 1) +" ...")
        if self.gravity_iterator_position >= len(self.list_of_gravities) - 1:
            print("No more gravities available. Exiting.")
            self.stop_signal = True

        self.gravity_iterator_position += 1
        new_gravity = self.list_of_gravities[self.gravity_iterator_position]
        self.spheres_model_part.ProcessInfo[GRAVITY_X]= new_gravity[0]
        self.spheres_model_part.ProcessInfo[GRAVITY_Y]= new_gravity[1]
        self.spheres_model_part.ProcessInfo[GRAVITY_Z]= new_gravity[2]

    def BreakSolutionStepsLoop(self):
        if self.stop_signal:
            return True


class SolutionFlex(SolutionDEM):

    def __init__(self):
        super(SolutionFlex, self).__init__()
        self.number_of_steps_until_flex_update = 10

    def Run(self):
        self.nvidia_flex_wrapper = FlexWrapper(self.spheres_model_part, self.rigid_face_model_part, self.creator_destructor)
        super(SolutionFlex, self).Run()

    def SolverSolve(self):
        if self.step < 2:
            self._CheckNvidiaParameters()
            self.nvidia_flex_wrapper.UpdateFlex(True, True)
        else:
            if not self.step % self.number_of_steps_until_flex_update:
                self._CheckNvidiaParameters()
                self.nvidia_flex_wrapper.UpdateFlex(True, False)

        self.nvidia_flex_wrapper.SolveTimeSteps(self.dt, 1) #DO NOT CHANGE THIS 1, OR INSTABILITIES MAY APPEAR
        self.nvidia_flex_wrapper.TransferDataFromFlexToKratos()
        if self.changing_gravity_option:
            if NvidiaFlexPreUtilities().CheckIfItsTimeToChangeGravity(self.spheres_model_part, self.time_at_last_gravity_change, self.velocity_threshold_for_gravity_change, self.min_time_between_gravity_changes, self.max_time_between_gravity_changes):
                self._ChangeGravity()

    def _CheckNvidiaParameters(self):
        min_time_step = 1e-3
        #TODO: Here we should find a tendency of stability based on the radius of the particles. This will require some additional research
        if False: #self.DEM_parameters["MaxTimeStep"].GetDouble()< min_time_step:
            Logger.PrintWarning("NVIDIA APP", "Too small time step. Please use values over", str(min_time_step), ". Exiting.")
            sys.exit()

if __name__=="__main__":

    if len(sys.argv) < 2:
        print("Too few arguments. Specify 'DEM' or 'Flex', like 'python3 " + str(sys.argv[0]) + " Flex'. Exiting.")
        sys.exit()

    if len(sys.argv) > 2:
        print("Too many arguments. Specify 'DEM' or 'Flex', like 'python3 " + str(sys.argv[0]) + " Flex'. Exiting.")
        sys.exit()

    if sys.argv[1] == "DEM":
        print("Running Case with DEM...")
        SolutionDEM().Run()
    elif sys.argv[1] == "Flex":
        print("Running Case with Flex...")
        SolutionFlex().Run()
    else:
        print("Argument not understood. Specify 'DEM' or 'Flex', like 'python3 " + str(sys.argv[0]) + " Flex'. Exiting.")
