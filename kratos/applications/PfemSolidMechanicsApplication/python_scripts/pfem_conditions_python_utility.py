from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
CheckForPreviousImport()

import math

class ConditionsUtility:
    #

    def __init__(self, model_part, domain_size, incr_disp, incr_load, rotation_dofs):

        self.model_part = model_part
        self.domain_size = domain_size

        # set time evolution
        self.incr_disp = False
        if(incr_disp == "True"):
            self.incr_disp = True

        self.incr_load = False
        if(incr_load == "True"):
            self.incr_load = True

        self.rotation_dofs = rotation_dofs;

    #
    def Initialize(self, time_step):
        self.SetIncrementalDisp(time_step)
        if(self.rotation_dofs):
            self.SetIncrementalRotation(time_step)

    #
    def SetWeight(self):
        for node in self.model_part.Nodes:
            gravetat = node.GetSolutionStepValue(VOLUME_ACCELERATION);
            gravetat[1] = -10;
            node.SetSolutionStepValue(VOLUME_ACCELERATION, gravetat);

        s1 = -10.0;
        s2 = -20.0;
        setWeight = SetMechanicalInitialStateProcess(self.model_part, True, s1, s2)
        setWeight.ExecuteInitialize() ;

    #
    def SetConstantWeight(self, s1, s2):
         for node in self.model_part.Nodes:
            gravetat = node.GetSolutionStepValue(VOLUME_ACCELERATION);
            gravetat[1] = -0;
         #setWeight = SetMechanicalInitialState()
         setWeight = SetMechanicalInitialStateProcess(self.model_part, False, s1, s2)
         setWeight.ExecuteInitialize();

    #
    def SetIncrementalDisp(self, time_step):

        for node in self.model_part.Nodes:
            ImposedDisp = node.GetSolutionStepValue(IMPOSED_DISPLACEMENT)
            Displacement = node.GetSolutionStepValue(DISPLACEMENT)
            Velocity = node.GetSolutionStepValue(VELOCITY)
            # WaterPressure = node.GetSolutionStepValue(WATER_PRESSURE)
            # For displacement imposition:
            if(node.IsFixed(DISPLACEMENT_X) == 1):
                ImposedDisp[0] = Displacement[0]
                Displacement[0] = 0
            if(node.IsFixed(DISPLACEMENT_Y) == 1):
                ImposedDisp[1] = Displacement[1];
                Displacement[1] = 0;
            if(node.IsFixed(DISPLACEMENT_Z) == 1):
                ImposedDisp[2] = Displacement[2];
                Displacement[2] = 0;

            ## THERE EXIST THE POSSIBILITY THAT WATER PRESSURE IS NOT DEFINED...)
            if(node.IsFixed(WATER_PRESSURE) == 1):
                WaterPressure = node.GetSolutionStepValue(WATER_PRESSURE);
                ImposedWaterPressure = WaterPressure;
                node.SetSolutionStepValue( WATER_PRESSURE, WaterPressure);
                node.SetSolutionStepValue( IMPOSED_WATER_PRESSURE, ImposedWaterPressure);

            # For velocity imposition instead of displacement
            #if(node.IsFixed(VELOCITY_X) == 1):
            #    ImposedDisp[0] = Velocity[0] * time_step;
            #    Velocity[0] = 0;
            #if(node.IsFixed(VELOCITY_Y) == 1):
            #    ImposedDisp[1] = Velocity[1] * time_step;
            #    Velocity[1] = 0;
            #if(node.IsFixed(VELOCITY_Z) == 1):
            #    ImposedDisp[2] = Velocity[2] * time_step;
            #    Velocity[2] = 0;

            # print " ImposedDisp  =", ImposedDisp
            # print " Displacement =", Displacement
            # print " Velocity     =", Velocity

            node.SetSolutionStepValue(IMPOSED_DISPLACEMENT, ImposedDisp);

            # set to buffer variables to zero
            node.SetSolutionStepValue(DISPLACEMENT, Displacement);
            node.SetSolutionStepValue(VELOCITY, Velocity);


    #
    def SetIncrementalRotation(self, time_step):


        for node in self.model_part.Nodes:
            ImposedRotation = node.GetSolutionStepValue(IMPOSED_ROTATION)
            Rotation = node.GetSolutionStepValue(ROTATION)

            # For displacement imposition:
            if(node.IsFixed(ROTATION_X) == 1):
                ImposedRotation[0] = Rotation[0]
                Rotation[0] = 0
            if(node.IsFixed(ROTATION_Y) == 1):
                ImposedRotation[1] = Rotation[1];
                Rotation[1] = 0;
            if(node.IsFixed(ROTATION_Z) == 1):
                ImposedRotation[2] = Rotation[2];
                Rotation[2] = 0;

            node.SetSolutionStepValue(IMPOSED_ROTATION, ImposedRotation)

            # set to buffer variables to zero
            node.SetSolutionStepValue(ROTATION, Rotation)


    def CorrectBoundaryConditions(self,incr_steps,time_step):


        for node in self.model_part.Nodes:
            CF = node.GetSolutionStepValue(CONTACT_FORCE);
            if ( abs(CF[0]) + abs(CF[1]) ):
                LineLoad = node.GetSolutionStepValue(LINE_LOAD);
                LineLoad = 0.0*LineLoad;
                node.SetSolutionStepValue(LINE_LOAD, LineLoad);
                if (node.HasDofFor(WATER_PRESSURE)):
                    if (node.IsFixed(WATER_PRESSURE)):
                        print( " CORRECTING BOUNDARY CONDITION ")
                        node.Free(WATER_PRESSURE);




    #
    def SetIncrementalLoad(self, incr_steps, time_step):
        if(self.incr_load):
            for node in self.model_part.Nodes:

                # line load conditions
                force = node.GetSolutionStepValue(LINE_LOAD);
                for dim in range(0,len(force)):
                    force[dim] = force[dim] / (time_step * (incr_steps))
                force = force * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(LINE_LOAD, force);

                # surface load conditions
                force = node.GetSolutionStepValue(SURFACE_LOAD);
                for dim in range(0,len(force)):
                    force[dim] = force[dim] / (time_step * (incr_steps))
                force = force * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(SURFACE_LOAD, force);

                # line or surface pressure conditions
                pressure = node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                pressure = pressure / (time_step * (incr_steps))
                pressure = pressure * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE, pressure);

                pressure = node.GetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
                pressure = pressure / (time_step * (incr_steps))
                pressure = pressure * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(NEGATIVE_FACE_PRESSURE, pressure);

                # point load conditions
                force = node.GetSolutionStepValue(POINT_LOAD);
                for dim in range(0,len(force)):
                    force[dim] = force[dim] / (time_step * (incr_steps))
                force = force * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(POINT_LOAD, force);

                # point moment conditions
                moment = node.GetSolutionStepValue(POINT_MOMENT);
                for comp in [0,1,2]:
                    moment[comp] = moment[comp] / (time_step * (incr_steps))
                    moment[comp] = moment[comp] * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(POINT_MOMENT, moment);

    #
    def RestartImposedDisp(self):

        if(self.incr_disp == False):
            for node in self.model_part.Nodes:
                ImposedDisp = node.GetSolutionStepValue(IMPOSED_DISPLACEMENT)

                # For displacement imposition:
                if(node.IsFixed(DISPLACEMENT_X) == 1):
                    ImposedDisp[0] = 0
                if(node.IsFixed(DISPLACEMENT_Y) == 1):
                    ImposedDisp[1] = 0;
                if(node.IsFixed(DISPLACEMENT_Z) == 1):
                    ImposedDisp[2] = 0;

                node.SetSolutionStepValue(IMPOSED_DISPLACEMENT, ImposedDisp);

            if(self.rotation_dofs == True):
                self.RestartImposedRotation();


    #
    def RestartImposedRotation(self):

        if(self.incr_disp == False):
            for node in self.model_part.Nodes:
                ImposedRotation = node.GetSolutionStepValue(IMPOSED_ROTATION)

                # For displacement imposition:
                if(node.IsFixed(ROTATION_X) == 1):
                    ImposedRotation[0] = 0
                if(node.IsFixed(ROTATION_Y) == 1):
                    ImposedRotation[1] = 0;
                if(node.IsFixed(ROTATION_Z) == 1):
                    ImposedRotation[2] = 0;

                node.SetSolutionStepValue(IMPOSED_ROTATION, ImposedRotation);

    #
