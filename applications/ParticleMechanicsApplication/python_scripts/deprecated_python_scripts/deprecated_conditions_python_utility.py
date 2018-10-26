from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()


class ConditionsUtility:
    #
    
    def __init__(self, model_part, domain_size, incr_disp, incr_load, incr_var):

    
        self.model_part = model_part
        self.domain_size = domain_size
        
        # set time evolution
        self.incr_disp = False
        if(incr_disp == "True"):
            self.incr_disp = True

        self.incr_load = False
        if(incr_load == "True"):
            self.incr_load = True

    #
    def Initialize(self, time_step):
        self.SetIncrementalDisp(time_step)

    #
    def SetIncrementalDisp(self, time_step):
        print("Currently the SetIncremenetalDisp function is not active.")
        
    #     for node in self.model_part.Nodes:
            
    #         ImposedDisp = node.GetSolutionStepValue(IMPOSED_DISPLACEMENT)
            
    #         Displacement = node.GetSolutionStepValue(DISPLACEMENT)
    #         Velocity = node.GetSolutionStepValue(VELOCITY)
            
    #         # For displacement imposition:
    #         if(node.IsFixed(DISPLACEMENT_X) == 1):
    #             ImposedDisp[0] = Displacement[0]
    #             Displacement[0] = 0
    #         if(node.IsFixed(DISPLACEMENT_Y) == 1):
    #             ImposedDisp[1] = Displacement[1]
    #             Displacement[1] = 0
    #         if(node.IsFixed(DISPLACEMENT_Z) == 1):
    #             ImposedDisp[2] = Displacement[2]
    #             Displacement[2] = 0

    #         # For velocity imposition instead of displacement
    #         if(node.IsFixed(VELOCITY_X) == 1):
    #             ImposedDisp[0] = Velocity[0] * time_step
    #             Velocity[0] = 0
    #         if(node.IsFixed(VELOCITY_Y) == 1):
    #             ImposedDisp[1] = Velocity[1] * time_step
    #             Velocity[1] = 0
    #         if(node.IsFixed(VELOCITY_Z) == 1):
    #             ImposedDisp[2] = Velocity[2] * time_step
    #             Velocity[2] = 0
                
    #         # print " ImposedDisp  =", ImposedDisp
    #         # print " Displacement =", Displacement
    #         # print " Velocity     =", Velocity

    #         node.SetSolutionStepValue(IMPOSED_DISPLACEMENT, ImposedDisp)

    #         # set to buffer variables to zero
    #         node.SetSolutionStepValue(DISPLACEMENT, Displacement)
    #         node.SetSolutionStepValue(VELOCITY, Velocity)



    #
    def SetIncrementalLoad(self, incr_steps, time_step):
        if(self.incr_load):
            for node in self.model_part.Nodes:

                # line load conditions
                force = node.GetSolutionStepValue(LINE_LOAD)
                for dim in range(0,len(force)):
                    force[dim] = force[dim] / (time_step * (incr_steps))
                force = force * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(LINE_LOAD, force)

                # surface load conditions
                force = node.GetSolutionStepValue(SURFACE_LOAD)
                for dim in range(0,len(force)):
                    force[dim] = force[dim] / (time_step * (incr_steps))
                force = force * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(SURFACE_LOAD, force)

                # line or surface pressure conditions
                pressure = node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE)
                pressure = pressure / (time_step * (incr_steps))
                pressure = pressure * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE, pressure)

                pressure = node.GetSolutionStepValue(NEGATIVE_FACE_PRESSURE)
                pressure = pressure / (time_step * (incr_steps))
                pressure = pressure * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(NEGATIVE_FACE_PRESSURE, pressure)

                # point load conditions
                force = node.GetSolutionStepValue(POINT_LOAD)
                for dim in range(0,len(force)):
                    force[dim] = force[dim] / (time_step * (incr_steps))
                force = force * time_step * (incr_steps + 1)
                node.SetSolutionStepValue(POINT_LOAD, force)
                

    #
    #def RestartImposedDisp(self):  #msi: there should be an Else where the imposed displacement comming from a imposed velocity is recalculated in case of modifying timestep
        #if(self.incr_disp == False):
            #for node in self.model_part.Nodes:
                #ImposedDisp = node.GetSolutionStepValue(IMPOSED_DISPLACEMENT)
 
                ## For displacement imposition:
                #if(node.IsFixed(DISPLACEMENT_X) == 1):
                    #ImposedDisp[0] = 0
                #if(node.IsFixed(DISPLACEMENT_Y) == 1):
                    #ImposedDisp[1] = 0
                #if(node.IsFixed(DISPLACEMENT_Z) == 1):
                    #ImposedDisp[2] = 0

                #node.SetSolutionStepValue(IMPOSED_DISPLACEMENT, ImposedDisp)
        
