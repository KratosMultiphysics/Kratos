from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPFEMSolid
import KratosMultiphysics.ContactMechanicsApplication as KratosContact


import os.path


def Factory( custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SoilSamplerUtility(Model, custom_settings["Parameters"])


class SoilSamplerUtility(KratosMultiphysics.Process):

    #
    def __init__(self, model_part, custom_settings):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
         "inner_separation_point": [0.03,0,0],
         "outer_separation_point": [0.0375,0,0],
         "bottom_point"          : [0.03375, -0.00375, 0.0], 
         "piston_position": [100000.0, 100000.0, 100000.0],
         "velocity": 0.0375
        }
        """)

        #LMV: OBS: be careful with the tip and so on

        settings = custom_settings
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model_part['Main_Domain']

        self.velocity = settings["velocity"].GetDouble()

        self.InnerSeparationPoint = KratosMultiphysics.Array3()
        self.OuterSeparationPoint = KratosMultiphysics.Array3()
        self.BottomPoint          = KratosMultiphysics.Array3()
        self.PistonPosition       = KratosMultiphysics.Array3()

        for ii in range(0,3):
            self.InnerSeparationPoint[ii] = settings["inner_separation_point"][ii].GetDouble()
            self.OuterSeparationPoint[ii] = settings["outer_separation_point"][ii].GetDouble()
            self.BottomPoint[ii]          = settings["bottom_point"][ii].GetDouble()
            self.PistonPosition[ii]       = settings["piston_position"][ii].GetDouble()



        self.InnerForce  = KratosMultiphysics.Array3()
        self.OuterForce  = KratosMultiphysics.Array3()
        self.TipForce    = KratosMultiphysics.Array3()
        self.PistonForce = KratosMultiphysics.Array3()

        self.InnerLength = 0.0
        self.OuterLength = 0.0

        self.Filling = 0.0

        # initialize figure path
        problem_path = os.getcwd()
        self.figure_path = os.path.join(problem_path, "TubeSampler.csv")


    def ExecuteBeforeOutputStep(self):

        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]


        for node in self.model_part.GetNodes(0):
            delta_disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            delta_disp = delta_disp -  node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1)
            for ii in range(0,2):
                delta_disp[ii] = delta_disp[ii] / delta_time
            if (node.SolutionStepsDataHas(KratosMultiphysics.VELOCITY)):
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, delta_disp)
            
    def ExecuteAfterOutputStep(self):

        for node in self.model_part.GetNodes(0):
            if ( node.SolutionStepsDataHas( KratosMultiphysics.CONTACT_FORCE)):
                CF = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)
                CF = 0.0*CF
                node.SetSolutionStepValue( KratosMultiphysics.CONTACT_FORCE, CF)
            if ( node.SolutionStepsDataHas( KratosMultiphysics.CONTACT_NORMAL)):
                CF = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_NORMAL)
                CF = 0.0*CF
                node.SetSolutionStepValue( KratosMultiphysics.CONTACT_NORMAL, CF)

            if ( node.SolutionStepsDataHas( KratosContact.CONTACT_STRESS)):
                CF = node.GetSolutionStepValue( KratosContact.CONTACT_STRESS)
                CF = 0.0*CF
                node.SetSolutionStepValue( KratosContact.CONTACT_STRESS, CF)


            if ( node.SolutionStepsDataHas( KratosContact.EFFECTIVE_CONTACT_FORCE)):
                CF = node.GetSolutionStepValue( KratosContact.EFFECTIVE_CONTACT_FORCE)
                CF = 0.0*CF
                node.SetSolutionStepValue( KratosContact.EFFECTIVE_CONTACT_FORCE, CF)

            if ( node.SolutionStepsDataHas( KratosContact.EFFECTIVE_CONTACT_STRESS)):
                CF = node.GetSolutionStepValue( KratosContact.EFFECTIVE_CONTACT_STRESS)
                CF = 0.0*CF
                node.SetSolutionStepValue( KratosContact.EFFECTIVE_CONTACT_STRESS, CF)


    def ExecuteFinalizeSolutionStep(self):

        time = self._GetStepTime()

        self._SetVariablesToZero()

        filling = self._GetInnerFilling() # distance from the toe to the maximum of the inner (r=0)

        self._ComputeForcesAndLengths()

        line = str(time) + " , " + str(filling) + " , "

        line = self._AppendForceAndLenght(line, self.InnerForce, self.InnerLength)
        line = self._AppendForceAndLenght(line, self.OuterForce, self.OuterLength)
        line = self._AppendForceAndLenght(line, self.TipForce)
        line = self._AppendForceAndLenght(line, self.PistonForce)

        line = line + " \n "

        figure_file = open(self.figure_path, "a")
        figure_file.write(line)
        figure_file.close()



    def _AppendForceAndLenght(self, line, Force, length = -55.5):

        for ii in range(0,3):
            line = line  + str(Force[ii]) + " , "

        if ( abs(length+55.5) > 1e-6):
            line = line  + str(length)  + " , "

        return line

    def _SetVariablesToZero(self):
        
        for ii in range(0,3):
            self.InnerForce[ii]  = 0.0
            self.OuterForce[ii]  = 0.0
            self.TipForce[ii]    = 0.0
            self.PistonForce[ii] = 0.0

        self.InnerLength = 0.0
        self.OuterLength = 0.0
        
        self.Filling = 0.0

    def  _GetStepTime(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.TIME]

    def _GetInnerFilling(self):

        YBest = -100000.0

        for node in self.model_part.GetNodes(0):
            if ( abs( node.X) < 1e-5):
                YThis = node.Y
                if ( YThis > YBest):
                    YBest = YThis

        time = self._GetStepTime()
        yBottom = self.BottomPoint[1] - time*self.velocity
        filling = YBest - yBottom

        return filling

    def _ComputeForcesAndLengths(self):

        time = self._GetStepTime()


        # Computation Inner Force

        y_lower  = self.InnerSeparationPoint[1] - time*self.velocity
        y_upper  = self.PistonPosition[1]
        x_center = self.BottomPoint[0]

        y_top = -100000.0

        for node in self.model_part.GetNodes(0):
            if (  (node.Y < y_upper) & (node.Y > y_lower) ):
                if ( node.X < x_center):
                    Force = node.GetSolutionStepValue( KratosMultiphysics.CONTACT_FORCE)
                    self.InnerForce = self.InnerForce + Force

                    if ( abs(Force[0]) > 1.0e-6):
                        if ( node.Y > y_top):
                            y_top = node.Y

        self.InnerLength = (y_top - y_lower)

        # Computation Outer Force

        y_lower = self.OuterSeparationPoint[1] - time * self.velocity
        x_center = self.BottomPoint[0]

        y_top = -100000.0

        for node in self.model_part.GetNodes(0):
            if (  node.Y > y_lower ):
                if (node.X > x_center):
                    Force = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)
                    self.OuterForce = self.OuterForce + Force

                    if (abs(Force[0]) > 1.0e-6):
                        if (node.Y > y_top):
                            y_top = node.Y

        self.OuterLength = (y_top - y_lower)

        # Computation Tip Force

        y_lower_left = self.InnerSeparationPoint[1] - time * self.velocity
        y_lower_right = self.OuterSeparationPoint[1] - time * self.velocity
        x_center = self.BottomPoint[0]

        for node in self.model_part.GetNodes(0):
            Force = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)
            if (abs(Force[0]) > 1.0e-6):

                if (node.X < x_center):
                    if ( node.Y < y_lower_left):
                        Force = node.GetSolutionStepValue( KratosMultiphysics.CONTACT_FORCE)
                        self.TipForce = self.TipForce + Force
                if (node.X > x_center):
                    if (node.Y < y_lower_right):
                        Force = node.GetSolutionStepValue( KratosMultiphysics.CONTACT_FORCE)
                        self.TipForce = self.TipForce + Force


        # Computation PistonForce

        y_lower = self.PistonPosition[1]
        x_center = self.BottomPoint[0]

        for node in self.model_part.GetNodes(0):
            if (node.X < x_center):
                if ( node.Y > y_lower):
                    Force = node.GetSolutionStepValue( KratosMultiphysics.CONTACT_FORCE)
                    self.PistonForce = self.PistonForce + Force


