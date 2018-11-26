from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPFEMSolid
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

import os.path

import numpy as np

def Factory( custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ConePenetrationUtility(Model, custom_settings["Parameters"])


class ConePenetrationUtility(KratosMultiphysics.Process):
    #

    def __init__(self, model_part, custom_settings):

        KratosMultiphysics.Process.__init__(self)


        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
         "cone_radius"        : 0.01784,
         "u2_initial_depth"   : 0.0332,
         "velocity"           : 0.02,
         "dissipation_set"    : true,
         "dissipation_depth"  : 0.2676
        }
        """)

        settings = custom_settings
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model_part['Main_Domain']
        
        self.Y0 = settings["u2_initial_depth"].GetDouble()
        self.Vy = settings["velocity"].GetDouble()
        self.radius = settings["cone_radius"].GetDouble()

        self.dissipation_set = settings["dissipation_set"].GetBool()
        if ( self.dissipation_set):
            self.dissipation_depth = settings["dissipation_depth"].GetDouble()
        
        
        problem_path = os.getcwd()
        self.figure_path = os.path.join(problem_path, "cone_penetration_data.csv")


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

            if ( node.SolutionStepsDataHas( KratosSolid.DISPLACEMENT_REACTION)):
                #print(node)
                for step in range(0,2):
                    CF = node.GetSolutionStepValue( KratosSolid.DISPLACEMENT_REACTION, step)
                    normal = node.GetSolutionStepValue( KratosMultiphysics.NORMAL, step);
                    #print('NODE', node.Id, ' normal', normal)
                    #print('  step ', step, 'REACTION', CF)
                    CF = 0.0*CF
                    #print('  step ', step, 'REACTION', CF)
                    node.SetSolutionStepValue( KratosSolid.DISPLACEMENT_REACTION, step, CF)
                #print('----------')

    #
    def ExecuteFinalizeSolutionStep(self):

        Q = self._GetResistance();
        Friction = self._GetFriction();

        U22 = self._GetPorePressureU2();
        U33 = self._GetPorePressureU3()
        U11 = self._GetPorePressureU1();

        time = self._GetStepTime(); # - self.GetStepDeltaTime();

        line_value = str(time) + " , " + str(Q) + " , " + str(Friction) + " , " + str(U22) +  " , " + str(U11) + " , " + str(U33)  + " \n"

        figure_file = open(self.figure_path, "a")
        figure_file.write(line_value)
        figure_file.close()



        
    #
    def _GetStepTime(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.TIME]

    #
    def _GetStepDeltaTime(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

    def _GetU2Position(self):
        avance = self.Vy*self._GetStepTime()
        if ( self.dissipation_set):
            if ( avance > self.dissipation_depth):
                avance = self.dissipation_depth
        YSearch = self.Y0 - avance;
        return YSearch


    def _GetPorePressureU2(self):
        YSearch = self._GetU2Position()
        U22 = self._GetPorePressureShaft( YSearch)
        return U22

    def _GetPorePressureU3(self):
        YSearch = self._GetU2Position()
        YSearch = YSearch + 7.5*self.radius
        U33 = self._GetPorePressureShaft( YSearch )
        return U33

    def _GetPorePressureU1(self):
        YSearch = self._GetU2Position()
        YSearch = YSearch - self.radius / np.tan(30.0*3.14159/180.0)/2.0
        U11 = self._GetPorePressureShaft( YSearch )
        return U11


    def _GetResistance(self):
        result = 0.0;
        YLim = self.Y0 - self.Vy * self._GetStepTime();
        YLim = self._GetU2Position()
        for node in self.model_part.GetNodes(0):
            if ( node.Y < YLim):
                Force = node.GetSolutionStepValue( KratosMultiphysics.CONTACT_FORCE);  
                result = result + Force[1];

        return result;

    #
    def _GetFriction(self):
        result = 0;
        YMin = self._GetU2Position()
        YMax = YMin + 7.5*self.radius;

        for node in self.model_part.GetNodes(0):
            if (node.Y >= YMin):
                if (node.Y <= YMax):
                    Force = node.GetSolutionStepValue( KratosMultiphysics.CONTACT_FORCE);
                    result = result + Force[1];

        return result;


    def _GetPorePressureShaft(self, YSearch):

        nodes = self.model_part.GetNodes();
        if ( nodes[3].HasDofFor( KratosMultiphysics.WATER_PRESSURE) ):
            variable = KratosMultiphysics.WATER_PRESSURE;
        elif (nodes[3].HasDofFor( KratosMultiphysics.PRESSURE ) ):
            variable = KratosMultiphysics.PRESSURE;
        else:
            return 0.0;

        YBestTop    =  100000000;
        YBestBottom = -100000000;
        nBestTop     = -100;
        nBestBottom  = -100;

        for node in self.model_part.GetNodes(0):
            Force = node.GetSolutionStepValue( KratosMultiphysics.CONTACT_FORCE);
            if ( abs(Force[0]) + abs(Force[1]) > 1e-8):
                YThis = node.Y;
                if (  (YThis-YSearch) > 0 ):
                    if ( abs(  YThis-YSearch )  < abs( YBestTop - YSearch) ):
                        nBestTop = node.Id;
                        YBestTop = YThis;
                elif ( (YThis-YSearch) <= 0):
                    if ( abs( YThis-YSearch) < abs( YBestBottom - YSearch) ):
                        nBestBottom = node.Id
                        YBestBottom = YThis; 

        if (  (nBestTop < 1) or (nBestBottom < 1) ):
            print( ' In the Usomething. NotFound Contacting nodes that are in the range of this kind of thing')
            return 0.0

        # Now i want to kwnow if there is an element that has these nodes
        ReallyFound = False; 
        for elem in self.model_part.GetElements(0):
            a = [1, 2, 3]
            conec = [];
            found = 0
            for thisNodeGeom in elem.GetNodes():
                thisNode =  thisNodeGeom.Id
                if (thisNode == nBestTop):
                    found = found + 1
                if (thisNode == nBestBottom ):
                    found = found + 1
            if ( found == 2):
                ReallyFound = True;
                break

        if ( ReallyFound == False):
            print( ' In the U Something. The two nodes do not share an element ')
            return 0.0

        DeltaY = abs(YBestBottom - YBestTop);
        NTop = 1-  abs(YBestTop -YSearch) / DeltaY;
        NBottom = 1 - NTop;

        if (NTop > 1.0 or NTop < 0):
            print( 'ULTRA MEGA STUPID ERROR ')

        if (NBottom > 1.0 or NBottom < 0):
            print( 'ULTRA MEGA STUPID ERROR ')


        uBottom = nodes[nBestBottom].GetSolutionStepValue(variable);
        uTop = nodes[nBestTop].GetSolutionStepValue(variable);
        ThisV = NTop*uTop + NBottom*uBottom;
        return ThisV;

