from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
CheckForPreviousImport()

#import matplotlib
import collections

import numpy as np
#from pylab import *

import matplotlib.pyplot as plt


class InterpretCPTData:
    #

    def __init__(self, model_part, problem_path,Y0, Vy):

        self.model_part = model_part

        self.problem_path = problem_path
        self.mesh_id = 0

        self.Y0 = Y0
        self.Vy = Vy
        self.VyInitial = Vy
        self.radius = 0.017841
        self.Init = True;
    #
    def Initialize(self, mesh_id):

        print(' Arrived at CPTData.Initialize() ')
        
        self.mesh_id = mesh_id
        figure_path = os.path.join(self.problem_path, "CPTInterpreter.csv")


        if(os.path.exists(figure_path) == False):
            # print file headers
            figure_file = open(figure_path, "w")
            line_header  = "Time q_t f_s U2 U1 U3 dU2 dU1 dU3" +"\n"
            figure_file.write(line_header)
            figure_file.close()


    #
    def GetStepTime(self):

        return self.model_part.ProcessInfo[TIME]

    #
    def GetStepDeltaTime(self):

        return self.model_part.ProcessInfo[DELTA_TIME]

    #
    def GetStepVariable(self, variable):

        variable_value = []
        for node in self.model_part.GetNodes(self.mesh_id):
            kratos_variable = globals()[variable]
            kratos_variable = variable
            nodal_value = node.GetSolutionStepValue(kratos_variable);
            variable_value = variable_value + nodal_value

        return variable_value

    #
    def GetElemStepVariable(self, variable):

       elem_value = 0.0;
       elems = self.model_part.GetElements(self.mesh_id);
       elem = elems[self.Element];
       proc_info = self.model_part.ProcessInfo

       elem_value = elem.GetValuesOnIntegrationPoints(variable, proc_info);

       return elem_value[0][0]

    #
    def GetNodalStepVariable(self, variable):

       nodes = self.model_part.GetNodes(self.mesh_id);
       if (self.model_part.NumberOfNodes( self.mesh_id) < self.Node):
           self.Node = 10;
       node = nodes[self.Node];
       nodal_value = node.GetSolutionStepValue(variable);
       return nodal_value;
       if ( node.HasDofFor(variable) ):
          nodal_value = node.GetSolutionStepValue(variable);
       else:
          nodal_value = 0.0;
       return nodal_value

    #
    def GetPorePressureShaftAlt(self, YSearch, BLOKI = False):

        # get WATER_PRESSURE as variable
        nodes = self.model_part.GetNodes(self.mesh_id)
        if ( nodes[1].HasDofFor( WATER_PRESSURE) ):
            variable = WATER_PRESSURE
            variable2 = EXCESS_WATER_PRESSURE
        elif (nodes[1].HasDofFor( PRESSURE ) ):
            variable = PRESSURE
        else:
            return 0.0;

        # declare and initialize aux variables for closest top and bottom nodes
        XBestTop    =  100000000; # stores the X-coordinate of the closest node from top
        XBestBottom = -100000000; # stores the X-coordinate of the closest node from bottom
        YBestTop    =  100000000; # stores the y-coordinate of the closest node from top
        YBestBottom = -100000000; # stores the y-coordinate of the closest node from bottom
        nBestTop     = -100; # stores the ID of the closest node from top
        nBestBottom  = -100; # stores the ID of the closest node from bottom

        # loop over nodes of model_part to find nodes with contact force
        # and find closest node from bottom and top to YSearch position
        for node in self.model_part.GetNodes(self.mesh_id):
            Force = node.GetSolutionStepValue( CONTACT_FORCE);
            if ( abs(Force[0]) + abs(Force[1]) > 1e-8):
                XThis = node.X;
                YThis = node.Y;
                if (  (YThis-YSearch) > 0 ):
                    if ( abs(  YThis-YSearch )  < abs( YBestTop - YSearch) ):
                        nBestTop = node.Id;
                        XBestTop = XThis;
                        YBestTop = YThis;
                elif ( (YThis-YSearch) <= 0):
                    if ( abs( YThis-YSearch) < abs( YBestBottom - YSearch) ):
                        nBestBottom = node.Id
                        XBestBottom = XThis; 
                        YBestBottom = YThis;

        # output closest top and bottom nodes
        print(' top node ID: ' + str(nBestTop) + ', X: ' + str(XBestTop) + ', Y: ' + str(YBestTop))
        print(' bottom node ID: ' + str(nBestBottom) + ', X: ' + str(XBestBottom) + ', Y: ' + str(YBestBottom))

        if (  (nBestTop < 1) or (nBestBottom < 1) ):
            print( ' In the Usomething. NotFound Contacting nodes that are in the range of this kind of thing')
            return 0.0, 0.0

        ### check element start
        if(1==0):
            # Now i want to kwnow if there is an element that has this nodes
            # which might not be the case since first shoulder node is not in contact
            ReallyFound = False; 
            for elem in self.model_part.GetElements( self.mesh_id):
                a = [0, 1, 2]
                conec = [];
                found = 0
                for ii in a:
                    thisNode =  elem.GetNode(ii).Id
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
        ### check element end

        # interpolate between top and bottom nodes and return
        DeltaY = abs(YBestBottom - YBestTop);
        NTop = 1-  abs(YBestTop -YSearch) / DeltaY;
        NBottom = 1 - NTop;

        if (NTop > 1.0 or NTop < 0):
            print( 'ULTRA MEGA STUPID ERROR ')

        if (NBottom > 1.0 or NBottom < 0):
            print( 'ULTRA MEGA STUPID ERROR ')

        uBottom = nodes[nBestBottom].GetSolutionStepValue(variable)
        uBottom2 = nodes[nBestBottom].GetSolutionStepValue(variable2)
        uTop = nodes[nBestTop].GetSolutionStepValue(variable)
        uTop2 = nodes[nBestTop].GetSolutionStepValue(variable2)
        ThisV = NTop*uTop + NBottom*uBottom
        ThisV2 = NTop*uTop2 + NBottom*uBottom2
        return ThisV, ThisV2


    #
    def GetPorePressureU22(self):
        #print( ' Getting U2')
        YSearch = self.Y0 + self.Vy*self.GetStepTime()
        U22, dU22 = self.GetPorePressureShaftAlt( YSearch )
        return U22, dU22


    #
    def GetPorePressureU3(self):
        #print( ' Getting U3')
        YSearch = self.Y0 + self.Vy*self.GetStepTime()
        YSearch = YSearch + 7.5*self.radius
        U33, dU33 = self.GetPorePressureShaftAlt( YSearch )
        return U33, dU33


    #
    def GetPorePressureU1(self):
        #print( ' Getting U1')
        YSearch = self.Y0 + self.Vy*self.GetStepTime();
        #YSearch = YSearch - 0.01544989
        YSearch = YSearch - 0.017841
        U11, dU11 = self.GetPorePressureShaftAlt( YSearch )
        return U11, dU11


    #
    def GetResistance(self):
        result = [];
        YLim = self.Y0 + self.Vy * self.GetStepTime();
        for node in self.model_part.GetNodes(self.mesh_id):
            if ( node.Y <= YLim):
                Force = node.GetSolutionStepValue( CONTACT_FORCE);  ## NOT WORKING
                result = result + Force;

        return result[1];


    #
    def GetFriction(self):
        result = [];
        YMin = self.Y0  +  self.Vy * (self.GetStepTime() );
        YMax = YMin + 7.5*self.radius;

        for node in self.model_part.GetNodes(self.mesh_id):
            if (node.Y >= YMin):
                if (node.Y <= YMax):
                    Force = node.GetSolutionStepValue( CONTACT_FORCE);
                    result = result + Force;

        if (len(result) == 0):
            return 0.0

        return result[1];


    #
    def UpdateCPTVelocity(self, StoppingTime, Amplitude = 0, Acceleration = 1.0):
        if (self.model_part.ProcessInfo[TIME] > StoppingTime and StoppingTime > 1e-6 ):
            self.Vy = self.VyInitial * StoppingTime / self.model_part.ProcessInfo[TIME]
            if ( Amplitude > 0):
               w = 8.0;
               w = Acceleration * abs(self.VyInitial) / Amplitude
               wall_movement = self.VyInitial * StoppingTime + Amplitude * np.sin( w * ( self.model_part.ProcessInfo[TIME] - StoppingTime) )
               self.Vy = wall_movement / self.model_part.ProcessInfo[TIME]
               #print( ' PositionFromTheOther')
               #print( wall_movement )
               #print( ' VelocityFromTheOther')
               #print( self.Vy )


    #
    def SetStepResult(self, footingProblem = False):

        print(' Arrived at CPTData.SetStepResult() ')

        if (self.Init == False):
            for node in self.model_part.GetNodes(self.mesh_id): ## ?? why I do a loop ?? ##
             WallInfo = node.GetSolutionStepValue( WALL_REFERENCE_POINT );
             if (WallInfo[0] != 0):
                 self.radius = 2.0*WallInfo[0];
                 self.Init = True;

        # get step results
        Q = self.GetResistance();
        Friction = self.GetFriction();

        U22, dU22 = self.GetPorePressureU22();
        U33, dU33 = self.GetPorePressureU3();
        U11, dU11 = self.GetPorePressureU1();

        time = self.GetStepTime(); # - self.GetStepDeltaTime();

        ### write CPTInterpreter.csv
        figure_path = os.path.join(self.problem_path, "CPTInterpreter.csv")
        figure_file = open(figure_path, "a")
        line_value = str(time) + " " + str(Q) + " " + str(Friction) + " " + str(U22) +  " " + str(U11) + " " + str(U33) + " " + str(dU22) +  " " + str(dU11) + " " + str(dU33)  + "\n"
        figure_file.write(line_value)
        figure_file.close()

        ### contacting nodes
        figure_path = os.path.join( self.problem_path, "AllContactingNodes.csv")
        figure_file = open(figure_path, "a")
        line_value = str(" WRITTING A NEW TIME ") + str(time) + "\n"
        figure_file.write(line_value)
        for node in self.model_part.GetNodes(self.mesh_id):
            CF = node.GetSolutionStepValue( CONTACT_FORCE )
            if ( abs(CF[0]) + abs(CF[1]) > 1e-6):
                #if ( node.HasDofFor( CONTACT_STRESS_X ) ):
                #CS = node.GetSolutionStepValue( CONTACT_STRESS )
                CF = node.GetSolutionStepValue( CONTACT_FORCE)
                pw = 0;
                if ( node.HasDofFor( WATER_PRESSURE ) ):
                    pw = node.GetSolutionStepValue( WATER_PRESSURE )
                elif ( node.HasDofFor( PRESSURE ) ):
                    pw = node.GetSolutionStepValue( PRESSURE )
                x = node.X
                y = node.Y

                #line_value = str(time) + " " + str(x) + " " + str(y) + " " +  str( CS[0] ) + " " + str( CS[1] ) + " " + str(CF[0]) + " " +  str(CF[1]) + " " + str(pw) + "\n"
                line_value = str(time) + " " + str(x) + " " + str(y) + " " + " " + str(CF[0]) + " " +  str(CF[1]) + " " + str(pw) + "\n"
                figure_file.write(line_value)
        figure_file.close()

        ### U2 Position
        Yu2Position = self.Y0 + self.Vy*self.GetStepTime();
        Xu2Position = self.radius;
        figure_path = os.path.join( self.problem_path, "AllNearU2PositionNodes.csv")
        figure_file = open(figure_path, "a")
        line_value = str(" WRITTING A NEW TIME ") + str(time) + "\n"
        figure_file.write(line_value)

        for node in self.model_part.GetNodes(self.mesh_id):
            x = node.X;
            y = node.Y;

            distance = ( x - Xu2Position ) * ( x - Xu2Position);
            distance = distance + ( y - Yu2Position ) * ( y - Yu2Position);
            if ( distance < self.radius * self.radius):
                pw = 0
                if ( node.HasDofFor( WATER_PRESSURE ) ):
                    pw = node.GetSolutionStepValue( WATER_PRESSURE )
                elif ( node.HasDofFor( PRESSURE ) ):
                    pw = node.GetSolutionStepValue( PRESSURE )
                line_value = str(time) + " " + str(x) + " " + str(y) + " "+ str(pw) + "\n"
                figure_file.write(line_value)
        figure_file.close()

        ### U3 Position
        Yu2Position = self.Y0 + self.Vy*self.GetStepTime() + 7.5 * self.radius;
        Xu2Position = self.radius;
        figure_path = os.path.join( self.problem_path, "AllNearU3PositionNodes.csv")
        figure_file = open(figure_path, "a")
        line_value = str(" WRITTING A NEW TIME ") + str(time) + "\n"
        figure_file.write(line_value)

        for node in self.model_part.GetNodes(self.mesh_id):
            x = node.X;
            y = node.Y;

            distance = ( x - Xu2Position ) * ( x - Xu2Position);
            distance = distance + ( y - Yu2Position ) * ( y - Yu2Position);
            if ( distance < self.radius * self.radius):
                pw = 0
                if ( node.HasDofFor( WATER_PRESSURE ) ):
                    pw = node.GetSolutionStepValue( WATER_PRESSURE )
                elif ( node.HasDofFor( PRESSURE ) ):
                    pw = node.GetSolutionStepValue( PRESSURE )
                line_value = str(time) + " " + str(x) + " " + str(y) + " "+ str(pw) + "\n"
                figure_file.write(line_value)
        figure_file.close()

        ### U1 Position
        Yu2Position = self.Y0 + self.Vy*self.GetStepTime() -0.017841;
        Xu2Position = self.radius;
        figure_path = os.path.join( self.problem_path, "AllNearU1PositionNodes.csv")
        figure_file = open(figure_path, "a")
        line_value = str(" WRITTING A NEW TIME ") + str(time) + "\n"
        figure_file.write(line_value)

        for node in self.model_part.GetNodes(self.mesh_id):
            x = node.X;
            y = node.Y;

            distance = ( x - Xu2Position ) * ( x - Xu2Position);
            distance = distance + ( y - Yu2Position ) * ( y - Yu2Position);
            if ( distance < self.radius * self.radius):
                pw = 0
                if ( node.HasDofFor( WATER_PRESSURE ) ):
                    pw = node.GetSolutionStepValue( WATER_PRESSURE )
                elif ( node.HasDofFor( PRESSURE ) ):
                    pw = node.GetSolutionStepValue( PRESSURE )
                line_value = str(time) + " " + str(x) + " " + str(y) + " "+ str(pw) + "\n"
                figure_file.write(line_value)

        figure_file.close()
        

        
