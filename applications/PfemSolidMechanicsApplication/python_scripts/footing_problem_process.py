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
    return FootingProblemProcess(Model, custom_settings["Parameters"])


class FootingProblemProcess(KratosMultiphysics.Process):

    #
    def __init__(self, model_part, custom_settings):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
          "tunnel_radius": 0.041,
          "boxHeight": 0.15
        }
        """)

        settings = custom_settings
        settings.ValidateAndAssignDefaults(default_settings)

        self.tunnel_radius = settings["tunnel_radius"].GetDouble()
        self.boxHeight = settings["boxHeight"].GetDouble()


        self.model_part = model_part['Main_Domain']
        print(self.model_part)

        # initialize figure path
        problem_path = os.getcwd()
        self.csv_path = os.path.join(problem_path, "footing.csv")
        csv_file = open(self.csv_path, "a")
        line = " 0.0 , 0.0 , 0.0 , 0.0  \n"
        csv_file.write(line)
        csv_file.close()

        self.csv_path_tunnel = os.path.join(problem_path, "tunnel.csv")
        csv_file = open(self.csv_path_tunnel, "a")
        line = " 0.0 , 0.0 ,   \n"
        csv_file.write(line)
        csv_file.close()

    def ExecuteFinalizeSolutionStep(self):
      
        model_part = self.model_part.GetSubModelPart("Solid_group_DeformableBodies-auto-1")

        time = self._GetStepTime()
        
        SomeProcessInfo = KratosMultiphysics.ProcessInfo()

        Force = 0.0*KratosMultiphysics.Array3()
        for node in model_part.GetNodes(0): 
            if ( abs(node.Y0) < 1.0e-5):
                if (node.X0 < 1.0001):
                    res = node.GetSolutionStepValue(KratosSolid.DISPLACEMENT_REACTION)
                    Force = Force + res

        line = str(time) + " , "

        line = self._AddToLine( line, Force)
        line = line + " \n"

        csv_file = open(self.csv_path, "a")
        csv_file.write(line)
        csv_file.close()



        #The same but this time for tunneling
        import math
        NormalForce = 0.0

        Radius = self.tunnel_radius * 1.3;
        for node in model_part.GetNodes(0): 
            if ( abs(node.Y0) < Radius and abs(node.X0) < Radius):
                thisContribution = 0.0;
                if ( node.SolutionStepsDataHas(KratosMultiphysics.CONTACT_FORCE)):
                    res = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)
                    thisContribution = ( res[0]*res[0]+ res[1]*res[1] )**0.5 
                elif (  node.SolutionStepsDataHas( KratosSolid.DISPLACEMENT_REACTION) ):
                    res = node.GetSolutionStepValue(KratosSolid.DISPLACEMENT_REACTION)
                    print('MoreWorkIsRequired')
                NormalForce = NormalForce + thisContribution
        line = str(time) + " , " + str(NormalForce) + " , " 
        line = line + " \n"

        csv_file = open(self.csv_path_tunnel, "a")
        csv_file.write(line)
        csv_file.close()

        #TryToPlotSomethingsimilar than the tunnel
        problem_path = os.getcwd()
        csv_path_tunnel = os.path.join(problem_path, "tunnel3.csv")
        csv_file = open(csv_path_tunnel, "a")

        time = self._GetStepTime()
        line = str(time) + " \n "
        csv_file.write(line)

        Radius = self.tunnel_radius;
        for node in model_part.GetNodes(0): 
            position = node.Y0*node.Y0 + node.X0 *node.X0
            position = position**0.5
            if ( position > 0.99*Radius and position < 1.01*Radius):
                CF  = 0.0*KratosMultiphysics.Array3();
                DR  = 0.0*KratosMultiphysics.Array3();
                PR  = 0.0*KratosMultiphysics.Array3();

                ToWrite = False;
                if ( node.SolutionStepsDataHas(KratosMultiphysics.CONTACT_FORCE)):
                    CF = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)
                    ToWrite = True;
                elif (  node.SolutionStepsDataHas( KratosSolid.DISPLACEMENT_REACTION)):
                    DR = node.GetSolutionStepValue(KratosSolid.DISPLACEMENT_REACTION)
                    ToWrite = True;
                elif (  node.SolutionStepsDataHas( KratosMultiphsycics.POSITIVE_FACE_PRESSURE_VECTOR) ):
                    PR = node.GetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE_VECTOR)
                    ToWrite = True;

                if ( ToWrite) :
                    X = node.X;
                    Y = node.Y;
                    disp = node.GetSolutionStepValue( KratosMultiphysics.DISPLACEMENT)
                    ux = disp[0];
                    uy = disp[1];
                    line = str(X) + " , " + str(Y) + " , " + str(ux) + " , " + str(uy) + " , "

                    if ( (abs(CF[0])+abs(CF[1])>1e-5)):
                        line = line +  str(CF[0]) + " , " + str(CF[1]) 
                    elif ( (abs(DR[0])+abs(DR[1])>1e-5)):
                        line = line + str(DR[0]) + " , " + str(DR[1]) 
                    elif ( (abs(PR[0])+abs(PR[1])>1e-5)):
                        line = line + str(PR[0]) + " , " + str(PR[1]) 
                    else:
                        line = line + str(0.0) + " , " + str(0.0)
                    line = line + " \n"
                    csv_file.write(line)


        csv_file.close()


        #TryToPlotSomethingsimilar than the tunnel
        problem_path = os.getcwd()
        csv_path_tunnel = os.path.join(problem_path, "tunnel4.csv")
        csv_file = open(csv_path_tunnel, "a")

        time = self._GetStepTime()
        line = str(time) + " \n "
        csv_file.write(line)

        for node in model_part.GetNodes(0): 
            if ( abs(node.Y0) > self.boxHeight - 1e-5):

                X = node.X;
                Y = node.Y;
                disp = node.GetSolutionStepValue( KratosMultiphysics.DISPLACEMENT)
                ux = disp[0];
                uy = disp[1];
                line = str(X) + " , " + str(Y) + " , " + str(ux) + " , " +  str(uy)
                line = line + " \n"
                csv_file.write(line)


        csv_file.close()



    def _AddToLine(self, line, vector):

        for i in vector:
            line = line + str(i) + " , "

        return line

    def  _GetStepTime(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.TIME]

