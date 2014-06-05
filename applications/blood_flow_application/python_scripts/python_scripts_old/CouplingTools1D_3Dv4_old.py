from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
CheckForPreviousImport()

import math
import time
import sys
import config


class TransferTools:

    def __init__(self, model_part_1d, model_part_3d):
        self.model_part_1d = model_part_1d
        self.model_part_3d = model_part_3d
    # 3

    def Initialize(self):
        self.inlets_1d = []
        self.outlets_1d = []
        self.inlets_3d = []
        self.outlets_3d = []
        self.inlet_areas_3d = []
        self.outlet_areas_3d = []
        self.inlet_velocity_directions = []
        self.flow_1d = 0.0
        self.flow_1d_pressure = 0.0
        self.flow_3d_in = 0.0
        self.flow_3d_in_pressure = 0.0
        self.flow_3d_out = 0.0
        self.flow_1d_out = 0.0
        self.flow_1d_in2 = 0.0
        # compute normals and 3d areas
        #NormalCalculationUtils().SwapNormals(self.model_part_3d)
        NormalCalculationUtils().CalculateOnSimplex(self.model_part_3d,3)
        
        #BodyNormalCalculationUtils().CalculateBodyNormals(
            #self.model_part_3d, 3)

        # detect 1d inlets

        for i in range(100, 101):
            nfound = 0
            aux = []
            for node in self.model_part_1d.Nodes:
                if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
                    aux.append(node)
                    # node.Fix(NODAL_AREA)
                    # node.Fix(FLOW)
                    print "1D_inlet has been assigned", node
                    # raw_input()
            if(len(aux) != 0):
                self.inlets_1d.append(aux)
            else:
                break

        if(len(self.inlets_1d) == 0):
            # print "number of 1d is zero"
            sys.exit(
                "number of 1d_inlets are zero!! Please check your config.py file!")

        # detect 3d inlets
        for i in range(100, 101):
            nfound = 0
            aux = []
            directions = []
            for node in self.model_part_3d.Nodes:
                if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
                    aux.append(node)
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    node.Fix(VELOCITY_Z)
                    tmp = node.GetSolutionStepValue(NORMAL)
                    normN = math.sqrt(tmp[0] ** 2 + tmp[1] ** 2 + tmp[2] ** 2)
                    tmp /= normN
                    directions.append(tmp)
                    print "3D_inlet has been assigned", node
            if(len(aux) != 0):
                self.inlets_3d.append(aux)
                self.inlet_velocity_directions.append(directions)
            else:
                break
            # raw_input()

        if(len(self.inlets_3d) == 0):
            sys.exit(
                "number of 3d_inlets are zero!! Please check your config.py file!")
            # print "number of 3d is zero"
            # err

        if(len(self.inlets_1d) != len(self.inlets_3d)):
            sys.exit(
                "number of 1d and 3d inlets is different!! Please check your config.py file!")
            # print "number of 1d and 3d inlets is different"
            # err

        # detect 1d outlets
        for i in range(1001, 1199):
            nfound = 0
            aux = []
            for node in self.model_part_1d.Nodes:
                if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
                    aux.append(node)
                    node.Fix(FLOW)
            if(len(aux) != 0):
                self.outlets_1d.append(aux)
                print "1D_outlet has been assigned(H)", node
                # raw_input()
            else:
                break

        if(len(self.outlets_1d) == 0):
            # print "number of 1d is zero"
            sys.exit(
                "number of 1d_outlets are zero!! Please check your config.py file!")

        # detect 3d outlets
        for i in range(1001, 1199):
            nfound = 0
            aux = []
            for node in self.model_part_3d.Nodes:
                if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
                    aux.append(node)
                    node.Fix(PRESSURE)
            if(len(aux) != 0):
                self.outlets_3d.append(aux)
            else:
                break

        if(len(self.outlets_3d) == 0):
            # print "number of 1d is zero"
            sys.exit(
                "number of 3d_outlets are zero!! Please check your config.py file!")

        for i in range(0, len(self.inlets_3d)):
            inlet_nodes_3d = self.inlets_3d[i]
            area3d = 0.0
            for node in inlet_nodes_3d:
                n = node.GetSolutionStepValue(NORMAL)
                a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                area3d += a
            self.inlet_areas_3d.append(area3d)
            # print "area3d_inlet", area3d

        for i in range(0, len(self.outlets_3d)):
            inlet_nodes_3d = self.outlets_3d[i]
            area3d = 0.0
            for node in inlet_nodes_3d:
                n = node.GetSolutionStepValue(NORMAL)
                a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                area3d += a
            self.outlet_areas_3d.append(area3d)
            print "area_inlet_3d", area3d
            #raw_input()

    def Initial_Contitions(self):
        initial_pressure = 0  # TODO
        print "Inicializo 3D"
        for i in range(0, len(self.inlets_1d)):
            inlet_nodes_1d = self.inlets_1d[i]
            print "3D-1D: inlet_nodes_1d [0].Id::::::>>>> ", inlet_nodes_1d[0].Id
            for i in range(0, len(self.outlets_1d)):
                outlet_nodes_1d = self.outlets_1d[i]
                pressinlet3D = outlet_nodes_1d[
                    0].GetSolutionStepValue(PRESSURE)
                print pressinlet3D
                print "3D-1D: outlet_nodes_1d[0].Id::::::>>>> ", outlet_nodes_1d[0].Id
                for i in range(0, len(self.inlets_3d)):
                    inlet_nodes_3d = self.inlets_3d[i]
                    area3d = self.inlet_areas_3d[i]
                    print "area_inlet_3d", area3d
                    directions = self.inlet_velocity_directions[i]
                    radio3d = math.sqrt(area3d * 3.1416)
                    vel1d = inlet_nodes_1d[
                        0].GetSolutionStepValue(FLOW) / area3d
                    k = 0
                    print vel1d
                    # Impongo velocidad y la presion
                    for node in inlet_nodes_3d:
                        n = node.GetSolutionStepValue(NORMAL)
                        a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                        orientation = directions[k]
                        node.SetSolutionStepValue(
                            VELOCITY, 0, directions[k] * vel1d)
                        # node.SetSolutionStepValue(PRESSURE, 0, pressinlet3D)
                        k = k + 1
            # raw_input()

    def Transfer1D_to_3D(self):

            # ARCHIVE TO SET :::::::::::::::::::::::::::>>>>>>>>>>>>>> VARIABLES
            # import config_full
        initial_pressure = config.systolic_pressure
        initial_pressure = 0

        print "Transfer1D_to_3D"
        for i in range(0, len(self.inlets_1d)):
            inlet_nodes_1d = self.inlets_1d[i]
            print "NODO 1D-3D:: inlet_nodes_1d[0].Id::::::inlet 1D coupled with the 3D inlet>>>> ",inlet_nodes_1d[0].Id
            # print "--"
            # raw_input()
            for i in range(0, len(self.inlets_3d)):
                inlet_nodes_3d = self.inlets_3d[i]
                area3d = self.inlet_areas_3d[i]
                directions = self.inlet_velocity_directions[i]
                radio3d = math.sqrt(area3d * 3.1416)
                vel1d = inlet_nodes_1d[0].GetSolutionStepValue(FLOW) / area3d
                # NOTA: REVISAR EL CAUDAL QUE ESTAMOS IMPONIENDO DE ENTRADA El hecho de considerar q/2 es debido a que el caudal que imponemos
                # pressinlet1D =
                # inlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
                k = 0
                # Impongo velocidad y la presion
                for node in inlet_nodes_3d:
                    n = node.GetSolutionStepValue(NORMAL)
                    a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                    orientation = directions[k]
                    node.SetSolutionStepValue(VELOCITY, 0, directions[k] * vel1d)
                    k = k + 1
            print "area3d",area3d ,"area1d",inlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA), "del nodo", inlet_nodes_1d[0].Id
            print "Estoy fijando la velocidad en el inlet del 3D que proviene del nodo " , inlet_nodes_1d[0].Id, " del 1D. La velocidad que estoy fijando es"
            print "vel1d",vel1d, "del nodo", inlet_nodes_1d[0].Id
            

        for i in range(0, len(self.outlets_1d)):
            outlet_nodes_1d = self.outlets_1d[i]
            pressinlet3D = outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
            # pressinlet3D = 0
            print "NODO 1D-3D: outlet_nodes_1d[0].Id::::::outlet 1D coupled with the 3D outlet>>>> ",outlet_nodes_1d[0].Id, "del nodo", outlet_nodes_1d[0].Id
            # print "pressinlet3D ",pressinlet3D
            # raw_input()
            # double beta = E0*thickness0*1.77245385/(1.0-nu0*nu0);
            beta = outlet_nodes_1d[0].GetSolutionStepValue(BETA)
            #beta = ((outlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * outlet_nodes_1d[0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)) / (
                #1 - (outlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO) * outlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO))))
            A = outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            print "AREA_1D", A
            A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
            press = initial_pressure + beta * \
                (math.sqrt(A) - math.sqrt(A0)) / \
                A0  # math.sqrt(A/A0)*beta - beta
            #press = initial_pressure + beta * \
                #(math.sqrt(A) - math.sqrt(A0)) / \
                #A0  # math.sqrt(A/A0)*beta - beta
            print "prress calucalteed",press 
            print "press solver", pressinlet3D
            for i in range(0, len(self.outlets_3d)):
                outlet_nodes_3d = self.outlets_3d[i]
                area3d = self.outlet_areas_3d[i]
                for node in outlet_nodes_3d:
                    node.SetSolutionStepValue(PRESSURE, 0, press)     
            print "area3d",area3d ,"area1d",outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            print "Estoy fijando la presion en el outlet del 3D que proviene del nodo " , outlet_nodes_1d[0].Id, " del 1D. La presion que estoy fijando es"
            print "pression",press, "del nodo", outlet_nodes_1d[0].Id

    def Transfer3D_to_1D(self):
        inlet_flow = (self.inlets_1d[0])[0].GetSolutionStepValue(FLOW)
        self.outlets_1d[0][0].SetSolutionStepValue(FLOW,0,inlet_flow)

        #outlet_area = self.outlets_3d[0][0].GetSolutionStepValue(NODAL_AREA)
        #outlet_beta = self.outlets_3d[0][0].GetSolutionStepValue(BETA)
        outlet_pressure = self.outlets_1d[0][0].GetSolutionStepValue(PRESSURE)

        inlet_beta = self.inlets_1d[0][0].GetSolutionStepValue(BETA)
        inlet_A0  = self.inlets_1d[0][0].GetValue(NODAL_AREA)
        Ainlet_to_prescribe = (((outlet_pressure * inlet_A0)/inlet_beta) + math.sqrt(inlet_A0))**2  # A0*(avg_press/beta + 1)**2
        self.outlets_1d[0][0].SetSolutionStepValue(NODAL_AREA,0,Ainlet_to_prescribe)

        print "inlet node = ",self.inlets_1d[0][0].Id
        print "outlet node = ",self.outlets_1d[0][0].Id
		
		
    #def Transfer3D_to_1D(self):
        #for i in range(0, len(self.inlets_3d)):
            #inlet_nodes_1d = self.inlets_1d[i]
            #inlet_nodes_3d = self.inlets_3d[i]
            #area3d = self.inlet_areas_3d[i]
            #outlet_nodes_1d = self.outlets_1d[i]
            #press_3d = 0.0
            #counter = 0.0
            #for node in inlet_nodes_3d:
                #press_3d += node.GetSolutionStepValue(PRESSURE)
                #counter += 1.0
            #avg_press = press_3d / counter
            ##beta = inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * inlet_nodes_1d[
                ##0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)
            ##beta = beta / \
                ##(1.0 - (inlet_nodes_1d[0].GetSolutionStepValue(
                    ##POISSON_RATIO) * inlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO)))
            ##print "beta", beta
            #beta = inlet_nodes_1d[0].GetSolutionStepValue(BETA)
            #print "beta", beta
            #print "inlet node is ", inlet_nodes_1d[0].Id
            #beta = inlet_nodes_1d[0].GetSolutionStepValue(BETA)
            #A0 = inlet_nodes_1d[0].GetValue(NODAL_AREA)
            #Ainlet_to_prescribe = (((avg_press * A0)/beta) + math.sqrt(A0))**2  # A0*(avg_press/beta + 1)**2
            #print "in Transfer AREA 3D_to_1D----->", Ainlet_to_prescribe, " to node" , inlet_nodes_1d[0].Id
            #print "Area del 1D Original", A0 
            #print "Area del 1D Modificada", inlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            
           
            ## A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
            ## outlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA ,0, A0)
            #inlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA, 0, Ainlet_to_prescribe)

            ## compute flow on inlet
            #flow = 0.0
            #for node in inlet_nodes_3d:
                #normal = node.GetSolutionStepValue(NORMAL)
                #vel = node.GetSolutionStepValue(VELOCITY)
                #flow += normal[0] * vel[0] + normal[
                    #1] * vel[1] + normal[2] * vel[2]
            ## self.flow_3d_in = str(flow)
            #print "TRASFER to 3D to 1D:::flow entering= ",flow
            #print "Flow de 1D inlet",  inlet_nodes_1d[0].GetSolutionStepValue(FLOW)
            ## flow_aux=  flow

        #for i in range(0, len(self.outlets_3d)):
            #outlet_nodes_1d = self.outlets_1d[i]
            #outlet_nodes_3d = self.outlets_3d[i]
            #area3d = self.outlet_areas_3d[i]

            ## compute flow on 3D outlet
            #flow = 0.0
            #for node in outlet_nodes_3d:
                #normal = node.GetSolutionStepValue(NORMAL)
                #vel = node.GetSolutionStepValue(VELOCITY)
                #flow += normal[0] * vel[0] + normal[
                    #1] * vel[1] + normal[2] * vel[2]

            ## assign flow to outlet
            #outlet_nodes_1d[0].SetSolutionStepValue(FLOW, 0, -flow)
            #print " velocity del 3D", vel
            #print " velocity del 1D", vel
            #print "in Transfer Flow 3D_to_1D----->", -flow, " to node" , outlet_nodes_1d[0].Id

        # raw_input()
#-------------------------------------------------------------------------

# Setting Contitions 3d

#-------------------------------------------------------------------------

    def Setting3d(self):
        # print self
     # print len(self.model_part_3d.Conditions)
        for node in self.model_part_3d.Nodes:
            node.Free(VELOCITY_X)
            node.Free(VELOCITY_Y)
            node.Free(VELOCITY_Z)
            node.Free(PRESSURE)
            node.SetSolutionStepValue(VISCOSITY, 0, 0.0035 / 1060.0)
            node.SetSolutionStepValue(DENSITY, 0, 1060.0)
        # set inlet
        for cond in self.model_part_3d.Conditions:
            if(cond.Properties.Id == 100):  # inlet
                for node in cond.GetNodes():
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    node.Fix(VELOCITY_Z)
                    node.SetSolutionStepValue(FLAG_VARIABLE, 0, 100.0)
            # if(cond.Properties.Id > 1000): ##outlet

                # for node in cond.GetNodes():

                    # node.Fix(PRESSURE)

                    # node.SetSolutionStepValue(FLAG_VARIABLE,0,cond.Properties.Id)

            # if(cond.Properties.Id > 1001): ##outlet

                # for node in cond.GetNodes():

                    # node.Fix(PRESSURE)

                    # node.SetSolutionStepValue(FLAG_VARIABLE,0,1002.0)

        # set sides (overwrites the inlet)
        for cond in self.model_part_3d.Conditions:
            if(cond.Properties.Id == 1):  # sides --> note that this is done in an outer separated loop!!
                for node in cond.GetNodes():
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    node.Fix(VELOCITY_Z)
                    node.SetSolutionStepValue(FLAG_VARIABLE, 0, 0.0)

        # set output (overwrites the others)

        for cond in self.model_part_3d.Conditions:
            if(cond.Properties.Id > 1000):  # outlet
                for node in cond.GetNodes():
                    # if(not node.IsFixed(VELOCITY_X)):
                    node.Fix(PRESSURE)
                    node.SetSolutionStepValue(
                        FLAG_VARIABLE, 0, cond.Properties.Id)

     # for cond in self.model_part_3d.Conditions:
            # if(cond.GetValue(IS_STRUCTURE) == 1):
                # for node in cond.GetNodes():
                    # node.Fix(VELOCITY_X)
                    # node.Fix(VELOCITY_Y)
                    # node.Fix(VELOCITY_Z)
                    # node.SetSolutionStepValue(VELOCITY_X,0,0.0)
        # this is only in this example...it should be set up by the problemtype
        # for node in self.model_part_3d.Nodes:
            # node.SetSolutionStepValue(VISCOSITY,0,0.0035/1060.0)
            # node.SetSolutionStepValue(DENSITY,0,1060.0)
            # if(node.IsFixed(VELOCITY_X) == True and
            # node.GetSolutionStepValue(VELOCITY_X) > 0.0001):
                # node.SetSolutionStepValue(FLAG_VARIABLE,0,100.0)
            # if(node.IsFixed(PRESSURE) == True):
                # node.SetSolutionStepValue(FLAG_VARIABLE,0,101.0)
            # node.SetSolutionStepValue(VELOCITY_X,0,0.0)
            # copy Y_WALL
        for node in self.model_part_3d.Nodes:
            y = node.GetSolutionStepValue(Y_WALL, 0)
            node.SetValue(Y_WALL, y)

        counter = 0.0

        for node in self.model_part_3d.Nodes:
            if(node.IsFixed(PRESSURE)):
                counter += 1.0


        # print "n pressure nodes ",counter
