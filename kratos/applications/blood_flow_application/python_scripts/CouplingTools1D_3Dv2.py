from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
CheckForPreviousImport()

import math
import time
import sys


class TransferTools:

    def __init__(self, model_part_1d, model_part_3d):
        self.model_part_1d = model_part_1d
        self.model_part_3d = model_part_3d
        # self.f = open("perdidas.txt","w")
        # self.f1d = open("1d.txt","w")
        # self.f3d = open("3d.txt","w")

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
        # self.f.write("EMPIEZO\n")
        # self.f1d.write("EMPIEZO\n")
        # self.f3d.write("EMPIEZO\n")

        # compute normals and 3d areas
        BodyNormalCalculationUtils().CalculateBodyNormals(
            self.model_part_3d, 3)

        # El nodo 23, seria el nodo auxiliar, para la entrada del 3D (1Dto3D ==> 11-23)
        # model_part1D.Nodes[23].SetSolutionStepValue(FLAG_VARIABLE,0,100)
        # El nodo 10, seria el nodo del 1D, de la salida del 3D (3Dto1D ==> (de momento esta desactivado, en el mpda ??)
        # model_part1D.Nodes[10].SetSolutionStepValue(FLAG_VARIABLE,0,101)

        # detect 1d inlets
        for i in range(100, 101):
            nfound = 0
            aux = []
            for node in self.model_part_1d.Nodes:
                if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
                    aux.append(node)
                    node.Fix(NODAL_AREA)
                    node.Fix(FLOW)
                    # print "1D_inlet has been assigned", node
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
                    # print "3D_inlet has been assigned", node
            if(len(aux) != 0):
                self.inlets_3d.append(aux)
                self.inlet_velocity_directions.append(directions)
            else:
                break

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

        # try len(self.inlets_1d) == len(self.inlets_3d):
            # print "found ",len(inlets_1d)," domains"
        # except RuntimeError:
            # print "number of 1d and 3d inlets is different"

        for i in range(0, len(self.inlets_3d)):
            inlet_nodes_3d = self.inlets_3d[i]
            area3d = 0.0
            for node in inlet_nodes_3d:
                n = node.GetSolutionStepValue(NORMAL)
                a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                area3d += a
            self.inlet_areas_3d.append(area3d)

        for i in range(0, len(self.outlets_3d)):
            inlet_nodes_3d = self.outlets_3d[i]
            area3d = 0.0
            for node in inlet_nodes_3d:
                n = node.GetSolutionStepValue(NORMAL)
                a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                area3d += a
            self.outlet_areas_3d.append(area3d)

        # print "inlet areas 1d = ",self.inlets_1d[0][0].GetValue(NODAL_AREA)
        # print "outlet areas 1d = ",self.outlets_1d[0][0].GetValue(NODAL_AREA)
        # print "inlets_1d = ",self.inlets_1d
        # print "inlet_velocity_directions = ",self.inlet_velocity_directions
        # print "outlets_3d = ",self.outlets_3d
        # print "outlet areas 3d = ",self.outlet_areas_3d
        # print "outlet areas 3d = ",self.outlet_areas_3d
        # raw_input()
        # for direc in self.inlet_velocity_directions[0]:
            # print direc[0]," ",direc[1]," ",direc[2]

        # for node in self.inlets_3d[0]:
            # print node.Id

    def Transfer1D_to_3D(self):
        # ARCHIVE TO SET :::::::::::::::::::::::::::>>>>>>>>>>>>>> VARIABLES
        # import config_full
        # initial_pressure=config_full.initial_pressure
        initial_pressure = 0

        for i in range(0, len(self.inlets_3d)):
            inlet_nodes_1d = self.inlets_1d[i]
            inlet_nodes_3d = self.inlets_3d[i]
            area3d = self.inlet_areas_3d[i]
            directions = self.inlet_velocity_directions[i]

            # self.f1d.write("FLOW 1D --> 3D (inlet_3d)--->   ")
            # self.f1d.write(str(inlet_nodes_1d[0].GetSolutionStepValue(FLOW)))
            # self.f1d.write("\n")

            self.flow_1d = str(inlet_nodes_1d[0].GetSolutionStepValue(FLOW))

            vel1d = inlet_nodes_1d[0].GetSolutionStepValue(FLOW) / area3d

            # print "inlet_nodes_1d[0].Id::::::>>>> ",inlet_nodes_1d[0].Id
            # print "velocity  1d to 3D ::::::> ", vel1d
            # print "flow         1d to 3D ::::::>",
            # inlet_nodes_1d[0].GetSolutionStepValue(FLOW)

            k = 0

            # Impongo velocidad
            for node in inlet_nodes_3d:
                orientation = directions[k]
                node.SetSolutionStepValue(VELOCITY, 0, directions[k] * vel1d)
                k = k + 1

        for i in range(0, len(self.outlets_3d)):
            outlet_nodes_1d = self.outlets_1d[i]
            outlet_nodes_3d = self.outlets_3d[i]
            # double beta = E0*thickness0*1.77245385/(1.0-nu0*nu0);
            beta = ((outlet_nodes_1d[0].GetValue(YOUNG_MODULUS) * outlet_nodes_1d[0].GetValue(THICKNESS) * math.sqrt(math.pi)) / (
                1 - (outlet_nodes_1d[0].GetValue(POISSON_RATIO) * outlet_nodes_1d[0].GetValue(POISSON_RATIO))))
            A = outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
            # press=outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
            # AQAUI HA TIRADO UN ERROR!!!!! en el paso 0,3380 (aproximadamente)
            press = initial_pressure + beta * \
                (math.sqrt(A) - math.sqrt(A0)) / \
                A0  # math.sqrt(A/A0)*beta - beta
            print "in Transfer PRESSURE 1D_to_3D :::::>>>>", press
            print "in Transfer Area 1D_to_3D :::::>>>>", A
            # Edu
            # press = beta/A0*(math.sqrt(A) - math.sqrt(A0))
            # print "in Transfer1D_to_3D(EDU)",press

            for node in outlet_nodes_3d:
                node.SetSolutionStepValue(PRESSURE, 0, press)
                node.SetSolutionStepValue(PRESSURE, 1, press)

    def Transfer3D_to_1D(self):
        for i in range(0, len(self.inlets_3d)):
            inlet_nodes_1d = self.inlets_1d[i]
            inlet_nodes_3d = self.inlets_3d[i]
            area3d = self.inlet_areas_3d[i]

            press_3d = 0.0
            counter = 0.0
            for node in inlet_nodes_3d:
                press_3d += node.GetSolutionStepValue(PRESSURE)
                counter += 1.0

            avg_press = press_3d / counter
            print "Trasfer inlet average pressure outlet 3D to 1D :::::>>>>> ", avg_press

            # TODO: make it to read from the input
            # print "inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS)",inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS)
            # print "inlet_nodes_1d[0].GetSolutionStepValue(THICKNESS)",inlet_nodes_1d[0].GetSolutionStepValue(THICKNESS)
            # print inlet_nodes_1d[0].Id
            beta = inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * inlet_nodes_1d[
                0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)
            A0 = inlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            # print "beta,", beta
            A = (avg_press * A0 / beta + math.sqrt(A0)
                 ) ** 2  # A0*(avg_press/beta + 1)**2
            print "in Transfer Area 3D_to_1D", A
            inlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA, 0, A)
            # print "inlet A (on node 23) ",A
            # print
            # "self.model_part_1d.Nodes[11].GetSolutionStepValue(NODAL_AREA) "
            # ,self.model_part_1d.Nodes[11].GetSolutionStepValue(NODAL_AREA)

            # assign flow to the outlet --> just for check
            flow = 0.0
            for node in inlet_nodes_3d:
                normal = node.GetSolutionStepValue(NORMAL)
                vel = node.GetSolutionStepValue(VELOCITY)
                flow += normal[0] * vel[0] + normal[
                    1] * vel[1] + normal[2] * vel[2]

            # self.flow_3d_in = str(flow)

            print "TRASFER to 3D to 1D:::flow entering= ", flow

            # self.f1d.write("FLOW3D -->1D (inlet_1D)--->   ")
            # self.f1d.write(str(flow))
            # self.f1d.write("\n")

            flow_aux = flow

        for i in range(0, len(self.outlets_3d)):
            outlet_nodes_1d = self.outlets_1d[i]
            outlet_nodes_3d = self.outlets_3d[i]
            area3d = self.outlet_areas_3d[i]

            # assign flow to the outlet
            flow = 0.0
            for node in outlet_nodes_3d:
                normal = node.GetSolutionStepValue(NORMAL)
                vel = node.GetSolutionStepValue(VELOCITY)
                flow += normal[0] * vel[0] + normal[
                    1] * vel[1] + normal[2] * vel[2]

            self.flow_3d_out = str(flow)
            self.flow_1d_in2 = str(flow)

            print "TRASFER to 3D to 1D:::flow exiting in 3D = ", flow

            # self.f1d.write(self.flow_3d_in+ " "+ self.flow_3d_out + "\n")

            # perdida= 0.0
            # if(flow_aux != 0):
                # perdida= ((flow_aux - flow) / flow_aux) *100
                # self.f.write(str(perdida))
                # self.f.write("\n")

            # if(perdida > 30):
                # self.f.write("--->")
                # self.f.write(str(perdida))
                # self.f.write("\n")
                # break

            # flow=2*flow-outlet_nodes_1d[0].GetSolutionStepValue(FLOW, 1)
            # print "con Riemman",flow

            outlet_nodes_1d[0].SetSolutionStepValue(FLOW, 0, flow)

#----------------------------------------------------------------------------------------------------------------------------------------
# Setting Contitions 3d
#-------------------------------------------------------------------------

    def Setting3d(self):

        print self
        print len(self.model_part_3d.Conditions)

        for node in self.model_part_3d.Nodes:
            node.Free(VELOCITY_X)
            node.Free(VELOCITY_Y)
            node.Free(VELOCITY_Z)
            node.Free(PRESSURE)
            node.SetSolutionStepValue(VISCOSITY, 0, 0.0035 / 1060.0)
            node.SetSolutionStepValue(DENSITY, 0, 1060.0)

        for cond in self.model_part_3d.Conditions:
            if(cond.Properties.Id == 100):  # inlet
                for node in cond.GetNodes():
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    node.Fix(VELOCITY_Z)
                    node.SetSolutionStepValue(FLAG_VARIABLE, 0, 100.0)
            if(cond.Properties.Id > 100):  # outlet
                for node in cond.GetNodes():
                    node.Fix(PRESSURE)
                    node.SetSolutionStepValue(FLAG_VARIABLE, 0, 1001.0)

        for cond in self.model_part_3d.Conditions:
            if(cond.Properties.Id == 1):  # sides --> note that this is done in an outer separated loop!!
                for node in cond.GetNodes():
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    node.Fix(VELOCITY_Z)
                    node.SetSolutionStepValue(FLAG_VARIABLE, 0, 0.0)

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
            # if(node.IsFixed(VELOCITY_X) == True and node.GetSolutionStepValue(VELOCITY_X) > 0.0001):
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

        print "n pressure nodes ", counter
