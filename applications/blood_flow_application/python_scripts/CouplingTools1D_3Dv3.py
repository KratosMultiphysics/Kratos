from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
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
        # self.fin = open("areasIN.txt","w")
        # self.fout = open("areasOUT.txt","w")
        # self.fout2 = open("areasOUT_2.txt","w")
        # self.press=open("press.txt","w")
        # self.flow=open("flow.txt","w")
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
        # self.fin.write("EMPIEZO\n")
        # self.fout.write("EMPIEZO\n")
        # self.fout2.write("EMPIEZO\n")
        # self.press.write("EMPIEZO\n")
        # self.flow.write("EMPIEZO|n")
        # self.f3d.write("EMPIEZO\n")
        # inicial = 1
        # compute normals and 3d areas
        BodyNormalCalculationUtils().CalculateBodyNormals(
            self.model_part_3d, 3)

        # detect 1d inlets
        for i in range(100, 101):
            nfound = 0
            aux = []
            for node in self.model_part_1d.Nodes:
                if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
                    aux.append(node)
                    node.Fix(NODAL_AREA)
                    node.Fix(FLOW)
                    print("1D_inlet has been assigned", node)
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
                    print("3D_inlet has been assigned", node)
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
                print("1D_outlet has been assigned(H)", node)
                # raw_input()
            else:
                break

        # raw_input()
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
        # raw_input()

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
            # print "area3d_outlet", area3d

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
        # for node in self.inlets_1d[0]:
            # print node.Id
            # raw_input
        # for i in range(0,len(self.inlets_1d)):
            # inlet_nodes_1d = self.inlets_1d[i]
            # for i in range(0,len(self.inlets_3d)):
            # inlet_nodes_3d = self.inlets_3d[i]
            # area3d = self.inlet_areas_3d[i]
            # directions = self.inlet_velocity_directions[i]
            # radio3d=math.sqrt(area3d*3.1416)
            # vel1d =  inlet_nodes_1d[0].GetSolutionStepValue(FLOW) / area3d
            # k = 0
            # for node in self.model_part_3d.Nodes:
                # node.SetSolutionStepValue(VELOCITY, 0, vel1d)
        # raw_input()

    def Initial_Contitions(self):

        initial_pressure = 0
        print("Inicializo 3D")
        for i in range(0, len(self.inlets_1d)):
            inlet_nodes_1d = self.inlets_1d[i]
            print("3D-1D: inlet_nodes_1d [0].Id::::::>>>> ", inlet_nodes_1d[0].Id)
            for i in range(0, len(self.outlets_1d)):
                outlet_nodes_1d = self.outlets_1d[i]
                pressinlet3D = outlet_nodes_1d[
                    0].GetSolutionStepValue(PRESSURE)
                # print pressinlet3D
                # print "3D-1D: outlet_nodes_1d[0].Id::::::>>>>
                # ",outlet_nodes_1d[0].Id
                for i in range(0, len(self.inlets_3d)):
                    inlet_nodes_3d = self.inlets_3d[i]
                    area3d = self.inlet_areas_3d[i]
                    directions = self.inlet_velocity_directions[i]
                    radio3d = math.sqrt(area3d * 3.1416)
                    vel1d = inlet_nodes_1d[
                        0].GetSolutionStepValue(FLOW) / area3d
                    k = 0
                    # print vel1d
                    # Impongo velocidad y la presion
                    for node in inlet_nodes_3d:
                        n = node.GetSolutionStepValue(NORMAL)
                        a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                        orientation = directions[k]
                        node.SetSolutionStepValue(
                            VELOCITY, 0, directions[k] * vel1d)
                        # node.SetSolutionStepValue(PRESSURE, 0, pressinlet3D)
                        k = k + 1
            input()

    def Transfer1D_to_3D(self):

            # ARCHIVE TO SET :::::::::::::::::::::::::::>>>>>>>>>>>>>> VARIABLES
            # import config_full
        initial_pressure = config.systolic_pressure
        initial_pressure = 0

        print("Transfer1D_to_3D")
        for i in range(0, len(self.inlets_1d)):
            inlet_nodes_1d = self.inlets_1d[i]
            # print "NODO 1D-3D:: inlet_nodes_1d[0].Id::::::>>>> ",inlet_nodes_1d[0].Id
            # print "--"
            # raw_input()
            for i in range(0, len(self.inlets_3d)):
                inlet_nodes_3d = self.inlets_3d[i]
                area3d = self.inlet_areas_3d[i]
                directions = self.inlet_velocity_directions[i]
                radio3d = math.sqrt(area3d * 3.1416)
                # print "radio-->3d-->",radio3d
                # print "area-->3d_inlet",area3d
                # print "area 1d-->", inlet_nodes_1d[0].GetValue(NODAL_AREA)
                # print "diferencia de Areas en el inlet 1d--3d:::::", inlet_nodes_1d[0].GetValue(NODAL_AREA)-area3d
                # print "Q", inlet_nodes_1d[0].GetSolutionStepValue(FLOW)
                # print "velocity1D", inlet_nodes_1d[0].GetSolutionStepValue(FLOW)/inlet_nodes_1d[0].GetValue(NODAL_AREA)
                # self.fin.write("Areas inlets ")
                # self.fin.write(str(inlet_nodes_1d[0].GetValue(NODAL_AREA)))
                # self.fin.write(" ")
                # self.fin.write(str(inlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)))
                # self.fin.write(" ")
                # self.fin.write(str(area3d))
                # self.fin.write(" ")
                # self.fin.write(str(inlet_nodes_1d[0].GetValue(NODAL_AREA)-area3d))
                # self.fin.write("\n")
                # self.press.write(str(outlet_nodes_1d[0].Id)+"-->  ")
                # self.press.write(str(press))
                # self.press.write("\n")
                # self.f1d.write(str(inlet_nodes_1d[0].GetSolutionStepValue(FLOW)))
                # self.f1d.write("\n")
                # self.flow_1d =
                # str(inlet_nodes_1d[0].GetSolutionStepValue(FLOW))
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
                    # velocidad=directions[k]*vel1d*(a/area3d)
                    node.SetSolutionStepValue(
                        VELOCITY, 0, directions[k] * vel1d)
                    # Edu
                    # node.SetSolutionStepValue( VELOCITY, 0, velocidad)
                    # node.SetSolutionStepValue( PRESSURE, 0, pressinlet1D)
                    # print node.Id
                    # print directions[k]*vel1d
                    # print velocidad
                    k = k + 1

                # raw_input()
                # print "velocity  1d to 3D ::::::>", vel1d
                # print "flow             1d to 3D ::::::>", inlet_nodes_1d[0].GetSolutionStepValue(FLOW)
                # print "pressure  1d to 3D ::::::> ", pressinlet1D
            # raw_input()
            # print "Acabo de inlets"
               # raw_input()
        # print "EMPIEZO de outlets"
        # raw_input()
        for i in range(0, len(self.outlets_1d)):
            # print "HOLA!"
            # raw_input()
            outlet_nodes_1d = self.outlets_1d[i]
            pressinlet3D = outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
            # pressinlet3D = 0
            # print "3D-1D: outlet_nodes_1d[0].Id::::::>>>> ",outlet_nodes_1d[0].Id
            # print "pressinlet3D ",pressinlet3D
            # raw_input()
            # double beta = E0*thickness0*1.77245385/(1.0-nu0*nu0);
            beta = ((outlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * outlet_nodes_1d[0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)) / (
                1 - (outlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO) * outlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO))))
            A = outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
            # press=outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
            press = initial_pressure + beta * \
                (math.sqrt(A) - math.sqrt(A0)) / \
                A0  # math.sqrt(A/A0)*beta - beta
            # print "in Transfer PRESSURE 1D_to_3D :::::>>>>", press
            # print "in Transfer Area 1D_to_3D(A) :::::>>>>", A
            # print "in Transfer Area 1D_to_3D(A0) :::::>>>>", A0
            # print "--->", (math.sqrt(A)-math.sqrt(A0))/A0
            # print "--beta", beta
            # print "", outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            # print "", outlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS)
            # print "", outlet_nodes_1d[0].GetSolutionStepValue(THICKNESS)
            # print "", outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
            # print "in Transfer1D_to_3D(EDU)",press
            # print "beta", beta
            # raw_input()
            # print "press ",press
            # press = pressinlet3D
            # press = -10.0
            # raw_input()
            # print "in Transfer Area 1D_to_3D :::::>>>>", A
            # Edu
            # press = beta/A0*(math.sqrt(A) - math.sqrt(A0))
            # print "--trasferencia al 3D"
            # self.press.write(str(outlet_nodes_1d[0].Id)+"-->  ")
            # self.press.write(str(press))
            # self.press.write("\n")
            # raw_input()
            for i in range(0, len(self.outlets_3d)):
                outlet_nodes_3d = self.outlets_3d[i]
                # print "outlet_nodes_1d", outlet_nodes_1d[0].Id
                area3d = self.outlet_areas_3d[i]
                # print "Area1D--->", outlet_nodes_1d[0].GetValue(NODAL_AREA)
                # print "Area3d--->", area3d
                # inlet_nodes_1d = self.outlets_1d[i]
                # print "1d-3d Areas Diference--->", outlet_nodes_1d[0].GetValue(NODAL_AREA)-area3d
                # self.fout.write(str(outlet_nodes_1d[0].GetValue(NODAL_AREA)))
                # self.fout.write(" ")
                # self.fout.write(str(outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)))
                # self.fout.write(" ")
                # self.fout.write(str(area3d))
                # self.fout.write(" ")
                # self.fout.write(str(outlet_nodes_1d[0].GetValue(NODAL_AREA)-area3d))
                # self.fout.write("\n")
                # for node in outlet_nodes_3d:
                # print node.Id
                for node in outlet_nodes_3d:
                    node.SetSolutionStepValue(PRESSURE, 0, pressinlet3D)
                    # print "pressure  1d to 3D ::::::> ", pressinlet3D
                # ModificacionEdu
                # for node in outlet_nodes_3d:
                    # old_pressure = node.GetSolutionStepValue(PRESSURE,1)
                    # effective_press = 0.5*(press+old_pressure)
                    # node.SetSolutionStepValue( PRESSURE, 0, effective_press)
                    # node.SetSolutionStepValue( PRESSURE, 0, press)
                    # print "node", node.Id
                    # node.SetSolutionStepValue( PRESSURE, 1, press)
            # raw_input()
        # print "Acabo de outles"
        # raw_input()

    def Transfer3D_to_1D(self):

        for i in range(0, len(self.inlets_3d)):
            inlet_nodes_1d = self.inlets_1d[i]
            inlet_nodes_3d = self.inlets_3d[i]
            area3d = self.inlet_areas_3d[i]
            outlet_nodes_1d = self.outlets_1d[i]
            press_3d = 0.0
            counter = 0.0
            for node in inlet_nodes_3d:
                press_3d += node.GetSolutionStepValue(PRESSURE)
                counter += 1.0
            avg_press = press_3d / counter

            # print "Trasfer inlet average PRESSURE outlet 3D to 1D :::::>>>>> ",avg_press
            # print "Inlet_nodes_1d[0].Id::::::>>>> ",inlet_nodes_1d[0].Id
            # inlet_nodes_1d[0].SetSolutionStepValue(PRESSURE ,0, avg_press)
            # raw_input()

            # TODO: make it to read from the input

            # print "inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS)",inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS)
            # print "inlet_nodes_1d[0].GetSolutionStepValue(THICKNESS)",inlet_nodes_1d[0].GetSolutionStepValue(THICKNESS)
            # print inlet_nodes_1d[0].Id
            beta = inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * inlet_nodes_1d[
                0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)
            A0 = inlet_nodes_1d[0].GetValue(NODAL_AREA)
            # print beta
            # beta = model_part1D.Nodes[63].GetSolutionStepValue(YOUNG_MODULUS) * model_part1D.Nodes[63].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)
            # A0 = model_part1D.Nodes[63].GetValue(NODAL_AREA)
            # print beta
            # print A0
            # print "(in transfer3d->1d, A0--->",A0
            # print "beta,", beta
            A = (avg_press * A0 / beta + math.sqrt(A0)
                 ) ** 2  # A0*(avg_press/beta + 1)**2
            print("in Transfer AREA 3D_to_1D----->", A)
            # self.fout2.write(str(inlet_nodes_1d[0].GetValue(NODAL_AREA)))
            # self.fout2.write(" ")
            # self.fout2.write(str(inlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)))
            # self.fout2.write(" ")
            # self.fout2.write(str(area3d))
            # self.fout2.write(" ")
            # self.fout2.write(str(inlet_nodes_1d[0].GetValue(NODAL_AREA)-area3d))
            # self.fout2.write("\n")
            # self.press.write(str(inlet_nodes_1d[0].Id)+"-->  ")
            # self.press.write(str(avg_press))
            # self.press.write("\n")
            A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
            outlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA, 0, A0)
            inlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA, 0, A0)
            # Ultimo cambio Edu
            # inlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA ,0, area3d)
            # print "inlet A (on node 23) ",A
            # print "self.model_part_1d.Nodes[11].GetSolutionStepValue(NODAL_AREA) " ,self.model_part_1d.Nodes[11].GetSolutionStepValue(NODAL_AREA)
            # raw_input()
            # assign flow to the outlet --> just for check
            flow = 0.0
            for node in inlet_nodes_3d:
                normal = node.GetSolutionStepValue(NORMAL)
                vel = node.GetSolutionStepValue(VELOCITY)
                flow += normal[0] * vel[0] + normal[
                    1] * vel[1] + normal[2] * vel[2]
            # self.flow_3d_in = str(flow)
            # print "TRASFER to 3D to 1D:::flow entering= ",flow
            # flow_aux=  flow
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

            # self.flow_3d_out = str(flow)
            # self.flow_1d_in2 = str(flow)
            # print "TRASFER to 3D to 1D:::flow exiting in 3D = ",flow
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
            # print "sin Riemman",-flow
            # print outlet_nodes_1d[0].GetSolutionStepValue(FLOW, 1)
            # print outlet_nodes_1d[0].Id
            # flow=2*flow-outlet_nodes_1d[0].GetSolutionStepValue(FLOW, 1)
            # print "con Riemman",-flow

        outlet_nodes_1d[0].SetSolutionStepValue(FLOW, 0, -flow)

                 # print outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)

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
