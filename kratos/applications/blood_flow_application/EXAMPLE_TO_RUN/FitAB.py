#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import math


class Fitter:

    def __init__(self, node_id):
        self.pins = []
        self.pouts = []
        self.Qs = []
        self.node_id = node_id
        self.Totaltimes = []

    def AddPin(self, pin):
        self.pins.append(pin)

    def AddPout(self, pout):
        self.pouts.append(pout)

    def AddQ(self, Q):
        self.Qs.append(Q)

    def AddTotal_time(self, total_time):
        self.Totaltimes.append(total_time)

    def DoFitting_1D(self, k, ffit_1d):
        print("DoFitting_1D")
        a11 = 0.0
        a12 = 0.0
        a22 = 0.0
        b1 = 0.0
        b2 = 0.0
        i = 0
          # print self.pins
        # print self.pouts

        for i in range(0, len(self.pins)):
            total_time = self.Totaltimes[i]
            node_id = self.node_id
            pout = self.pouts[i]
            pin = self.pins[i]
            dpi = (self.pins[i] - self.pouts[i])
            Qi = self.Qs[i]
            a11 += Qi ** 2
            a12 += (math.fabs(Qi)) * (Qi ** 2)
            a22 += Qi ** 4
            b1 += dpi * Qi
            b2 += dpi * (math.fabs(Qi)) * Qi
            ToWrite = str(node_id) + "	 "
            ToWrite += str(total_time) + "	 "
            ToWrite += str(pin) + "		"
            ToWrite += str(pout) + "		"
            ToWrite += str(dpi) + "		"
            ToWrite += str(Qi) + "\n"
            ffit_1d[k].write(ToWrite)
        a21 = a12

        # print "a11=",a11
        # print "a12=",a12
        # print "a22=",a22
        # print "b1 = ",b1
        # print "b2 = ",b2
        B = (b2 - a21 / a11 * b1) / (a22 - a12 * a21 / a11)
        A = (b1 / a11) - ((a12 * B) / a11)
        ToWrite = str(A) + " " + str(B) + "\n"
        ffit_1d[k].write(ToWrite)
        return [self.node_id, A, B]

    def DoFitting_3D_2(self, k, ffit_3d):
        print("DoFitting_3D_2: (Q positivo)")
        a11 = 0.0
        a12 = 0.0
        a22 = 0.0
        b1 = 0.0
        b2 = 0.0
        # print self.pins
        # print self.pouts
        # print range(0,len(self.pins))
        begin = int(0.2 * len(self.pins))
        for i in range(begin, len(self.pins)):
            # print i
            total_time = self.Totaltimes[i]
            node_id = self.node_id
            dpi = (self.pins[i] - self.pouts[i])
            Qi = self.Qs[i]
            pout = self.pouts[i]
            pin = self.pins[i]
            a11 += Qi ** 2
            a12 += (math.fabs(Qi)) * (Qi ** 2)
            a22 += Qi ** 4
            b1 += dpi * Qi
            b2 += dpi * (math.fabs(Qi)) * Qi
            ToWrite = str(node_id) + "	 "
            ToWrite += str(total_time) + "	 "
            ToWrite += str(pin) + "		"
            ToWrite += str(pout) + "		"
            ToWrite += str(dpi) + "		"
            ToWrite += str(Qi) + "\n"
            ffit_3d[k].write(ToWrite)
        a21 = a12
        B = (b2 - a21 / a11 * b1) / (a22 - a12 * a21 / a11)
        A = b1 / a11 - a12 / a11 * B
        ToWrite = str(A) + "\n "
        ffit_3d[k].write(ToWrite)
        ToWrite = str(B) + "\n "
        ffit_3d[k].write(ToWrite)
        # A=0
        # B=0
        return [self.node_id, A, B]

    def DoFitting_3D(self, k, ffit_3d):
        print("DoFitting_3D: mod(qi), pout-pi")
        a11 = 0.0
        a12 = 0.0
        a22 = 0.0
        b1 = 0.0
        b2 = 0.0
        # print self.pins
        # print self.pouts
        # print range(0,len(self.pins))
        begin = int(0.2 * len(self.pins))
        for i in range(begin, len(self.pins)):
            # print i
            total_time = self.Totaltimes[i]
            node_id = self.node_id
            dpi = (self.pins[i] - self.pouts[i])
            Qi = (-1) * self.Qs[i]
            pout = self.pouts[i]
            pin = self.pins[i]
            a11 += Qi ** 2
            a12 += (math.fabs(Qi)) * (Qi ** 2)
            a22 += Qi ** 4
            b1 += dpi * Qi
            b2 += dpi * (math.fabs(Qi)) * Qi
            ToWrite = str(node_id) + "	 "
            ToWrite += str(total_time) + "	 "
            ToWrite += str(pin) + "		"
            ToWrite += str(pout) + "		"
            ToWrite += str(dpi) + "		"
            ToWrite += str(Qi) + "\n"
            ffit_3d[k].write(ToWrite)
        a21 = a12
        B = (b2 - a21 / a11 * b1) / (a22 - a12 * a21 / a11)
        A = b1 / a11 - a12 / a11 * B
        ToWrite = str(A) + "\n "
        ffit_3d[k].write(ToWrite)
        ToWrite = str(B) + "\n "
        ffit_3d[k].write(ToWrite)
        # A=0
        # B=0
        return [self.node_id, A, B]

    # def DoFitting_3D_3(self,k,ffit_3d):
        #a11 = 0.0
        #a12 = 0.0
        #a22 = 0.0
        # b1=0.0
        #b2 = 0.0
        # print "DoFitting_3D_3: (fabs(pou-pi))and fabs(q)"
        # print self.pins
        # print self.pouts
        # print range(0,len(self.pins))
        # for i in range(0,len(self.pins)):
            # print i
            # total_time=self.Totaltimes[i]
            #node_id = self.node_id
            #dpi = math.fabs(self.pouts[i] - self.pins[i])
            #Qi = math.fabs(self.Qs[i])
            # pout=self.pouts[i]
            # pin=self.pins[i]
            #a11 += Qi**2
            #a12 += Qi**3
            #a22 += Qi**4
            #b1 += dpi*Qi
            #b2 += dpi*Qi*Qi
            #ToWrite = str(node_id) + "	 "
            #ToWrite += str(total_time) + "	 "
            #ToWrite += str(pin) + "		"
            #ToWrite += str(pout) +  "		"
            #ToWrite += str(dpi) + "		"
            #ToWrite += str(Qi) + "\n"
            # ffit_3d[k].write(ToWrite)
        #a21 = a12
        # B=(b2-a21/a11*b1)/(a22-a12*a21/a11)
        # A=b1/a11-a12/a11*B
        # print A
        # print B
        # print "--------------------------"
        #ToWrite = str(A) + "\n "
        # ffit_3d[k].write(ToWrite)
        #ToWrite = str(B) + "\n "
        # ffit_3d[k].write(ToWrite)
        # A=(b1-(b2*(a12/a22)))*(a22/((a11*a22)-(a12*a12)))
        # B=(b2/a22)-((a21*A)/a22)
        # print A
        # print B
        # raw_input()
        # A=0
        # B=0
        # return [self.node_id,A,B]

    # def DoFitting_3D_4(self,k,ffit_3d):
        #a11 = 0.0
        #a12 = 0.0
        #a22 = 0.0
        # b1=0.0
        #b2 = 0.0
        # print "DoFitting_3D_4: mod(qi), fabs(pi-pou)"
        # print self.pins
        # print self.pouts
        # print range(0,len(self.pins))
        # for i in range(0,len(self.pins)):
            # print i
            # total_time=self.Totaltimes[i]
            #node_id = self.node_id
            ##dpi = math.fabs(self.pouts[i] - self.pins[i])
            #dpi = self.pouts[i]-self.pins[i]
            #Qi = self.Qs[i]
            # pout=self.pouts[i]
            # pin=self.pins[i]
            #a11 += Qi**2
            #a12 += Qi**3
            #a22 += Qi**4
            #b1 += dpi*Qi
            #b2 += dpi*Qi**2
            #ToWrite = str(node_id) + "	 "
            #ToWrite += str(total_time) + "	 "
            #ToWrite += str(pin) + "		"
            #ToWrite += str(pout) +  "		"
            #ToWrite += str(dpi) + "		"
            #ToWrite += str(Qi) + "\n"
            # ffit_3d[k].write(ToWrite)
        #a21 = a12
        # B=(b2-a21/a11*b1)/(a22-a12*a21/a11)
        # A=b1/a11-a12/a11*B
        #ToWrite = str(A) + "\n "
        # ffit_3d[k].write(ToWrite)
        #ToWrite = str(B) + "\n "
        # ffit_3d[k].write(ToWrite)
        # A=0
        # B=0
        # return [self.node_id,A,B]

    # def DoFitting_3D_5(self,k,ffit_3d,Aprime,Bprime):
        #a11 = 0.0
        #a12 = 0.0
        #a22 = 0.0
        # b1=0.0
        #b2 = 0.0
        # print "DoFitting_3D_5:"
        # print self.pins
        # print self.pouts
        # print range(0,len(self.pins))
        # for i in range(0,len(self.pins)):
            # print i
            # total_time=self.Totaltimes[i]
            #node_id = self.node_id
            ##dpi = math.fabs(self.pouts[i] - self.pins[i])
            #dpi = self.pouts[i]-self.pins[i]
            #Qi = self.Qs[i]
            # pout=self.pouts[i]
            # pin=self.pins[i]
            #a11 += Qi**2
            #a12 += Qi**3
            #a22 += Qi**4
            #b1 += dpi*Qi
            #b2 += dpi*Qi*Qi
            #ToWrite  = str(node_id) + "	 "
            #ToWrite += str(total_time) + "	 "
            #ToWrite += str(pin) + "		"
            #ToWrite += str(pout) +  "		"
            #ToWrite += str(dpi) + "		"
            #ToWrite += str(Qi) + "\n"
            # ffit_3d[k].write(ToWrite)
        #a21 = a12
        # B=(b2-a21/a11*b1)/(a22-a12*a21/a11)
        # B=B/Bprime
        # A=((b1/a11-a12/a11*B)/Aprime)
        #ToWrite = str(A) + "\n "
        # ffit_3d[k].write(ToWrite)
        #ToWrite = str(B) + "\n "
        # ffit_3d[k].write(ToWrite)
        #ToWrite = str(Aprime) + "\n "
        # ffit_3d[k].write(ToWrite)
        #ToWrite = str(Bprime) + "\n "
        # ffit_3d[k].write(ToWrite)
        # A=0
        # B=0
        # return [self.node_id,A,B]
    def DoFitting_3D_Gauss(self, k, ffit_3d):
        a11 = 0.0
        a12 = 0.0
        a22 = 0.0
        a13 = 0.0
        a33 = 0.0
        a23 = 0.0
        a32 = 0.0
        b1 = 0.0
        b2 = 0.0
        b3 = 0.0
        Qmean = 0.0
        counter = 0.0
        A_Matriz = [[0 for x in range(3)] for x in range(3)]
        B_Matriz = [[0 for x in range(1)] for x in range(3)]

        print("DoFitting_3D_Gauss:")
        # print self.pins
        # print self.pouts
        # print range(0,len(self.pins))
        # Calculo Q medio
        for i in range(0, len(self.pins)):
            Qmean += self.Qs[i]
            counter += 1

        Qmean = Qmean / counter

        for i in range(0, len(self.pins)):
            # print i
            total_time = self.Totaltimes[i]
            node_id = self.node_id
            dpi = (self.pouts[i] - self.pins[i])
            Qi = self.Qs[i]
            pout = self.pouts[i]
            pin = self.pins[i]
            a11 += Qi ** 2
            a12 += (math.fabs(Qi)) * (Qi ** 2)
            a13 += Qmean * Qi
            a23 += Qmean * Qi * (math.fabs(Qi))
            a22 += Qi ** 4
            a33 += Qmean * Qmean
            b1 += dpi * Qi
            b2 += dpi * math.fabs(Qi) * Qi
            b3 += dpi * Qmean
            ToWrite = str(node_id) + "	 "
            ToWrite += str(total_time) + "	 "
            ToWrite += str(pin) + "		"
            ToWrite += str(pout) + "		"
            ToWrite += str(dpi) + "		"
            ToWrite += str(Qi) + "\n"
            ffit_3d[k].write(ToWrite)
        a21 = a12
        a32 = a23
        a31 = a13
        A_Matriz[0][0] = a11
        A_Matriz[0][1] = a12
        A_Matriz[0][2] = a13
        A_Matriz[1][0] = a21
        A_Matriz[1][1] = a22
        A_Matriz[1][2] = a23
        A_Matriz[2][0] = a31
        A_Matriz[2][1] = a32
        A_Matriz[2][2] = a33
        B_Matriz[0] = b1
        B_Matriz[1] = b2
        B_Matriz[2] = b3
        print(Qmean)
        print(A_Matriz)
        print(B_Matriz)
        [x] = doMatrixSolve(A_Matriz, B_Matriz)
        print(x)
        # B=(b2-a21/a11*b1)/(a22-a12*a21/a11)
        # B=B/BprimedoMatrixSolve
        # A=((b1/a11-a12/a11*B)/Aprime)
        #ToWrite = str(A) + "\n "
        # ffit_3d[k].write(ToWrite)
        #ToWrite = str(B) + "\n "
        # ffit_3d[k].write(ToWrite)
        #ToWrite = str(Aprime) + "\n "
        # ffit_3d[k].write(ToWrite)
        #ToWrite = str(Bprime) + "\n "
        # ffit_3d[k].write(ToWrite)
        # A=0
        # B=0
        return [self.node_id, x]


# A must be square n by n matrix, B is 1 by n vector
def doMatrixSolve(A, B, tol=10 ** (-7)):
    x_new = [[1. for i in range(len(B[0]))]]
    converged = False
    while not converged:
        x_old = [x_new[0][:]]
        x_new = gaussSiedel(A, B, x_old)

        converged = True  # succeeds by default
        for i in range(len(x_old[0])):
            if abs(2 * (x_new[0][i] - x_old[0][i]) / (x_new[0][i] + x_old[0][i])) >= tol:
                converged = False
                break
    return x_new


def gaussSiedel(A, B, x):
    n = len(x[0])
    y = [[0. for l in range(n)]]
    for i in range(len(x[0])):
        val = 0  # dummy variable
        for j in range(i):
            val += A[i][j] * y[0][j]
            # print 'a',i,j
        for j in range(i + 1, n):
            val += A[i][j] * x[0][j]
            # print 'b',i,j
        y[0][i] = (B[i] - val) / A[i][i]
    # print y
    return y
