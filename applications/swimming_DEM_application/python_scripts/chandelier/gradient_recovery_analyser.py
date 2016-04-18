import math 
import cmath
import mpmath
import matplotlib.pyplot as plt

class sinus_field():
    def __init__(self):
        pass
    def Calculate(self, x, y, z):
        self.f = math.sin(x + y + z)
        self.fx = math.cos(x + y + z)
        self.fy = math.cos(x + y + z)
        self.fz = math.cos(x + y + z)
        self.fxx = -math.sin(x + y + z)
        self.fyy = -math.sin(x + y + z)
        self.fzz = -math.sin(x + y + z)
        
class sinus_vector_field():
    def __init__(self):
        pass
    def Calculate(self, x, y, z):
        self.f = [math.sin(y*z), math.sin(x*z), math.sin(x*y)]
        self.fx = [0, z * math.cos(x*z), y * math.cos(x*y)]
        self.fy = [z * math.cos(y*z), 0, x * math.cos(x*y)]
        self.fz = [y * math.cos(y*z), x * math.cos(x*z), 0]
        self.fxx = [0, - z ** 2 * math.sin(x*z), - y ** 2 * math.sin(x*y)]
        self.fyy = [- z ** 2 * math.sin(y*z), 0, - x ** 2 * math.sin(x*y)]
        self.fzz = [- y ** 2 * math.sin(y*z), - x ** 2 * math.sin(x*z), 0]                

class GradientRecoveryAnalyser:
    def __init__(self, scalar_field, vector_field):
        self.f = scalar_field
        self.v = vector_field
    
    def CalculateDerivatives(self, x, y, z):
        f = self.f
        v = self.v
        f.Calculate(x, y, z)
        v.Calculate(x, y, z)
        self.value = f.f
        self.gx = f.fx
        self.gy = f.fy
        self.gz = f.fz
        self.v0 = v.f[0]
        self.v1 = v.f[1] 
        self.v2 = v.f[2] 
        self.lx = v.fxx[0] + v.fyy[0] + v.fzz[0] 
        self.ly = v.fxx[1] + v.fyy[1] + v.fzz[1] 
        self.lz = v.fxx[2] + v.fyy[2] + v.fzz[2] 
        
