import math 
import cmath
import mpmath
import matplotlib.pyplot as plt

def prod(my_list, scalar):    
    my_list[:] = [scalar * x for x in my_list] 

class taylor_green_vector_field():
    def __init__(self):
        self.U0 = 1.0
        self.k = 10

    def Calculate(self, x, y, z):      
        U0 = self.U0
        k = self.k    
        
        # non-dimensionalizing input
        z = (z * k) % 10 # There are 10 identical cells, one on top of the other
        x = x * k 
        sx = math.sin(x)
        cx = math.cos(x) 
        sz = math.sin(z)
        cz = math.cos(z) 
        self.f =   [  sx * cz, 0.0, -cx * sz]
        self.fx =  [  cx * cz, 0.0,  sx * sz]
        self.fy =  [0.0, 0.0, 0.0]
        self.fz =  [ -sx * sz, 0.0,  cx * cz]
        self.fxy = [0.0, 0.0, 0.0]        
        self.fxz = [ -cx * sz, 0.0,  sx * cz]
        self.fyz = [0.0, 0.0, 0.0]
        self.fxx = [- sx * cz, 0.0,  cx * sz]
        self.fyy = [0.0, 0.0, 0.0]
        self.fzz = [- sx * cz, 0.0, -cx * sz]
        
        # dimensionalizing output
        U0_inv = 1. / U0
        k_inv = 1. / k
        prod(self.f, U0_inv)
        prod(self.fx, U0_inv * k_inv)
        prod(self.fz, U0_inv * k_inv)
        prod(self.fxz, U0_inv * k_inv ** 2)
        prod(self.fxx, U0_inv * k_inv ** 2)
        prod(self.fzz, U0_inv * k_inv ** 2)       
        self.fyx = self.fxy
        self.fzy = self.fyz
        self.fzx = self.fxz
