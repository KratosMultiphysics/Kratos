import math

def FindRootBisection(function, x0, tol = 1e-12):
    x = x0
    dx = 2 * tol + 1
    residual = 2 * tol
    
    while abs(residual) > tol:
        residual = function.f(x)
        x -= residual / function.df(x)
    return x

class pol:
    def __init__(self):
        pass
    def SetObjective(self, x):
        self.x = x
    def f(self, y):
        return y ** 2 - self.x
    def df(self, y):
        return 2 * y

my_pol = pol()
my_pol.SetObjective(3)
x = FindRootBisection(my_pol, 1)
print('x = ', x)