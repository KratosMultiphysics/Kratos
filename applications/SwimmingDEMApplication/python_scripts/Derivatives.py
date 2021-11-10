import numpy as np
import sympy as sp

class TimeDependantFluidFractionSolution():
    def __init__(self):
        self.time  = sp.Symbol('time')
        self.x1 = sp.Symbol('x1')
        self.x2 = sp.Symbol('x2')
        self.L  = sp.Symbol('L')
        self.delta_alpha = sp.Symbol('delta_alpha')
        self.omega = sp.Symbol('omega')
        self.squeeze_amplitude = sp.Symbol('squeeze_amplitude')
        self.x10 = sp.Symbol('x10')
        self.x20 = sp.Symbol('x20')
        self.R = sp.Symbol('R')
        self.u_char = sp.Symbol('u_char')
        self.c = (1 + self.squeeze_amplitude * sp.sin(self.omega * self.time))

    def VelocityField(self):

        self.alpha = 1 - self.delta_alpha * sp.exp(1 - 1 / (1 - ((self.c*(self.x1 - self.x10)/self.R)**2 + (((self.x2 - self.x20)/self.R)/self.c)**2)))
        self.f_x1 = self.u_char * self.x1**2 * (1 - self.x1)**2
        self.f_x2 = self.u_char * self.x2**2 * (1 - self.x2)**2
        self.g_t = sp.cos(sp.pi * self.time) * sp.exp(-self.time)
        self.df_x1 = sp.diff(self.f_x1, self.x1)
        self.df_x2 = sp.diff(self.f_x2, self.x2)
        self.u1_codina = self.f_x1 * self.df_x2 * self.g_t
        self.u2_codina = -self.df_x1 * self.f_x2 * self.g_t

        self.u1 = self.u1_codina / self.alpha
        self.u2 = self.u2_codina / self.alpha

    def DerivativeCalculator(self):
        self.VelocityField()
        print('u1 : \n', self.u1)
        print('u2 : \n', self.u2)
        self.du1dt = sp.diff(self.u1, self.time)
        print('du1dt : \n', self.du1dt)
        self.du2dt = sp.diff(self.u2, self.time)
        print('du2dt : \n', self.du2dt)
        print('alpha : \n', self.alpha)
        self.dalphat = sp.diff(self.alpha, self.time)
        print('dalphat : \n', self.dalphat)
        self.alpha1 = sp.diff(self.alpha, self.x1)
        print('alpha1 : \n', self.alpha1)
        self.alpha2 = sp.diff(self.alpha, self.x2)
        print('alpha2 : \n', self.alpha2)
        self.du11 = sp.diff(self.u1, self.x1)
        print('du11 : \n', self.du11)
        self.du12 = sp.diff(self.u1, self.x2)
        print('du12 : \n', self.du12)
        self.du111 = sp.diff(self.du11, self.x1)
        print('du111 : \n', self.du111)
        self.du112 = sp.diff(self.du11, self.x2)
        print('du112 : \n', self.du112)
        self.du121 = sp.diff(self.du12, self.x1)
        print('du121 : \n', self.du121)
        self.du122 = sp.diff(self.du12, self.x2)
        print('du122 : \n', self.du122)
        self.du21 = sp.diff(self.u2, self.x1)
        print('du21 : \n', self.du21)
        self.du22 = sp.diff(self.u2, self.x2)
        print('du22 : \n', self.du22)
        self.du211 = sp.diff(self.du21, self.x1)
        print('du211 : \n', self.du211)
        self.du212 = sp.diff(self.du21, self.x2)
        print('du212 : \n', self.du212)
        self.du221 = sp.diff(self.du22, self.x1)
        print('du221 : \n', self.du221)
        self.du222 = sp.diff(self.du22, self.x2)
        print('du222 : \n', self.du222)

class PorousCodina2007():
    def __init__(self):
        self.time  = sp.Symbol('time')
        self.x1 = sp.Symbol('x1')
        self.x2 = sp.Symbol('x2')
        self.L  = sp.Symbol('L')
        self.delta_alpha = sp.Symbol('delta_alpha')
        self.omega = sp.Symbol('omega')
        self.squeeze_amplitude = sp.Symbol('squeeze_amplitude')
        self.x10 = sp.Symbol('x10')
        self.x20 = sp.Symbol('x20')
        self.R = sp.Symbol('R')
        self.u_char = sp.Symbol('u_char')
        self.c = (1 + self.squeeze_amplitude * sp.sin(self.omega))

    def VelocityField(self):

        self.alpha = 1 - self.delta_alpha * sp.exp(1 - 1 / (1 - ((self.c*(self.x1 - self.x10)/self.R)**2 + (((self.x2 - self.x20)/self.R)/self.c)**2)))
        self.u1_codina = sp.sin(sp.pi * self.x1 - 0.7) * sp.sin(sp.pi * self.x2 + 0.2)
        self.u2_codina = sp.cos(sp.pi * self.x1 - 0.7) * sp.cos(sp.pi * self.x2 + 0.2)
        self.pressure = sp.sin(sp.pi * self.x1)*sp.cos(sp.pi*self.x2) + (sp.cos(1)-1)*sp.sin(1)

        self.u1 = self.u1_codina / self.alpha
        self.u2 = self.u2_codina / self.alpha

    def DerivativeCalculator(self):
        self.VelocityField()
        print('u1 : \n', self.u1)
        print('u2 : \n', self.u2)
        self.du1dt = sp.diff(self.u1, self.time)
        print('du1dt : \n', self.du1dt)
        self.du2dt = sp.diff(self.u2, self.time)
        print('du2dt : \n', self.du2dt)
        print('alpha : \n', self.alpha)
        self.dalphat = sp.diff(self.alpha, self.time)
        print('dalphat : \n', self.dalphat)
        self.alpha1 = sp.diff(self.alpha, self.x1)
        print('alpha1 : \n', self.alpha1)
        self.alpha2 = sp.diff(self.alpha, self.x2)
        print('alpha2 : \n', self.alpha2)
        self.du11 = sp.diff(self.u1, self.x1)
        print('du11 : \n', self.du11)
        self.du12 = sp.diff(self.u1, self.x2)
        print('du12 : \n', self.du12)
        self.du111 = sp.diff(self.du11, self.x1)
        print('du111 : \n', self.du111)
        self.du112 = sp.diff(self.du11, self.x2)
        print('du112 : \n', self.du112)
        self.du121 = sp.diff(self.du12, self.x1)
        print('du121 : \n', self.du121)
        self.du122 = sp.diff(self.du12, self.x2)
        print('du122 : \n', self.du122)
        self.du21 = sp.diff(self.u2, self.x1)
        print('du21 : \n', self.du21)
        self.du22 = sp.diff(self.u2, self.x2)
        print('du22 : \n', self.du22)
        self.du211 = sp.diff(self.du21, self.x1)
        print('du211 : \n', self.du211)
        self.du212 = sp.diff(self.du21, self.x2)
        print('du212 : \n', self.du212)
        self.du221 = sp.diff(self.du22, self.x1)
        print('du221 : \n', self.du221)
        self.du222 = sp.diff(self.du22, self.x2)
        print('du222 : \n', self.du222)
        print('pressure : \n', self.pressure)
        self.dp1 = sp.diff(self.pressure, self.x1)
        print('dp1 : \n', self.dp1)
        self.dp2 = sp.diff(self.pressure, self.x2)
        print('dp2 : \n', self.dp2)

class Codina2007():
    def __init__(self):
        self.time  = sp.Symbol('time')
        self.x1 = sp.Symbol('x1')
        self.x2 = sp.Symbol('x2')
        self.L  = sp.Symbol('L')
        self.delta_alpha = sp.Symbol('delta_alpha')
        self.omega = sp.Symbol('omega')
        self.squeeze_amplitude = sp.Symbol('squeeze_amplitude')
        self.x10 = sp.Symbol('x10')
        self.x20 = sp.Symbol('x20')
        self.R = sp.Symbol('R')
        self.u_char = sp.Symbol('u_char')
        self.c = (1 + self.squeeze_amplitude * sp.sin(self.omega))

    def VelocityField(self):

        self.alpha = 1
        self.u1_codina = sp.sin(sp.pi * self.x1 - 0.7) * sp.sin(sp.pi * self.x2 + 0.2)
        self.u2_codina = sp.cos(sp.pi * self.x1 - 0.7) * sp.cos(sp.pi * self.x2 + 0.2)
        self.pressure = sp.sin(sp.pi * self.x1)*sp.cos(sp.pi*self.x2) + (sp.cos(1)-1)*sp.sin(1)

        self.u1 = self.u1_codina / self.alpha
        self.u2 = self.u2_codina / self.alpha

    def DerivativeCalculator(self):
        self.VelocityField()
        print('u1 : \n', self.u1)
        print('u2 : \n', self.u2)
        self.du1dt = sp.diff(self.u1, self.time)
        print('du1dt : \n', self.du1dt)
        self.du2dt = sp.diff(self.u2, self.time)
        print('du2dt : \n', self.du2dt)
        print('alpha : \n', self.alpha)
        self.dalphat = sp.diff(self.alpha, self.time)
        print('dalphat : \n', self.dalphat)
        self.alpha1 = sp.diff(self.alpha, self.x1)
        print('alpha1 : \n', self.alpha1)
        self.alpha2 = sp.diff(self.alpha, self.x2)
        print('alpha2 : \n', self.alpha2)
        self.du11 = sp.diff(self.u1, self.x1)
        print('du11 : \n', self.du11)
        self.du12 = sp.diff(self.u1, self.x2)
        print('du12 : \n', self.du12)
        self.du111 = sp.diff(self.du11, self.x1)
        print('du111 : \n', self.du111)
        self.du112 = sp.diff(self.du11, self.x2)
        print('du112 : \n', self.du112)
        self.du121 = sp.diff(self.du12, self.x1)
        print('du121 : \n', self.du121)
        self.du122 = sp.diff(self.du12, self.x2)
        print('du122 : \n', self.du122)
        self.du21 = sp.diff(self.u2, self.x1)
        print('du21 : \n', self.du21)
        self.du22 = sp.diff(self.u2, self.x2)
        print('du22 : \n', self.du22)
        self.du211 = sp.diff(self.du21, self.x1)
        print('du211 : \n', self.du211)
        self.du212 = sp.diff(self.du21, self.x2)
        print('du212 : \n', self.du212)
        self.du221 = sp.diff(self.du22, self.x1)
        print('du221 : \n', self.du221)
        self.du222 = sp.diff(self.du22, self.x2)
        print('du222 : \n', self.du222)
        print('pressure : \n', self.pressure)
        self.dp1 = sp.diff(self.pressure, self.x1)
        print('dp1 : \n', self.dp1)
        self.dp2 = sp.diff(self.pressure, self.x2)
        print('dp2 : \n', self.dp2)

class CodinaVelocityField():
    def __init__(self):
        self.time  = sp.Symbol('time')
        self.x1 = sp.Symbol('x1')
        self.x2 = sp.Symbol('x2')
        self.u_char = sp.Symbol('u_char')

    def VelocityField(self):
        self.alpha = 1.0
        self.f_x1 = self.u_char * self.x1**2 * (1 - self.x1)**2
        self.f_x2 = self.u_char * self.x2**2 * (1 - self.x2)**2
        self.g_t = sp.cos(sp.pi * self.time) * sp.exp(-self.time)
        self.df_x1 = sp.diff(self.f_x1, self.x1)
        self.df_x2 = sp.diff(self.f_x2, self.x2)
        self.u1_codina = self.f_x1 * self.df_x2 * self.g_t
        self.u2_codina = -self.df_x1 * self.f_x2 * self.g_t

        self.u1 = self.u1_codina / self.alpha
        self.u2 = self.u2_codina / self.alpha

    def DerivativeCalculator(self):
        self.VelocityField()
        print('u1 : \n', self.u1)
        print('u2 : \n', self.u2)
        self.du1dt = sp.diff(self.u1, self.time)
        print('du1dt : \n', self.du1dt)
        self.du2dt = sp.diff(self.u2, self.time)
        print('du2dt : \n', self.du2dt)
        print('alpha : \n', self.alpha)
        self.dalphat = sp.diff(self.alpha, self.time)
        print('dalphat : \n', self.dalphat)
        self.alpha1 = sp.diff(self.alpha, self.x1)
        print('alpha1 : \n', self.alpha1)
        self.alpha2 = sp.diff(self.alpha, self.x2)
        print('alpha2 : \n', self.alpha2)
        self.du11 = sp.diff(self.u1, self.x1)
        print('du11 : \n', self.du11)
        self.du12 = sp.diff(self.u1, self.x2)
        print('du12 : \n', self.du12)
        self.du111 = sp.diff(self.du11, self.x1)
        print('du111 : \n', self.du111)
        self.du112 = sp.diff(self.du11, self.x2)
        print('du112 : \n', self.du112)
        self.du121 = sp.diff(self.du12, self.x1)
        print('du121 : \n', self.du121)
        self.du122 = sp.diff(self.du12, self.x2)
        print('du122 : \n', self.du122)
        self.du21 = sp.diff(self.u2, self.x1)
        print('du21 : \n', self.du21)
        self.du22 = sp.diff(self.u2, self.x2)
        print('du22 : \n', self.du22)
        self.du211 = sp.diff(self.du21, self.x1)
        print('du211 : \n', self.du211)
        self.du212 = sp.diff(self.du21, self.x2)
        print('du212 : \n', self.du212)
        self.du221 = sp.diff(self.du22, self.x1)
        print('du221 : \n', self.du221)
        self.du222 = sp.diff(self.du22, self.x2)
        print('du222 : \n', self.du222)

class PorousSinusoidalFluidFractionSolution():
    def __init__(self):
        self.time  = sp.Symbol('time')
        self.x1 = sp.Symbol('x1')
        self.x2 = sp.Symbol('x2')
        self.k = sp.Symbol('k')
        self.L  = sp.Symbol('L')
        self.delta_alpha = sp.Symbol('delta_alpha')
        self.x10 = sp.Symbol('x10')
        self.x20 = sp.Symbol('x20')
        self.u_char = sp.Symbol('u_char')

    def VelocityField(self):

        self.alpha = 1.0/2.0 + self.delta_alpha * sp.sin(self.k * self.x1) * sp.cos(self.k * self.x2)
        self.f_x1 = self.u_char * self.x1**2 * (1 - self.x1)**2
        self.f_x2 = self.u_char * self.x2**2 * (1 - self.x2)**2
        self.g_t = sp.cos(sp.pi * 1.0) * sp.exp(-1)
        self.df_x1 = sp.diff(self.f_x1, self.x1)
        self.df_x2 = sp.diff(self.f_x2, self.x2)
        self.u1_codina = self.f_x1 * self.df_x2 * self.g_t
        self.u2_codina = -self.df_x1 * self.f_x2 * self.g_t

        self.u1 = self.u1_codina / self.alpha
        self.u2 = self.u2_codina / self.alpha

    def DerivativeCalculator(self):
        self.VelocityField()
        print('u1 : \n', self.u1)
        print('u2 : \n', self.u2)
        self.du1dt = sp.diff(self.u1, self.time)
        print('du1dt : \n', self.du1dt)
        self.du2dt = sp.diff(self.u2, self.time)
        print('du2dt : \n', self.du2dt)
        print('alpha : \n', self.alpha)
        self.dalphat = sp.diff(self.alpha, self.time)
        print('dalphat : \n', self.dalphat)
        self.alpha1 = sp.diff(self.alpha, self.x1)
        print('alpha1 : \n', self.alpha1)
        self.alpha2 = sp.diff(self.alpha, self.x2)
        print('alpha2 : \n', self.alpha2)
        self.du11 = sp.diff(self.u1, self.x1)
        print('du11 : \n', self.du11)
        self.du12 = sp.diff(self.u1, self.x2)
        print('du12 : \n', self.du12)
        self.du111 = sp.diff(self.du11, self.x1)
        print('du111 : \n', self.du111)
        self.du112 = sp.diff(self.du11, self.x2)
        print('du112 : \n', self.du112)
        self.du121 = sp.diff(self.du12, self.x1)
        print('du121 : \n', self.du121)
        self.du122 = sp.diff(self.du12, self.x2)
        print('du122 : \n', self.du122)
        self.du21 = sp.diff(self.u2, self.x1)
        print('du21 : \n', self.du21)
        self.du22 = sp.diff(self.u2, self.x2)
        print('du22 : \n', self.du22)
        self.du211 = sp.diff(self.du21, self.x1)
        print('du211 : \n', self.du211)
        self.du212 = sp.diff(self.du21, self.x2)
        print('du212 : \n', self.du212)
        self.du221 = sp.diff(self.du22, self.x1)
        print('du221 : \n', self.du221)
        self.du222 = sp.diff(self.du22, self.x2)
        print('du222 : \n', self.du222)

class HyperbolicTangentFluidFractionSolution():
    def __init__(self):
        self.time  = sp.Symbol('time')
        self.x1 = sp.Symbol('x1')
        self.x2 = sp.Symbol('x2')
        self.height = sp.Symbol('height')
        self.b = sp.Symbol('b')
        self.L  = sp.Symbol('L')
        self.c = sp.Symbol('c')
        self.mean_alpha = sp.Symbol('mean_alpha')
        self.x10 = sp.Symbol('x10')
        self.x20 = sp.Symbol('x20')
        self.u_char = sp.Symbol('u_char')

    def VelocityField(self):

        self.alpha = self.mean_alpha * (1 + self.height * sp.tanh(self.b * (self.x1 - 0.5))) - self.c
        self.f_x1 = self.u_char * self.x1**2 * (1 - self.x1)**2
        self.f_x2 = self.u_char * self.x2**2 * (1 - self.x2)**2
        self.g_t = sp.cos(sp.pi * 1.0) * sp.exp(-1)
        self.df_x1 = sp.diff(self.f_x1, self.x1)
        self.df_x2 = sp.diff(self.f_x2, self.x2)
        self.u1_codina = self.f_x1 * self.df_x2 * self.g_t
        self.u2_codina = -self.df_x1 * self.f_x2 * self.g_t

        self.u1 = self.u1_codina / self.alpha
        self.u2 = self.u2_codina / self.alpha

    def DerivativeCalculator(self):
        self.VelocityField()
        print('u1 : \n', self.u1)
        print('u2 : \n', self.u2)
        self.du1dt = sp.diff(self.u1, self.time)
        print('du1dt : \n', self.du1dt)
        self.du2dt = sp.diff(self.u2, self.time)
        print('du2dt : \n', self.du2dt)
        print('alpha : \n', self.alpha)
        self.dalphat = sp.diff(self.alpha, self.time)
        print('dalphat : \n', self.dalphat)
        self.alpha1 = sp.diff(self.alpha, self.x1)
        print('alpha1 : \n', self.alpha1)
        self.alpha2 = sp.diff(self.alpha, self.x2)
        print('alpha2 : \n', self.alpha2)
        self.du11 = sp.diff(self.u1, self.x1)
        print('du11 : \n', self.du11)
        self.du12 = sp.diff(self.u1, self.x2)
        print('du12 : \n', self.du12)
        self.du111 = sp.diff(self.du11, self.x1)
        print('du111 : \n', self.du111)
        self.du112 = sp.diff(self.du11, self.x2)
        print('du112 : \n', self.du112)
        self.du121 = sp.diff(self.du12, self.x1)
        print('du121 : \n', self.du121)
        self.du122 = sp.diff(self.du12, self.x2)
        print('du122 : \n', self.du122)
        self.du21 = sp.diff(self.u2, self.x1)
        print('du21 : \n', self.du21)
        self.du22 = sp.diff(self.u2, self.x2)
        print('du22 : \n', self.du22)
        self.du211 = sp.diff(self.du21, self.x1)
        print('du211 : \n', self.du211)
        self.du212 = sp.diff(self.du21, self.x2)
        print('du212 : \n', self.du212)
        self.du221 = sp.diff(self.du22, self.x1)
        print('du221 : \n', self.du221)
        self.du222 = sp.diff(self.du22, self.x2)
        print('du222 : \n', self.du222)

HyperbolicTangentFluidFractionSolution().DerivativeCalculator()