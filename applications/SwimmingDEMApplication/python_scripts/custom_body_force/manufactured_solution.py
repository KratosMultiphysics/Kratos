import KratosMultiphysics

def CreateManufacturedSolution(custom_settings):
    return ManufacturedSolution(custom_settings)

class ManufacturedSolution():
    def __init__(self, settings):
        '''
        This is a base class to build manufactured fluid solutions.
        At least, it should return the body force and the velocity.
        The input viscosity is the DYNAMIC viscosity
        NOTE: the operators are implemented for the 2D case. It could be extended to the 3D case.
        '''

        default_settings = KratosMultiphysics.Parameters("""
            {
                "viscosity"   : 1.0e-2,
                "density"     : 1.0
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.rho = settings["density"].GetDouble()
        self.nu = settings["viscosity"].GetDouble() / self.rho

    # Public methods

    def BodyForce(self, x1, x2, x3, t):
        return [self.body_force1(x1, x2, t), self.body_force2(x1, x2, t), self.body_force3(x1, x2, t)]

    def Velocity(self, x1, x2, x3, t):
        return [self.u1(x1, x2, t), self.u2(x1, x2, t), self.u3(x1, x2, t)]

    def Pressure(self, x1, x2, x3, t):
        return self.p(x1, x2, t)

    # Operators

    def body_force1(self, x1, x2, t):
        return self.du1dt(x1, x2, t) + self.convective1(x1, x2, t) + 1 / self.rho * self.press_grad1(x1, x2, t) - self.nu * self.laplacian1(x1, x2, t)

    def body_force2(self, x1, x2, t):
        return self.du2dt(x1, x2, t) + self.convective2(x1, x2, t) + 1 / self.rho * self.press_grad2(x1, x2, t) - self.nu * self.laplacian2(x1, x2, t)

    def body_force3(self, x1, x2, t):
        return 0.0

    def convective1(self, x1, x2, t):
        return self.u1(x1, x2, t) * self.du11(x1, x2, t) + self.u2(x1, x2, t) * self.du12(x1, x2, t)

    def convective2(self, x1, x2, t):
        return self.u1(x1, x2, t) * self.du21(x1, x2, t) + self.u2(x1, x2, t) * self.du22(x1, x2, t)

    def laplacian1(self, x1, x2, t):
        return self.du111(x1, x2, t) + self.du122(x1, x2, t)

    def laplacian2(self, x1, x2, t):
        return self.du211(x1, x2, t) + self.du222(x1, x2, t)

    def press_grad1(self, x1, x2, t):
        return self.dp1(x1, x2, t)

    def press_grad2(self, x1, x2, t):
        return self.dp2(x1, x2, t)

    # Velocity and derivatives

    def u1(self, x1, x2, t):
        """ Velocity
        """
        raise Exception("Method not implemented")

    def u2(self, x1, x2, t):
        """ Velocity
        """
        raise Exception("Method not implemented")

    def u3(self, x1, x2, t):
        return 0.0

    def du1dt(self, x1, x2, t):
        raise Exception("Method not implemented")

    def du2dt(self, x1, x2, t):
        raise Exception("Method not implemented")

    def du11(self, x1, x2, t):
        raise Exception("Method not implemented")

    def du12(self, x1, x2, t):
        raise Exception("Method not implemented")

    def du21(self, x1, x2, t):
        raise Exception("Method not implemented")

    def du22(self, x1, x2, t):
        raise Exception("Method not implemented")

    def du111(self, x1, x2, t):
        raise Exception("Method not implemented")

    def du122(self, x1, x2, t):
        raise Exception("Method not implemented")

    def du211(self, x1, x2, t):
        raise Exception("Method not implemented")

    def du222(self, x1, x2, t):
        raise Exception("Method not implemented")

    # Pressure and derivatives

    def p(self, x1, x2, t):
        '''
        By default, pressure is 0
        '''
        return 0.0

    def dp1(self, x1, x2, t):
        '''
        By default, pressure is 0
        '''
        return 0.0

    def dp2(self, x1, x2, t):
        '''
        By default, pressure is 0
        '''
        return 0.0
