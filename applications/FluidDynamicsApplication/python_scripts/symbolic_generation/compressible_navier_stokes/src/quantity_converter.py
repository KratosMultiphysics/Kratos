import sympy

import KratosMultiphysics.sympy_fe_utilities as KratosSympy

class QuantityConverter:
    """
    Converts between conservative and primitive variables.

    When represented as vectors, their meanings are:
    - `V`: Primitive magnitudes: `[p, vel, T]'`
    - `U`: Conservative magnitudes: `[rho, mom, e_tot]'`

    When represented as classes, they are named `primitives` and `conservatives`
    """
    def __init__(self):
        raise RuntimeError("This class is purely static.")

    @classmethod
    def Primitives(cls, U, params):
        vel = cls.velocity(U)
        T = cls.temperature(U, params, vel=vel)
        p = cls.pressure(U, params, T=T)

        vel.simplify()
        T.simplify()
        p.simplify()

        return (p, vel, T)

    @classmethod
    def PrimitivesGradients(cls, U, grad_U, params):
        T = cls.temperature(U, params)

        grad_V = cls.velocity_gradient(U, grad_U)
        grad_T = cls.temperature_gradient(U, grad_U, params, T=T)
        grad_P = cls.pressure_gradient(U, grad_U, params, T=T, grad_T=grad_T)

        grad_V.simplify()
        grad_T.simplify()
        grad_P.simplify()

        return (grad_P.T, grad_V.T, grad_T.T)

    @classmethod
    def density(cls, primitives, params):
        "Returns the density as a function of primitive variables"
        return primitives.P / (cls.gas_constant_R(params) * primitives.T)

    @classmethod
    def momentum(cls, primitives, params=None, density=None):
        "Returns the momentum as a function of primitive variables"
        if density is None and params is None:
            raise RuntimeError("Either params or density must be specified")

        if density is not None:
            density = cls.density(primitives, params)

        return primitives / density

    @classmethod
    def total_energy(cls, primitives, params, rho=None):
        """Returns the density as a function of primitive variables"""
        if rho is not None:
            rho = cls.density(primitives, params)

        v_2 = sum([vi*vi for vi in primitives.V])
        T = primitives.T

        return rho * (0.5 * v_2 + params.c_v*T)

    @classmethod
    def velocity(cls, U):
        return sympy.Matrix(U[1:-1]) / U[0]

    @classmethod
    def temperature(cls, U, params, vel=None):
        rho = U[0]
        if vel is None:
            vel = cls.velocity(U)
        e_tot = U[-1]

        e_kinetic = 0.5 * rho * sum([v**2 for v in vel])
        return (e_tot - e_kinetic) / (rho * params.c_v)

    @classmethod
    def pressure(cls, U, params, T=None):
        rho = U[0]
        R = cls.gas_constant_R(params)
        if T is None:
            T = cls.temperature(U, params)
        return rho * R * T

    @classmethod
    def _ValidateUAndDUShapes(cls, U, DU):
        """Asserts to prevent cryptic sympy error messages for a common error"""
        class LazyErrorMsg:
            def __init__(self, expected, obtained):
                self.rows_expected = expected[0]
                self.cols_expected = expected[1]
                self.rows = obtained[0]
                self.cols = obtained[1]

            def __str__(self):
                return "[{}, {}] does not match expected [{},{}]".format(
                    self.rows, self.cols, self.rows_expected, self.cols_expected
                )

        blocksize = U.shape[0]
        dim = blocksize - 2
        msg = LazyErrorMsg((blocksize, dim), DU.shape)
        assert DU.shape[0] == blocksize, msg
        assert DU.shape[1] == dim, msg

    @classmethod
    def velocity_gradient(cls, U, DU):
        """
        Velocity gradient. Gradients defined as:

        grad_f := df_j/dx_i

        """
        cls._ValidateUAndDUShapes(U, DU)
        dim = U.shape[0] - 2

        rho = U[0]
        mom = sympy.Matrix(dim, 1, lambda i,_: U[1+i])

        grad_rho = sympy.Matrix(1,   dim, lambda _,j: DU[0, j])
        grad_mom = sympy.Matrix(dim, dim, lambda i,j: DU[1+i, j])

        return (grad_mom*rho - mom*grad_rho) / rho**2


    @classmethod
    def temperature_gradient(cls, U, DU, params, vel=None, T=None, grad_vel=None):
        """
        Temperature gradient.  Gradients defined as:

        grad_f := df_j/dx_i

        Explanation
        ---
        Temperature is defined as:
        ```
            T = (e_t - e_k) / (rho * c_v)
            [with e_k = 1/2 * rho * V²]
        ```
        where `e_k` is the kinetic energy.

        The gradient is:
        ```
            grad(T) = grad(e_t - e_k) / (rho * c_v)   -   (e_t - e_k)/(rho²*c_v) * grad(rho)
        ```
        Rearranging yields:
        ```
            grad(T) = (grad(e_t) - grad(e_k)) / (rho * c_v)   -   T/rho * grad(rho)
            [with grad(e_k) = rho*transp(V)*grad(V) + 1/2 * V² * grad(rho)]
        ```

        """
        cls._ValidateUAndDUShapes(U, DU)

        if vel is None:
            vel = cls.velocity(U)
        V2 = vel.transpose() * vel

        if grad_vel is None:
            grad_vel = cls.velocity_gradient(U, DU)

        if T is None:
            T = cls.temperature(U, params, vel)

        rho = U[0]
        grad_rho = sympy.Matrix(DU[0, :])
        grad_e_tot = sympy.Matrix(DU[-1, :])

        grad_e_k = rho*vel.transpose()*grad_vel + 0.5 * V2 * grad_rho
        return (grad_e_tot - grad_e_k)/(rho*params.c_v) - T/rho * grad_rho

    @classmethod
    def pressure_gradient(cls, U, DU, params, T=None, grad_T=None):
        """
        Pressure gradient.  Gradients defined as:

        grad_f := df_j/dx_i

        """
        cls._ValidateUAndDUShapes(U, DU)

        if T is None:
            T = cls.temperature(U, params)

        if grad_T is None:
            grad_T = cls.temperature_gradient(U, DU, params, T=T)

        rho = U[0]
        grad_rho = sympy.Matrix(DU[0, :])

        R = cls.gas_constant_R(params)

        return (grad_rho*T + rho*grad_T) * R

    @classmethod
    def density_gradient(cls, primitives, params):
        R =  cls.gas_constant_R(params)
        return (primitives.T*primitives.grad_P - primitives.grad_T*primitives.P) / (R * primitives.T**2)

    @classmethod
    def gas_constant_R(cls, params):
        return (params.gamma - 1) * params.c_v

    @classmethod
    def dVdU(cls, U, primitives, params):
        """
        Returns the derivative of the primitive vector V with respect to the conservative vector U.

        returns `dV_i/dU_j`
        """
        blocksize = primitives.ndims+2
        V = primitives.AsVector()
        QuantityConverter.SubstitutePrimitivesWithConservatives(V, primitives, U, None, params, substitute_gradients=False)
        D = sympy.zeros(blocksize, blocksize)
        for i in range(blocksize):
            for j in range(blocksize):
                D[i,j] = sympy.diff(V[i], U[j])
        return D

    @classmethod
    def SubstitutePrimitivesWithConservatives(cls, expr, primitives, U, grad_U, params, substitute_gradients=True):
        """
        Substitutes the primitive variables in the expression `expr` with their definitions as functions of constervative variables `U`
        """
        (P, V, T) = QuantityConverter.Primitives(U, params)

        KratosSympy.SubstituteScalarValue(expr, primitives.P, P)
        KratosSympy.SubstituteMatrixValue(expr, primitives.V, V)
        KratosSympy.SubstituteScalarValue(expr, primitives.T, T)

        if not substitute_gradients:
            return

        (grad_P, grad_V, grad_T) = QuantityConverter.PrimitivesGradients(U, grad_U, params)

        KratosSympy.SubstituteMatrixValue(expr, primitives.grad_P, grad_P)
        KratosSympy.SubstituteMatrixValue(expr, primitives.grad_V, grad_V)
        KratosSympy.SubstituteMatrixValue(expr, primitives.grad_T, grad_T)