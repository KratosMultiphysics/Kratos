import sympy
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .src.defines import CompressibleNavierStokesDefines as defs
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes. \
    src.quantity_converter import QuantityConverter


def ComputeDiffusiveFlux(primitives, params):
    """Calculate the diffusive flux matrix without shock capturing contribution."""

    # Auxiliary variables
    dim = params.dim
    vel = primitives.V

    # Calculate the viscous stress tensor
    mu = params.mu  # Dynamic viscosity
    beta = 0.0      # Null bulk viscosity (Stoke's assumption)
    tau_stress = CalculateViscousStressTensor(mu, beta, primitives, dim)

    # Calculate the heat flux vector
    lamb = params.lamb  # Thermal conductivity
    heat_flux = CalculateHeatFluxVector(lamb, primitives)

    # Define and fill the diffusive flux matrix
    G = sympy.zeros(dim+2, dim)
    G[1:-1, :] = -tau_stress.T
    G[-1, :] = heat_flux.T - vel.T*tau_stress

    G.simplify()
    return G


def ComputeDiffusiveFluxWithShockCapturing(primitives, params, sc_params):
    """
    Calculate the diffusive flux matrix with a physics-based shock
    capturing contribution. See:

    2018. Fernandez, Nguyen and Peraire
    A physics-based shock capturing methods for large-eddy simulation..

    """

    # Auxiliary variables
    dim = params.dim
    vel = primitives.V

    # Calculate the density flux
    grad_rho = QuantityConverter.density_gradient(primitives, params)
    mass_flux = CalculuateMassFluxVector(sc_params.alpha, grad_rho)

    # Calculate the viscous stress tensor
    mu = params.mu          # Dynamic viscosity
    mu += sc_params.mu      # Artificial dynamic viscosity
    beta = 0.0              # Null physical bulk viscosity (Stoke's assumption)
    beta += sc_params.beta  # Artificial bulk viscosity
    tau_stress = CalculateViscousStressTensor(mu, beta, primitives, dim)

    # Calculate the heat flux vector
    lamb = params.lamb      # Thermal conductivity
    lamb += sc_params.lamb  # Artificial thermal conductivity for shock capturing
    heat_flux = CalculateHeatFluxVector(lamb, primitives)

    # Define and fill the isotropic shock capturing diffusive flux matrix
    G = sympy.zeros(dim+2, dim)
    G[0,:] = mass_flux.T
    G[1:-1, :] = -tau_stress.T
    G[-1, :] = heat_flux.T - vel.T*tau_stress

    G.simplify()
    return G


def CalculuateMassFluxVector(alpha, grad_rho):
    """
    Auxiliary function to calculate mass flux vector f_rho = alpha * gradient(u)

    Mass diffusivity (alpha) is 0 in the Navier-Stokes equations but some shock
    capturing methods introduce a positive, non-zero value
    """
    return alpha * grad_rho


def CalculateViscousStressTensor(mu, beta, primitives, dim):
    """
    Auxiliary function to calculate the viscous stress tensor for the given
    dynamic and bulk viscosity values
    """


    # Calculate the viscous stress tensor
    # Note that the second viscosity coefficient is computed
    # as the bulk viscosity minus 2/3 of the dynamic one
    kappa = (beta - 2/3*mu)

    dynamic = mu * (primitives.grad_V + primitives.grad_V.transpose())
    bulk    = - kappa * sympy.trace(primitives.grad_V) * sympy.eye(dim)

    return dynamic + bulk


def CalculateHeatFluxVector(lamb, primitives):
    """Auxiliary function to calculate the heat flux vector with Fourier's law"""

    # Calculate the heat flux vector (Fourier's law q = -lambda * grad(theta))¡
    return - lamb * primitives.grad_T


def WriteInVoigtNotation(dim, tensor):
    """Auxiliary function to represent a 2nd order tensor in Voigt notation"""

    voigt_size = 3 * dim - 3
    voigt_tensor = defs.ZeroMatrix(voigt_size, 1)
    voigt_tensor[0] = tensor[0,0]
    voigt_tensor[1] = tensor[1,1]
    if dim == 2:
        voigt_tensor[2] = tensor[0,1]
    else:
        voigt_tensor[2] = tensor[2,2]
        voigt_tensor[3] = tensor[1,2]
        voigt_tensor[4] = tensor[0,2]
        voigt_tensor[5] = tensor[0,1]

    return voigt_tensor


def RevertVoigtNotation(dim, voigt_tensor):
    """Auxiliary function to set a 2nd order tensor from its Voigt representation"""

    # Revert the Voigt notation
    tensor = defs.ZeroMatrix(dim, dim)
    for d in range(dim):
        tensor[d, d] = voigt_tensor[d]
    if dim == 2:
        tensor[0, 1] = voigt_tensor[2]
        tensor[1, 0] = voigt_tensor[2]
    else:
        tensor[1, 2] = voigt_tensor[3]
        tensor[2, 1] = voigt_tensor[3]
        tensor[0, 2] = voigt_tensor[4]
        tensor[2, 0] = voigt_tensor[4]
        tensor[0, 1] = voigt_tensor[5]
        tensor[1, 0] = voigt_tensor[5]

    return tensor


def PrintDiffusiveFluxMatrix(G, params):
    """Auxiliary function to print the diffusive flux matrix (G)"""

    dim = params.dim
    print("The diffusive matrix is:\n")
    for ll in range(dim+2):
        for mm in range(dim):
            print("G[", ll, ",", mm, "]=", G[ll, mm], "\n")

    return 0
