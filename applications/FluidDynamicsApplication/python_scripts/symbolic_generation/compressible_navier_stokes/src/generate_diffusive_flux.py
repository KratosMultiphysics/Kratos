from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .src.defines import CompressibleNavierStokesDefines as defs


def ComputeDiffusiveFlux(dofs, dUdx, params):
    """Calculate the diffusive flux matrix without shock capturing contribution."""

    # Auxiliary variables
    dim = params.dim
    rho = dofs[0]
    mom = []
    vel = []
    for i in range(dim):
        mom.append(dofs[i + 1])
        vel.append(dofs[i + 1] / rho)
    e_tot = dofs[dim + 1]

    # Calculate the viscous stress tensor
    mu = params.mu  # Dynamic viscosity
    beta = 0.0      # Null bulk viscosity (Stoke's assumption)
    tau_stress = CalculateViscousStressTensor(mu, beta, rho, mom, dim, dUdx)

    # Calculate the heat flux vector
    c_v = params.c_v    # Specific heat at constant volume
    lamb = params.lamb  # Thermal conductivity
    heat_flux = CalculateHeatFluxVector(c_v, lamb, rho, mom, e_tot, dim, dUdx)

    # Define and fill the diffusive flux matrix
    G = defs.Matrix('G', dim + 2, dim, real=True)
    for j in range(dim):
        G[0, j] = 0.0
        G[dim + 1, j] = heat_flux[j]
        for i in range(dim):
            G[i + 1, j] = -tau_stress[j, i]
            G[dim + 1, j] -= vel[i] * tau_stress[i, j]

    return G


def ComputeDiffusiveFluxWithShockCapturing(dofs, dUdx, params, sc_params):
    """
    Calculate the diffusive flux matrix with a physics-based shock
    capturing contribution. See:

    2018. Fernandez, Nguyen and Peraire
    A physics-based shock capturing methods for large-eddy simulation..

    """

    # Auxiliary variables
    dim = params.dim
    rho = dofs[0]
    mom = []
    vel = []
    for i in range(dim):
        mom.append(dofs[i + 1])
        vel.append(dofs[i + 1] / rho)
    e_tot = dofs[dim + 1]

    # Calculate the density flux
    mass_flux = CalculuateMassFluxVector(sc_params.alpha, dUdx)

    # Calculate the viscous stress tensor
    mu = params.mu          # Dynamic viscosity
    mu += sc_params.mu      # Artificial dynamic viscosity
    beta = 0.0              # Null physical bulk viscosity (Stoke's assumption)
    beta += sc_params.beta  # Artificial bulk viscosity
    tau_stress = CalculateViscousStressTensor(mu, beta, rho, mom, dim, dUdx)

    # Calculate the heat flux vector
    c_v = params.c_v        # Specific heat at constant volume
    lamb = params.lamb      # Thermal conductivity
    lamb += sc_params.lamb  # Artificial thermal conductivity for shock capturing
    heat_flux = CalculateHeatFluxVector(c_v, lamb, rho, mom, e_tot, dim, dUdx)

    # Define and fill the isotropic shock capturing diffusive flux matrix
    G = defs.Matrix('G', dim + 2, dim, real=True)
    G[0, :] = mass_flux
    for j in range(dim):
        G[dim + 1, j] = heat_flux[j]
        for i in range(dim):
            G[i + 1, j] = -tau_stress[j,i]
            G[dim + 1, j] -= vel[i] * tau_stress[i,j]

    return G


def CalculuateMassFluxVector(alpha, dUdx):
    """
    Auxiliary function to calculate mass flux vector f_rho = alpha * gradient(u)

    Mass diffusivity (alpha) is 0 in the Navier-Stokes equations but some shock
    capturing methods introduce a positive, non-zero value
    """
    return alpha * dUdx[0, :]


def CalculateViscousStressTensor(mu, beta, rho, mom, dim, dUdx):
    """
    Auxiliary function to calculate the viscous stress tensor for the given
    dynamic and bulk viscosity values
    """

    # Calculate velocity divergence
    # Note that this is computed as div(mom/rho) = (dx(mom)*rho - mom*dx(rho))/rho**2
    div_vel = 0.0
    for d in range(dim):
        div_vel += (dUdx[d + 1, d] * rho - mom[d] * dUdx[0, d])
    div_vel /= rho**2

    # Calculate the viscous stress tensor
    # Note that the derivatives in here involve grad(mom/rho) = (dx(mom)*rho - mom*dx(rho))/rho**2
    tau_stress = defs.Matrix('tau_stress', dim, dim, real=True)
    for d1 in range(dim):
        for d2 in range(dim):
            dv1_dx2 = (dUdx[d1 + 1, d2] * rho - mom[d1] * dUdx[0, d2]) / rho**2
            dv2_dx1 = (dUdx[d2 + 1, d1] * rho - mom[d2] * dUdx[0, d1]) / rho**2
            tau_stress[d1, d2] = mu * (dv1_dx2 + dv2_dx1)
            if d1 == d2:
                # Note that in here the second viscosity coefficient is computed
                # as the bulk viscosity minus 2/3 of the dynamic one
                tau_stress[d1, d2] += (beta - 2.0 * mu / 3.0) * div_vel

    return tau_stress


def CalculateHeatFluxVector(c_v, lamb, rho, mom, e_tot, dim, dUdx):
    """Auxiliary function to calculate the heat flux vector with Fourier's law"""

    # Calculate the heat flux vector (Fourier's law q = -lambda * grad(theta))
    # Note that the temperature is expressed in terms of the total energy
    heat_flux = []
    for d in range(dim):
        aux_1 = (dUdx[dim + 1, d]*rho - e_tot * dUdx[0, d]) / rho**2
        aux_2 = 0.0
        for i in range(dim):
            aux_2 += mom[i] * dUdx[i + 1, d] / rho**2
            aux_2 -= mom[i]**2 * dUdx[0, d] / rho**3
        heat_flux.append(- (lamb / c_v) * (aux_1 - aux_2))

    return heat_flux


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
