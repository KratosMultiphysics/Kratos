from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

from sympy import *
import pprint

## Computation of the flux
def ComputeDiffusiveFlux(dofs, dUdx, params):
    print("\nCompute Diffusive Flux\n")

    ## Auxiliary variables
    dim = params["dim"]
    gamma = params["gamma"]
    rho = dofs[0]
    mom = []
    vel = []
    for i in range(dim):
        mom.append(dofs[i + 1])
        vel.append(dofs[i + 1] / rho)
    e_tot = dofs[dim + 1]

    ## Calculate the viscous stress tensor
    mu = params["mu"] # Dynamic viscosity
    tau_stress = CalculateViscousStressTensor(mu, rho, mom, dim, dUdx)

    ## Calculate the heat flux vector
    c_v = params["c_v"]	# Specific heat at constant volume
    lamb = params["lambda"] # Thermal conductivity
    heat_flux = CalculateHeatFluxVector(c_v, lamb, rho, mom, e_tot, dim, dUdx)

    ## Define and fill the diffusive flux matrix
    G = DefineMatrix('G', dim + 2, dim)
    for j in range(dim):
        G[0,j] = 0.0
        G[dim + 1, j] = heat_flux[j]
        for i in range(dim):
            G[i + 1, j] = -tau_stress[j,i]
            G[dim + 1, j] -= vel[i] * tau_stress[i,j]

    return G

def ComputeDiffusiveFluxIsotropicShockCapturing(dofs, dUdx, params, nu_sc, k_sc):
    print("\nCompute Diffusive Flux")
    print("\tAdd diffusive contribution contribution")

    ## Auxiliary variables
    dim = params["dim"]
    rho = dofs[0]
    mom = []
    vel = []
    for i in range(dim):
        mom.append(dofs[i + 1])
        vel.append(dofs[i + 1] / rho)
    e_tot = dofs[dim + 1]

    ## Calculate the viscous stress tensor
    mu = params["mu"] # Dynamic viscosity
    tau_stress = CalculateViscousStressTensor(mu, rho, mom, dim, dUdx)

    ## Calculate the shock capturing viscous stress tensor
    mu_sc = rho * nu_sc
    tau_stress_mu_sc = CalculateViscousStressTensor(mu_sc, rho, mom, dim, dUdx)

    ## Calculate the heat flux vector
    c_v = params["c_v"]	# Specific heat at constant volume
    lamb = params["lambda"] # Thermal conductivity
    heat_flux = CalculateHeatFluxVector(c_v, lamb, rho, mom, e_tot, dim, dUdx)

    ## Calculate the shock capturing flux vector
    alpha_sc = k_sc # TODO: Remove when using alpha_sc --> This is way misleading as alpha is thermal diffusivity and k (or lamda) is normally the conductivity
    lamb_sc = alpha_sc * rho * (gamma * c_v) # TODO: Check that this is c_p (Camilo has c_v here but this is not the definition of the conductivity)
    heat_flux_lamb_sc = CalculateHeatFluxVector(c_v, lamb_sc, rho, mom, e_tot, dim, dUdx)

    ## Define and fill the isotropic shock capturing diffusive flux matrix
    print("\tAdd isotropic shock capturing")
    G = DefineMatrix('G', dim + 2, dim)
    for j in range(dim):
        G[0,j] = 0.0
        G[dim + 1, j] = heat_flux[j] + heat_flux_lamb_sc
        # G[dim + 1, j] = (1 + rho * c_v * k_sc / lamb) * heat_flux[j]
        for i in range(dim):
            G[i + 1, j] = -(tau_stress[j,i] + tau_stress_mu_sc[j,i])
            # G[dim + 1, j] -= vel[i] * tau_stress[i,j] # Note that we do not add the viscous shock capturing to the energy equation
            G[dim + 1, j] -= vel[i] * (tau_stress[i,j] + tau_stress_mu_sc[i,j])
            # Old implementation
            # G[i + 1, j] = - (1 + rho * nu_sc / mu) * tau_stress[j,i]
            # G[dim + 1, j] -= vel[i] * (1 + rho * nu_sc / mu) *  tau_stress[i,j]

    return G

def ComputeDiffusiveFluxAnisotropicShockCapturing(dofs, dUdx, params, nu_sc, k_sc, nu_st, k_st, lin_m, lin_m_norm):
    print("\nCompute Diffusive Flux")
    print("\tAdd diffusive contribution contribution")

    ## Auxiliary variables
    # Note that in here the momentum is not linearized.
    # This is OK for the explicit case but must be modified to avoid differenctiating these terms in the implicit case
    dim = params["dim"]
    rho = dofs[0]
    mom = []
    vel = []
    mom_prod = 0.0
    for i in range(dim):
        mom.append(dofs[i + 1])
        vel.append(dofs[i + 1] / rho)
        mom_prod += dofs[i + 1]**2
    e_tot = dofs[dim + 1]

    ## Calculate the viscous stress tensor
    mu = params["mu"] # Dynamic viscosity
    tau_stress = CalculateViscousStressTensor(mu, rho, mom, dim, dUdx)

    ## Calculate the heat flux vector
    c_v = params["c_v"]	# Specific heat at constant volume
    lamb = params["lambda"] # Thermal conductivity
    gamma = params["gamma"] # Heat capacity ratio
    heat_flux = CalculateHeatFluxVector(c_v, lamb, rho, mom, e_tot, dim, dUdx)

    ## Define and fill the anisotropic shock capturing diffusive flux matrices
    print("\tAdd anisotropic shock capturing contribution")
    # Heat flux shock capturing terms
    I = eye(dim)
    q_sc_mat = zeros(dim,dim)
    S = zeros(dim, dim)
    for d1 in range(dim):
        for d2 in range(dim):
            S[d1,d2] = lin_m[d1] * lin_m[d2] / lin_m_norm
    O = I - S
    # q_sc_mat = I + (rho * c_v * k_sc / lamb) * O + (rho * c_v * k_st / lamb) * S
    # heat_flux_sc = []
    # for i in range(dim):
    #     aux = 0.0
    #     for j in range(dim):
    #         aux += q_sc_mat[i,j] * heat_flux[j]
    #     heat_flux_sc.append(aux)


    alpha_sc = k_sc # TODO: Remove when using alpha_sc --> This is way misleading as alpha is thermal diffusivity and k (or lamda) is normally the conductivity --> This is way misleading as alpha is thermal diffusivity and k (or lamda) is normally the conductivity
    lamb_sc = alpha_sc * rho * (gamma * c_v) # TODO: Check that this is c_p (Camilo has c_v here but this is not the definition of the conductivity)
    heat_flux_lamb_sc = CalculateHeatFluxVector(c_v, lamb_sc, rho, mom, e_tot, dim, dUdx)

    alpha_st = k_st # TODO: Remove when using alpha_sc --> This is way misleading as alpha is thermal diffusivity and k (or lamda) is normally the conductivity
    lamb_st = alpha_st * rho * (gamma * c_v) # TODO: Check that this is c_p (Camilo has c_v here but this is not the definition of the conductivity)
    heat_flux_lamb_st = CalculateHeatFluxVector(c_v, lamb_st, rho, mom, e_tot, dim, dUdx)

    heat_flux_sc = []
    for i in range(dim):
        aux = 0.0
        for j in range(dim):
            aux += I[i,j] * heat_flux[j] + O[i,j] * heat_flux_lamb_sc[j] + S[i,j] * heat_flux_lamb_st[j]
        heat_flux_sc.append(aux)

    # Viscous shock capturing terms
    # Calculate shock caturing viscous contribution
    mu_sc = rho * nu_sc
    tau_stress_mu_sc = CalculateViscousStressTensor(mu_sc, rho, mom, dim, dUdx)
    # Calculate stabilization viscous contribution
    mu_st = rho * nu_st
    tau_stress_mu_st = CalculateViscousStressTensor(mu_st, rho, mom, dim, dUdx)

    # Write viscous contributions in Voigt notation
    tau_stress_voigt = WriteInVoigtNotation(dim, tau_stress)
    tau_stress_mu_sc_voigt = WriteInVoigtNotation(dim, tau_stress_mu_sc)
    tau_stress_mu_st_voigt = WriteInVoigtNotation(dim, tau_stress_mu_st)

    # # Write the viscous stress in Voigt notation
    # voigt_size = 3 * dim - 3
    # tau_stress_voigt = zeros(voigt_size, 1)
    # tau_stress_voigt[0] = tau_stress[0,0]
    # tau_stress_voigt[1] = tau_stress[1,1]
    # if dim == 2:
    #     tau_stress_voigt[2] = tau_stress[0,1]
    # else:
    #     tau_stress_voigt[2] = tau_stress[2,2]
    #     tau_stress_voigt[3] = tau_stress[1,2]
    #     tau_stress_voigt[4] = tau_stress[0,2]
    #     tau_stress_voigt[5] = tau_stress[0,1]

    # Compute the anisotropic shock capturing tensors
    voigt_size = 3 * dim - 3
    I_4_voigt = eye(voigt_size)
    S_4_voigt = zeros(voigt_size, voigt_size)
    for d1 in range(dim):
        for d2 in range(dim):
            S_4_voigt[d1,d2] = lin_m[d1] * lin_m[d2] / lin_m_norm
    if dim == 2:
        S_4_voigt[2,2] = lin_m[0] * lin_m[1] / lin_m_norm
    else:
        S_4_voigt[3,3] = lin_m[1] * lin_m[2] / lin_m_norm
        S_4_voigt[4,4] = lin_m[0] * lin_m[2] / lin_m_norm
        S_4_voigt[5,5] = lin_m[0] * lin_m[1] / lin_m_norm
    O_4_voigt = I_4_voigt - S_4_voigt

    # tau_sc_mat = I_4_voigt + (rho * nu_sc / mu) * O_4_voigt + (rho * nu_st / mu) * S_4_voigt
    # tau_stress_voigt_sc = tau_sc_mat * tau_stress_voigt

    # # Revert the Voigt notation
    # tau_stress_sc = zeros(dim, dim)
    # for d in range(dim):
    #     tau_stress_sc[d,d] = tau_stress_voigt_sc[d]
    # if dim == 2:
    #     tau_stress_sc[0,1] = tau_stress_voigt_sc[2]
    #     tau_stress_sc[1,0] = tau_stress_voigt_sc[2]
    # else:
    #     tau_stress_sc[1,2] = tau_stress_voigt_sc[3]
    #     tau_stress_sc[2,1] = tau_stress_voigt_sc[3]
    #     tau_stress_sc[0,2] = tau_stress_voigt_sc[4]
    #     tau_stress_sc[2,0] = tau_stress_voigt_sc[4]
    #     tau_stress_sc[0,1] = tau_stress_voigt_sc[5]
    #     tau_stress_sc[1,0] = tau_stress_voigt_sc[5]

    # Calculate the modified viscous flux
    tau_stress_voigt_sc = I_4_voigt * tau_stress_voigt + O_4_voigt * tau_stress_mu_sc_voigt + S_4_voigt * tau_stress_mu_st_voigt
    tau_stress_sc = RevertVoigtNotation(dim, tau_stress_voigt_sc)

    ## Output the shock capturing flux matrix
    G = DefineMatrix('G', dim + 2, dim)
    for j in range(dim):
        G[0,j] = 0.0
        G[dim + 1, j] = heat_flux_sc[j]
        for i in range(dim):
            G[i + 1, j] = - tau_stress_sc[j,i]
            G[dim + 1, j] -= vel[i] * tau_stress[i,j] # Note that we do not add the viscous shock capturing to the energy equation
            # G[dim + 1, j] -= vel[i] * tau_stress_sc[i,j]

    return G

## Auxiliary function to calculate the viscous stress tensor for a given viscosity value
def CalculateViscousStressTensor(mu, rho, mom, dim, dUdx):
    ## Calculate velocity divergence
    ## Note that this is computed as div(mom/rho) = (dx(mom)*rho - mom*dx(rho))/rho**2
    div_vel = 0.0
    for d in range(dim):
        div_vel += (dUdx[d + 1, d] * rho - mom[d] * dUdx[0, d])
    div_vel /= rho**2

    ## Calculate the viscous stress tensor
    ## Note that the derivatives in here involve grad(mom/rho) = (dx(mom)*rho - mom*dx(rho))/rho**2
    tau_stress = DefineMatrix('tau_stress', dim, dim)
    for d1 in range(dim):
        for d2 in range(dim):
            dv1_dx2 = (dUdx[d1 + 1, d2] * rho - mom[d1] * dUdx[0,d2]) / rho**2
            dv2_dx1 = (dUdx[d2 + 1, d1] * rho - mom[d2] * dUdx[0,d1]) / rho**2
            tau_stress[d1, d2] = mu * (dv1_dx2 + dv2_dx1)
            if d1 == d2:
                tau_stress[d1, d2] -= 2.0 * mu * div_vel / 3.0

    return tau_stress

def CalculateHeatFluxVector(c_v, lamb, rho, mom, e_tot, dim, dUdx):
    ## Calculate the heat flux vector (Fourier's law q = -lambda * grad(theta))
    ## Note that the temperature is expressed in terms of the total energy
    heat_flux = []
    for d in range(dim):
        aux_1 = (dUdx[dim + 1, d]*rho - e_tot * dUdx[0,d]) / rho**2
        aux_2 = 0.0
        for i in range(dim):
            aux_2 += mom[i] * dUdx[i + 1, d] / rho**2
            aux_2 -= mom[i]**2 * dUdx[0, d] / rho**3
        heat_flux.append(- (lamb / c_v) * (aux_1 - aux_2))

    return heat_flux

def WriteInVoigtNotation(dim, tensor):
    voigt_size = 3 * dim - 3
    voigt_tensor = zeros(voigt_size, 1)
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
    # Revert the Voigt notation
    tensor = zeros(dim, dim)
    for d in range(dim):
        tensor[d,d] = voigt_tensor[d]
    if dim == 2:
        tensor[0,1] = voigt_tensor[2]
        tensor[1,0] = voigt_tensor[2]
    else:
        tensor[1,2] = voigt_tensor[3]
        tensor[2,1] = voigt_tensor[3]
        tensor[0,2] = voigt_tensor[4]
        tensor[2,0] = voigt_tensor[4]
        tensor[0,1] = voigt_tensor[5]
        tensor[1,0] = voigt_tensor[5]
    
    return tensor

## Printing the Diffusive Matrix G
def printG(G,params):
    dim = params["dim"]
    print("The diffusive matrix is:\n")
    for ll in range(dim+2):
        for mm in range(dim):
            print("G[",ll,",",mm,"]=",G[ll,mm],"\n")

    return 0
