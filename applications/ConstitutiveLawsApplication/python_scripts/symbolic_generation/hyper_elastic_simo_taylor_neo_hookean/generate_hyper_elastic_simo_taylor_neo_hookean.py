import sympy
import KratosMultiphysics.sympy_fe_utilities as KratosSympy

# Symbolic generation settings
mode = "c"
# dim_vect = [2]
dim_vect = [2, 3]
aux_filename = "hyper_elastic_simo_taylor_neo_hookean_"

# Set the symbols for the required material parameters
Kappa = sympy.Symbol("Kappa") # Equivalent bulk modulus
Mu = sympy.Symbol("Mu") # 2nd Lame constant (shear modulus)

for dim in dim_vect:
    if dim == 2:
        strain_size = 3
    elif dim == 3:
        strain_size = 6
    else:
        raise ValueError("Wrong dimension {}.".format(dim))

    # Definition of a vector with the Voigt notation symbols to be used in the cpp
    E_voigt = KratosSympy.DefineVector("rStrain", strain_size)

    # Define the required strain measures
    E = KratosSympy.ConvertVoigtStrainToMatrix(E_voigt) # Green-Lagrange strain tensor
    C = 2*E + sympy.eye(dim, dim) # Right Cauchy-Green strain tensor
    J = sympy.sqrt(sympy.det(C).doit()) # Jacobian of the deformation
    devC = (J**(-sympy.Rational(2,dim)))*C # Deviatoric Right Cauchy-Green strain tensor

    # Set the definition of the Helmholtz free energy functional corresponding to the Simo-Taylor model of 1992
    # The energy is based in a volumetric/deviatoric decomposition of the Helmholtz energy functional such that a monotonic isocoric response is obtained
    # More information can be found in Simo et al. (DOI: 10.1016/0045-7825(85)90033-7) and Simo et al. (DOI: 10.1016/0045-7825(91)90100-K) or in Scovazzi et al. 2015 (DOI: 10.1002/nme.5138)
    psi = (Kappa/4)*(J**2-1.0) - (Kappa/2)*sympy.ln(J) # Volumetric contribution
    psi += (Mu/2)*(sympy.Trace(devC)-dim) # Deviatoric contribution

    # Bonet's book material to test
    # lamb = Kappa - 2.0*Mu/3.0
    # psi = 0.5*Mu*(sympy.Trace(C) - dim) - Mu*sympy.ln(J) + lamb*0.5*((sympy.ln(J))**2)

    # Calculate the second Piola-Kirchhoff stress vector differentiation
    S_voigt = sympy.Matrix(sympy.zeros(strain_size, 1))
    for i in range(strain_size):
        S_voigt[i, 0] = sympy.diff(psi, E_voigt[i, 0])

    # Calculate the constitutive matrix differentiation
    D_voigt = sympy.Matrix(sympy.zeros(strain_size,strain_size))
    for i in range(strain_size):
        for j in range(strain_size):
            D_voigt[i, j] = sympy.simplify(sympy.diff(S_voigt[i], E_voigt[j]))

    # Output the computed PK2 stress and constitutive matrx
    S_out = KratosSympy.OutputVector_CollectingFactors(S_voigt, "rStressVector", mode)
    D_out = KratosSympy.OutputMatrix_CollectingFactors(D_voigt, "rConstitutiveMatrix", mode)
    output_filename = aux_filename + ("3d.cpp" if dim == 3 else "plane_strain_2d.cpp")
    template_filename = aux_filename + ("3d_template.cpp" if dim == 3 else "plane_strain_2d_template.cpp")
    outstring = open(template_filename).read()

    # Replace the computed PK2 stress and constitutive matrix in the template outstring
    outstring = outstring.replace(f"//substitute_PK2_stress", S_out)
    outstring = outstring.replace(f"//substitute_PK2_constitutive_matrix", D_out)

    # Write the modified template
    with open(output_filename, 'w') as f:
        f.write(outstring)
