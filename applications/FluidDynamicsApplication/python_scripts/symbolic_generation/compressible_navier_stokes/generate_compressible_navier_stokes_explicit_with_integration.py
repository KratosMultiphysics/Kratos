import pathlib
import os

import KratosMultiphysics
import sympy
import KratosMultiphysics.sympy_fe_utilities as KratosSympy

from params_dict import FormulationParameters, ShockCapturingParameters, ShockCapturingNodalParameters
from geometry_data import GeometryDataFactory
import generate_convective_flux
import generate_diffusive_flux
import generate_source_term
import generate_stabilization_matrix

class SymbolicGenerator:
    def __init__(self, settings):
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())
        self.geometry = GeometryDataFactory(self.settings["geometry"].GetString())
        self._GenerateFiles()

    @classmethod
    def GetDefaultParameters(cls):
        return KratosMultiphysics.Parameters("""
        {
            "mode": "c",
            "explicit": true,
            "do_simplifications": false,
            "geometry": "triangle",
            "shock_capturing": true,
            "subscales": {
                "ASGS" : true,
                "OSS" : true
            },
            "template_filename" : "PLEASE PROVIDE A template_filename",
            "output_filename"   : "symbolic_generator_name_not_provided.cpp",
            "echo_level" : 1
        }""")

    def _print(self, lvl, *text):
        if self.settings["echo_level"].GetInt() >= lvl:
            print(*text)

    def _GenerateFiles(self):
        with open(self.settings["template_filename"].GetString(), "r") as template:
            self.outstring = template.read()

        # Checking if outfile is valid
        self.outputfile = open(self.settings["output_filename"].GetString(), "w")

    def _ComputeNonLinearOperator(self, A, H, S, Ug):
        L = KratosSympy.Matrix(KratosSympy.zeros(self.geometry.blocksize, 1))
        for j in range(self.geometry.ndims):
            # Convective operator product (A x grad(U))
            A_j = A[j]
            H_j = H.col(j)
            L += A_j * H_j
            # Diffusive flux
            # Note that the diffusive flux is not added as it will involve 2nd order derivatives that vanish when introducing the linear FE discretization
        # Source term addition
        L -= S * Ug

        return L

    def _ComputeNonLinearAdjointOperator(self, A, H, Q, S, Ug, V):
        L_adj = KratosSympy.Matrix(KratosSympy.zeros(self.geometry.blocksize, 1))
        for j in range(self.geometry.ndims):
            Q_j = Q.col(j)
            H_j = H.col(j)
            # Convective operator product
            A_j_trans = A[j].transpose()
            L_adj += A_j_trans * Q_j
            aux_conv = KratosSympy.Matrix(KratosSympy.zeros(self.geometry.blocksize, self.geometry.blocksize))
            for m in range(self.geometry.blocksize):
                for n in range(self.geometry.blocksize):
                    A_j_trans_mn = A_j_trans[m,n]
                    for l in range(self.geometry.blocksize):
                        aux_conv[m,n] += sympy.diff(A_j_trans_mn, Ug[l]) * H_j[l]
            L_adj += aux_conv * V
            # Diffusive operator product
            # Note that the adjoint diffusive flux is not added as it will involve 2nd order derivatives that vanish when introducing the linear FE discretization
        # Source term addition
        L_adj += S.transpose() * V

        return L_adj

    def _ComputeVariationalFormulation(self, A, acc, G, H, L_adj, Q, S, Ug, V):
        # Mass (inertial) term - FE scale (only computed in the implicit case)
        if not self.is_explicit:
            n1 = - V.transpose()*acc

        # Convective term - FE scale
        conv_flux = KratosSympy.zeros(self.geometry.blocksize, 1)
        for j in range(self.geometry.ndims):
            conv_flux += A[j] * H.col(j)
        n2 = - V.transpose() * conv_flux

        # Diffusive term - FE scale
        n3 = KratosSympy.Matrix(KratosSympy.zeros(1,1))
        for j in range(self.geometry.ndims):
            for k in range(self.geometry.blocksize):
                n3[0,0] += Q[k,j] * G[k,j]

        # Source term - FE scale
        n4 = V.transpose() * (S * Ug)

        # VMS_adjoint - Subscales
        subscales = KratosSympy.DefineVector('subscales', self.geometry.blocksize)
        n5 = L_adj.transpose() * subscales

        # Variational formulation (Galerkin functional)
        self._print(1, "\nCompute variational formulation\n")
        if not self.is_explicit:
            rv = n1 + n2 + n3 + n4 + n5 # Implicit case (includes the inertial term n1)
        else:
            rv = n2 + n3 + n4 + n5         # Explicit case (without inertial term n1)

        return (rv, subscales)

    def _ComputeProjectionsAtGaussPoint(self, acc, bdf, dUdt, f, forcing_terms, H, i_gauss, mg, projections, res, rg, U, Ug, Un, Unn):
        ## Get Gauss point geometry data
        N = KratosSympy.DefineVector('N', self.geometry.nnodes)
        for i_node in range(self.geometry.nnodes):
            N[i_node] = self.geometry.N()[i_gauss, i_node]

        ## Data interpolation at the gauss point
        U_gauss = U.transpose() * N
        f_gauss = forcing_terms["force"].transpose() * N
        r_gauss = (forcing_terms["thermal"].transpose()*N)[0]
        mass_gauss = (forcing_terms["mass"].transpose()*N)[0]

        if self.is_explicit:
            # In the explicit case, the acceleration is linearised taking the previous step one
            # Note that in the explicit case this acceleration is only used in the calculation of the stabilization terms
            acc_gauss = dUdt.transpose()*N
        else:
            # In the implicit case, calculate the time derivatives with the BDF2 formula
            acc_gauss = (bdf[0] * U + bdf[1] * Un + bdf[2] * Unn).transpose()*N

        ## Gradients computation
        grad_U = KratosSympy.DfjDxi(self.geometry.DN(), U).transpose()

        ## Substitute the symbols in the residual
        res_gauss = res.copy()
        KratosSympy.SubstituteMatrixValue(res_gauss, Ug, U_gauss)
        KratosSympy.SubstituteMatrixValue(res_gauss, acc, acc_gauss)
        KratosSympy.SubstituteMatrixValue(res_gauss, H, grad_U)
        KratosSympy.SubstituteMatrixValue(res_gauss, f, f_gauss)
        KratosSympy.SubstituteScalarValue(res_gauss, rg, r_gauss)
        KratosSympy.SubstituteScalarValue(res_gauss, mg, mass_gauss)

        ## Add the projection contributions
        for i_node in range(self.geometry.nnodes):
            # Note that the weights will be added later on in the cpp
            projections["rho"][i_node] += N[i_node] * res_gauss[0]
            for d in range(self.geometry.ndims):
                projections["momentum"][i_node * self.geometry.ndims + d] += N[i_node] * res_gauss[1 + d]
            projections["energy"][i_node] += N[i_node] * res_gauss[self.geometry.ndims + 1]

    def _OutputProjections(self, outstring, res_rho_proj, res_mom_proj, res_tot_ener_proj):
        mode = self.settings["mode"].GetString()
        res_rho_proj_out = KratosSympy.OutputVector_CollectingFactors(res_rho_proj, "rho_proj", mode)
        res_mom_proj_out = KratosSympy.OutputVector_CollectingFactors(res_mom_proj, "mom_proj", mode)
        res_tot_ener_proj_out = KratosSympy.OutputVector_CollectingFactors(res_tot_ener_proj, "tot_ener_proj", mode)
        dim = self.geometry.ndims
        outstring = outstring.replace("//substitute_rho_proj_{}D".format(dim), res_rho_proj_out)
        outstring = outstring.replace("//substitute_mom_proj_{}D".format(dim), res_mom_proj_out)
        outstring = outstring.replace("//substitute_tot_ener_proj_{}D".format(dim), res_tot_ener_proj_out)
        return outstring

    def _SubstituteSubscales(self, res, res_proj, rv, subscales, subscales_type, Tau):
        rv_gauss = rv.copy()
        if subscales_type == "ASGS":
            subs = Tau * res
        elif subscales_type == "OSS":
            subs =  Tau * (res - res_proj)
        else:
            raise ValueError("Unrecognized subscales type: {}".format(subscales_type))

        KratosSympy.SubstituteMatrixValue(rv_gauss, subscales, subs)
        return rv_gauss

    def _ComputeResidualAtGaussPoint(self, acc, bdf, dUdt, f, forcing_terms, H, i_gauss, mg, params, Q, res_proj, ResProj, rg, rv_gauss, sc_nodes, sc_params, subscales_type, Tau, U, Ug, Un, Unn, V, w):
        self._print(1, "\tGauss point: " + str(i_gauss))

        ## Get Gauss point geometry data
        Ng = sympy.Matrix(self.geometry.nnodes, 1, lambda i,_: self.geometry.N()[i_gauss, i])

        ## Data interpolation at the gauss point
        U_gauss = U.transpose() * Ng
        w_gauss = w.transpose() * Ng
        f_gauss    = forcing_terms["force"].transpose() * Ng
        r_gauss    = (forcing_terms["thermal"].transpose()*Ng)[0]
        mass_gauss = (forcing_terms["mass"].transpose()*Ng)[0]
        alpha_sc_gauss = (sc_nodes.alpha.transpose()*Ng)[0]
        mu_sc_gauss    = (sc_nodes.mu.transpose()*Ng)[0]
        beta_sc_gauss  = (sc_nodes.beta.transpose()*Ng)[0]
        lamb_sc_gauss  = (sc_nodes.lambda_.transpose()*Ng)[0]
        if not self.is_explicit:
            # In the implicit case, calculate the time derivatives with the BDF2 formula
            acc_gauss = (bdf[0] * U + bdf[1] * Un + bdf[2] * Unn).transpose()*Ng
        else:
            # In the explicit case, the acceleration is linearised taking the previous step one
            # Note that in the explicit case this acceleration is only used in the calculation of the stabilization terms
            acc_gauss = dUdt.transpose()*Ng

        ## Gauss pt. stabilization matrix calculation
        if self.settings["shock_capturing"].GetBool():
            tau_gauss = generate_stabilization_matrix.ComputeStabilizationMatrixOnGaussPoint(params, U_gauss, f_gauss, r_gauss, mu_sc_gauss, lamb_sc_gauss)
        else:
            tau_gauss = generate_stabilization_matrix.ComputeStabilizationMatrixOnGaussPoint(params, U_gauss, f_gauss, r_gauss)

        ## If OSS, residual projections interpolation
        if subscales_type == "OSS":
            res_proj_gauss = ResProj.transpose() * Ng

        ## Gradients computation
        grad_U = KratosSympy.DfjDxi(self.geometry.DN(), U).transpose()
        grad_w = KratosSympy.DfjDxi(self.geometry.DN(), w).transpose()

        self._print(1, "\t- Substitution in the variational formulation")
        KratosSympy.SubstituteMatrixValue(rv_gauss, Ug, U_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, acc, acc_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, H, grad_U)
        KratosSympy.SubstituteMatrixValue(rv_gauss, V, w_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, Q, grad_w)
        KratosSympy.SubstituteMatrixValue(rv_gauss, Tau, tau_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, f, f_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, rg, r_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, mg, mass_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.alpha, alpha_sc_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.mu, mu_sc_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.beta, beta_sc_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.lambda_, lamb_sc_gauss)
        if subscales_type == "OSS":
            KratosSympy.SubstituteMatrixValue(rv_gauss, res_proj, res_proj_gauss)

        ## Accumulate in the total value
        return rv_gauss

    def _SubscalesTypes(self):
        return [name for (name, enabled) in self.settings["subscales"].items() if enabled]

    def _ComputeLHSandRHS(self, rv_tot, U, w):
        ## Set the DOFs and test function matrices to do the differentiation
        dofs = KratosSympy.Matrix(KratosSympy.zeros(self.geometry.ndofs, 1))
        testfunc = KratosSympy.Matrix(KratosSympy.zeros(self.geometry.ndofs, 1))
        for i in range(0, self.geometry.nnodes):
            for j in range(0,self.geometry.blocksize):
                dofs[i*self.geometry.blocksize + j] = U[i,j]
                testfunc[i*self.geometry.blocksize +j] = w[i,j]

        ## Compute LHS and RHS
        do_simplifications =  self.settings["do_simplifications"].GetBool()
        self._print(1, "\n- Compute RHS")
        rhs = KratosSympy.Compute_RHS(rv_tot.copy(), testfunc, do_simplifications)

        if not self.is_explicit:
            self._print(1, "\n- Compute LHS")
            lhs = KratosSympy.Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS
        else:
            lhs = None

        return(lhs, rhs)

    def _OutputLHSandRHS(self, lhs, rhs, outstring, subscales_type):
        ## Reading and filling the template file
        self._print(1, "\n- Substituting outstring in {}\n".format(self.settings["template_filename"].GetString()))
        mode = self.settings["mode"].GetString()

        rhs_out = KratosSympy.OutputVector_CollectingFactors(rhs, "rRightHandSideBoundedVector", mode)
        outstring = outstring.replace("//substitute_rhs_{}D_{}".format(self.geometry.ndims, subscales_type), rhs_out)

        if not self.is_explicit:
            lhs_out = KratosSympy.OutputMatrix_CollectingFactors(lhs, "lhs", mode)
            outstring = outstring.replace("//substitute_lhs_{}D_{}".format(self.geometry.ndims, subscales_type), lhs_out)

        ## In the explicit element case the container values are referenced in the cpp to limit the container accesses to one per element
        if self.is_explicit:
            ## Substitute the solution values container accesses
            for i_node in range(self.geometry.nnodes):
                for j_block in range(self.geometry.blocksize):
                    to_substitute = 'U({},{})'.format(i_node, j_block)
                    substituted_value = 'U_{}_{}'.format(i_node, j_block)
                    outstring = outstring.replace(to_substitute, substituted_value)

            ## Substitute the solution values time derivatives container accesses
            for i_node in range(self.geometry.nnodes):
                for j_block in range(self.geometry.blocksize):
                    to_substitute = 'dUdt({},{})'.format(i_node, j_block)
                    substituted_value = 'dUdt_{}_{}'.format(i_node, j_block)
                    outstring = outstring.replace(to_substitute, substituted_value)

            ## Substitute the shape function gradients container accesses
            for i_node in range(self.geometry.nnodes):
                for j_dim in range(self.geometry.ndims):
                    to_substitute = 'DN({},{})'.format(i_node, j_dim)
                    substituted_value = 'DN_DX_{}_{}'.format(i_node, j_dim)
                    outstring = outstring.replace(to_substitute, substituted_value)

            ## Substitute the residuals projection container accesses
            for i_node in range(self.geometry.nnodes):
                for j_block in range(self.geometry.blocksize):
                    to_substitute = 'ResProj({},{})'.format(i_node, j_block)
                    substituted_value = 'ResProj_{}_{}'.format(i_node, j_block)
                    outstring = outstring.replace(to_substitute, substituted_value)

        return outstring

    def Generate(self):
        self._print(1, "\nComputing geometry: {}\n".format(self.settings["geometry"].GetString()))

        dim = self.geometry.ndims
        n_nodes = self.geometry.nnodes
        block_size = self.geometry.blocksize

        self.is_explicit = self.settings["explicit"].GetBool()
        params = FormulationParameters(self.geometry)

        # Unknown fields definition (Used later for the gauss point interpolation)
        U = KratosSympy.DefineMatrix('U', n_nodes, block_size)               # Vector of Unknowns (Density, Velocity[dim], Total Energy)
        ResProj = KratosSympy.DefineMatrix('ResProj', n_nodes, block_size)   # Vector of residuals projection

        if self.is_explicit:
            dUdt = KratosSympy.DefineMatrix('dUdt', n_nodes, block_size)     # Vector of Unknowns time derivatives (Density, Velocity[dim], Total Energy)
            Un = None
            Unn = None
        else:
            dUdt = None
            Un = KratosSympy.DefineMatrix('Un', n_nodes, block_size)         # Vector of Unknowns one step back
            Unn = KratosSympy.DefineMatrix('Unn', n_nodes, block_size)       # Vector of Unknowns two steps back

        # Test functions defintiion
        w = KratosSympy.DefineMatrix('w', n_nodes, block_size)     # Variables field test

        # External terms definition
        forcing_terms = {
            "mass":    KratosSympy.DefineVector('m_ext', n_nodes),        # Mass source term
            "thermal": KratosSympy.DefineVector('r_ext', n_nodes),        # Thermal sink/source term
            "force":   KratosSympy.DefineMatrix('f_ext', n_nodes, dim)    # Forcing term
        }

        # Nodal artificial magnitudes
        sc_nodes = ShockCapturingNodalParameters(self.geometry)

        # Backward differantiation coefficients
        bdf = None if self.is_explicit else [sympy.Symbol('bdf'+ i) for i in range(3)]

        ### Construction of the variational equation
        Ug  = KratosSympy.DefineVector('Ug',block_size) # Dofs vector
        H   = KratosSympy.DefineMatrix('H',block_size, dim) # Gradient of U
        mg  = sympy.Symbol('mg') # Mass source term
        f   = KratosSympy.DefineVector('f', dim) # Body force vector
        rg  = sympy.Symbol('rg') # Thermal source/sink term
        V   = KratosSympy.DefineVector('V',block_size) # Test function
        Q   = KratosSympy.DefineMatrix('Q',block_size, dim) # Gradient of V
        acc = KratosSympy.DefineVector('acc',block_size) # Derivative of Dofs/Time
        G   = KratosSympy.DefineMatrix('G',block_size, dim) # Diffusive Flux matrix
        res_proj = KratosSympy.DefineVector('res_proj',block_size) # Residuals projection for the OSS

        ## Calculate the Gauss point residual
        ## Matrix Computation
        S = generate_source_term.ComputeSourceMatrix(Ug, mg, f, rg, params)
        A = generate_convective_flux.ComputeEulerJacobianMatrix(Ug, params)
        if self.settings["shock_capturing"].GetBool():
            sc_params = ShockCapturingParameters()
            G = generate_diffusive_flux.ComputeDiffusiveFluxWithShockCapturing(Ug, H, params, sc_params)
        else:
            G = generate_diffusive_flux.ComputeDiffusiveFlux(Ug, H, params)
        Tau = generate_stabilization_matrix.ComputeStabilizationMatrix(params)

        ## Non-linear operator definition
        self._print(1, "\nCompute non-linear operator\n")
        L = self._ComputeNonLinearOperator(A, H, S, Ug)

        ## FE residuals definition
        # Note that we include the DOF time derivatives in both the implicit and the explicit cases
        # It is required to include it in both cases to calculate the subscale inertial component
        # In the implicit case it is computed with the BDF formulas
        # In the explicit case it is linearised by using the values already stored in the database
        res = - acc - L

        ## Non-linear adjoint operator definition
        self._print(1, "\nCompute non-linear adjoint operator\n")
        L_adj = self._ComputeNonLinearAdjointOperator(A, H, Q, S, Ug, V)

        ## Variational Formulation - Final equation
        (rv, subscales) = self._ComputeVariationalFormulation(A, acc, G, H, L_adj, Q, S, Ug, V)

        #### OSS Residual projections calculation ####
        # Calculate the residuals projection
        self._print(1, "\nCalculate the projections of the residuals")
        projections = {
            "rho"      : KratosSympy.Matrix(KratosSympy.zeros(n_nodes,1)),
            "momentum" : KratosSympy.Matrix(KratosSympy.zeros(n_nodes*dim,1)),
            "energy"   : KratosSympy.Matrix(KratosSympy.zeros(n_nodes,1))
        }
        for i_gauss in range(self.geometry.ngauss):
            self._print(1, "\tGauss point: " + str(i_gauss))
            self._ComputeProjectionsAtGaussPoint(acc, bdf, dUdt, f, forcing_terms, H, i_gauss, mg, projections, res, rg, U, Ug, Un, Unn)

        ## Output the projections
        self.outstring = self._OutputProjections(self.outstring, *projections.values())

        #### Algebraic form calculation ####
        for subscales_type in self._SubscalesTypes():
            ### Substitution of the discretized values at the gauss points
            ## Loop and accumulate the residual in each Gauss point
            rv_tot = KratosSympy.Matrix(KratosSympy.zeros(1,1))

            self._print(1, "\nSubscales type: " + subscales_type)
            self._print(1, "\n- Substitution of the discretized values at the gauss points")
            for i_gauss in range(self.geometry.ngauss):
                ## Substitute the subscales model
                rv_gauss = self._SubstituteSubscales(res, res_proj, rv, subscales, subscales_type, Tau)
                rv_tot += self._ComputeResidualAtGaussPoint(acc, bdf, dUdt, f, forcing_terms, H, i_gauss, mg, params, Q, res_proj, ResProj, rg, rv_gauss, sc_nodes, sc_params, subscales_type, Tau, U, Ug, Un, Unn, V, w)


            (lhs, rhs) = self._ComputeLHSandRHS(rv_tot, U, w)
            self.outstring = self._OutputLHSandRHS(lhs, rhs, self.outstring, subscales_type)

    def Write(self):
        filename = self.settings["output_filename"].GetString()
        self._print(1, "\nWriting {}\n".format(filename))

        self.outputfile.write(self.outstring)
        self.outputfile.close()

        self._print(1, "{} generated\n".format(filename))


if __name__ == "__main__":
    parameters = KratosMultiphysics.Parameters("""
    {
        "2D" :
        {
            "geometry": "triangle",
            "template_filename" : "compressible_navier_stokes_explicit_cpp_template_with_integration.cpp",
            "output_filename"   : "compressible_explicit_navier_stokes.cpp"
        },
        "3D" :
        {
            "geometry": "tetrahedron",
            "template_filename" : "compressible_explicit_navier_stokes.cpp",
            "output_filename"   : "compressible_explicit_navier_stokes.cpp"
        }
    }
    """)

    path = pathlib.Path(__file__).parent
    os.chdir(path)

    generator_2d = SymbolicGenerator(parameters["2D"])
    generator_2d.Generate()
    generator_2d.Write()

    generator_3d = SymbolicGenerator(parameters["3D"])
    generator_3d.Generate()
    generator_3d.Write()
