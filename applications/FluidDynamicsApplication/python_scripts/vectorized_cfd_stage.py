import KratosMultiphysics as KM
import KratosMultiphysics.analysis_stage as analysis_stage
import numpy as np
import KratosMultiphysics.FluidDynamicsApplication.cfd_utils as cfd_utils
import scipy.sparse.linalg
import KratosMultiphysics.scipy_conversion_tools

class VectorizedCFDStage(analysis_stage.AnalysisStage):
    def __init__(self, model, project_parameters):
        super().__init__(model,project_parameters)

        self.cfd_utils = cfd_utils.CFDUtils()

        # Either retrieve the model part from the model or create a new one
        settings = self.project_parameters["solver_settings"]
        model_part_name = settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            self.model_part = self.model.CreateModelPart(model_part_name)

        self.dim = settings["domain_size"].GetInt()
        if(self.dim==2):
            self.n_in_el = 3
        else:
            self.n_in_el = 4

        if self.dim == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')

        self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, self.dim)

        self.AddVariables()

        self.automatic_dt = settings["time_stepping"]["automatic_time_step"].GetBool()
        if self.automatic_dt:
            self.target_cfl = settings["time_stepping"]["CFL_number"].GetDouble() # CFL number from which the DT will be computed
            self.target_fourier = settings["time_stepping"]["Viscous_Fourier_number"].GetDouble() # Viscous Fourier number from which the DT will be computed
            self.dt_max = settings["time_stepping"]["maximum_delta_time"].GetDouble() # Maximum value of the automatically computed time step
        self.dt = settings["time_stepping"]["time_step"].GetDouble() if not self.automatic_dt else None # Initialize dt to None if automatic, to be computed later

        self.convection_integration_order = 2 # TODO: think on exposing this to the json

    def ComputeLumpedMass(self):
        Mscalar = np.zeros((len(self.model_part.Nodes)))
        Mel = np.einsum("e,i->ei",self.elemental_volumes,self.N)
        self.cfd_utils.AssembleVector(self.connectivity,Mel,Mscalar) #here i have one mass per node, but we will need one per component
        M = np.empty((Mscalar.shape[0],self.dim))
        for i in range(0,self.dim):
            M[:,i] = Mscalar[:]
        M = M.ravel()
        return M

    def ApplyVelocityDirichletConditions(self, vold, v_end_of_step, dt_factor, v):
        v.ravel()[self.fix_vel_indices] = (1-dt_factor)*vold.ravel()[self.fix_vel_indices] + dt_factor*v_end_of_step

    def ComputeH(self,DN):
        inverse_h2 = np.max(np.linalg.norm(DN,axis=2),axis=1)
        h = np.sqrt(1./inverse_h2)
        return h

    def ComputeTau(self,h,v_elemental,viscosity,dt,rho):
        c_1 = 4.0
        c_2 = 2.0
        vavg=np.sum(v_elemental,axis=1) / (self.dim+1.0)
        vnorms = np.linalg.norm(vavg,axis=1)
        tau_1 = 1.0/(rho/dt + c_2*rho*vnorms/h + c_1*viscosity/h**2)
        # tau_2 = h**2/(c_1 * tau_1)
        return tau_1
        # return tau_1, tau_2

    def PrepareModelPart(self):
        super().PrepareModelPart()
        self.model_part.ProcessInfo[KM.STEP] = 0

        vol_geometries = []
        for g in self.model_part.Geometries:
            if(len(g) == 3):
                vol_geometries.append(g.Id)

        self.vol_mp = self.model_part.CreateSubModelPart("fluid_computational_model_part")
        self.vol_mp.AddGeometries(vol_geometries)

        ##TODO: change viscosity
        self.nu = 1.0e-1
        self.rho = 1.0e3
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KM.BODY_FORCE, 0, np.array([0.0,-10.0, 0.0]))

    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.REACTION)
        self.model_part.AddNodalSolutionStepVariable(KM.REACTION_WATER_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        self.model_part.AddNodalSolutionStepVariable(KM.BODY_FORCE)
        self.model_part.AddNodalSolutionStepVariable(KM.EXTERNAL_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.CONV_PROJ) #FIXME: remove after debugging
        self.model_part.AddNodalSolutionStepVariable(KM.DIVPROJ) #FIXME: remove after debugging
        self.model_part.AddNodalSolutionStepVariable(KM.ACCELERATION) #FIXME: remove after debugging
        self.model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT) #FIXME: remove after debugging

        self.model_part.SetBufferSize(2)

    def AddDofs(self):
        dofs_and_reactions_to_add = []
        dofs_and_reactions_to_add.append(["VELOCITY_X", "REACTION_X"])
        dofs_and_reactions_to_add.append(["VELOCITY_Y", "REACTION_Y"])
        dofs_and_reactions_to_add.append(["VELOCITY_Z", "REACTION_Z"])
        dofs_and_reactions_to_add.append(["PRESSURE", "REACTION_WATER_PRESSURE"])
        KM.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.model_part)

        KM.Logger.PrintInfo(self.__class__.__name__, "Fluid solver DOFs added correctly.")

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super()._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["gravity", "initial_conditions_process_list", "boundary_conditions_process_list", "auxiliar_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format (or no processes are specified)
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        info_msg  = "Using the old way to create the processes, this was removed!\n"
                        info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                        info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                        info_msg += "for a description of the new format"
                        raise Exception("FluidDynamicsAnalysis: " + info_msg)
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not allowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                info_msg  = "Using the old way to create the gid-output, this was removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                raise Exception("FluidDynamicsAnalysis: " + info_msg)
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _GetOrderOfProcessesInitialization(self):
        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        return ["gravity",
                "initial_conditions_process_list",
                "boundary_conditions_process_list",
                "auxiliar_process_list"]

    def _GetSimulationName(self):
        return "Fluid Dynamics Analysis"


    def ElemData(self,v,connectivities,out):
        expected_shape = (*connectivities.shape, *v.shape[1:])
        if out.shape != expected_shape:
            raise ValueError(f"Shape of 'out' {out.shape} does not match expected shape {expected_shape}")
        np.take(v, connectivities,axis=0, out=out)
        return out

    def AdvanceInTime(self):
        dt = self._ComputeDeltaTime()
        self.time += dt

        self.model_part.CloneTimeStep(self.time)
        self.model_part.ProcessInfo[KM.STEP] += 1

        return self.time

    def Initialize(self):
        super().Initialize()

        self.nnodes = len(self.model_part.Nodes)

        self.v_adaptor_n = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.VELOCITY,data_shape=[self.dim],step_index=1)
        self.v_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.VELOCITY,data_shape=[self.dim],step_index=0)
        self.p_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.PRESSURE,0)
        self.p_adaptor_n = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, KM.PRESSURE, 1) #TODO: most probably this is not required and we can just copy current step data (but just in case while debugging)
        self.b_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, KM.BODY_FORCE, data_shape=[self.dim], step_index=0)

        self.v_adaptor.CollectData()
        v = self.v_adaptor.data.reshape((len(self.model_part.Nodes),self.dim))

        self.p_adaptor.CollectData()
        p = self.p_adaptor.data

        self.connectivity_adaptor = KM.TensorAdaptors.ConnectivityIdsTensorAdaptor(self.vol_mp.Geometries)
        self.connectivity_adaptor.CollectData()
        self.connectivity = self.connectivity_adaptor.data
        self.connectivity -= 1 #have indices to start in
        self.connectivity = self.connectivity.astype(np.int64) #############OUCH!!! TODO: this by defaults should be probably int64
        # max_idx = self.connectivity.max()
        # if max_idx <= np.iinfo(np.int32).max: #use 32 bit integers for connectivities
        #     self.connectivity = self.connectivity.astype(np.int32)
        # else:
        #     print(f"Warning: Max index {max_idx} is too large for uint32. Staying at 64-bit.")

        #preallocation of elemental arrays of velocities, pressures and body force
        self.vec_elemental_data = np.empty((*self.connectivity.shape, v.shape[1]))
        self.vec_elemental_data_b = np.empty((*self.connectivity.shape, v.shape[1]))
        self.scalar_elemental_data = np.empty(self.connectivity.shape)

        #preallocation of local array of shape functions
        self.N = np.ones(self.dim+1)/(self.dim+1)

        #obtain the shape function derivatives
        geometry_adaptor_DN = KM.TensorAdaptors.GeometriesTensorAdaptor(
            self.vol_mp.Geometries,
            KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.ShapeFunctionDerivatives,
            KM.GeometryData.IntegrationMethod.GI_GAUSS_1)
        geometry_adaptor_DN.CollectData()

        self.DN = np.squeeze(geometry_adaptor_DN.data).copy() #this has shape nel*1*nnodes_in_el*dim - the copy is important as we need to own the data

        #obtain elemental volumes #TODO: an adaptor should be available for these
        vols = []
        for geom in self.vol_mp.Geometries:
            vols.append(geom.DomainSize())
        self.elemental_volumes = np.array(vols)

        #compute h per every element
        self.h = self.ComputeH(self.DN)

        #compute lumped mass matrix and its inverse
        self.M=self.ComputeLumpedMass() #TODO: decide if we need to store M
        self.Minv = 1./self.M

        # Preallocate all outputs
        nelem = self.DN.shape[0]
        self.grad_v   = np.empty((nelem, self.dim, self.dim))   # gradient of velocity
        # self.div_v = np.empty((nelem))                            # divergence of velocity
        # self.proj_div_v = np.empty((nelem))                            # projeciton of the velocity divergence
        self.v_gauss  = np.empty((nelem, self.dim))      # velocity at Gauss point
        self.b_gauss  = np.empty((nelem, self.dim))      # body force at Gauss point
        self.aux_res_el = np.empty((nelem, self.n_in_el, self.dim))
        self.res        = np.empty((nelem, 3, self.dim))
        self.Projection         = np.empty((self.nnodes, self.dim))
        self.Proj_el = np.empty((nelem, 3, self.dim))
        self.div_proj_el = np.empty((nelem, self.n_in_el))
        self.div_proj = np.empty(self.nnodes)

        # SCRATCH ARRAYS for memory reuse
        self.nodal_vector_scratch = np.empty((self.nnodes, self.dim))
        self.elemental_scalar_scratch = np.empty((nelem, self.n_in_el))
        self.nodal_vector_scratch = np.empty((self.nnodes, self.dim))
        self.elemental_grad_scratch = np.empty(self.DN.shape)
        self.Lel_scratch = np.empty((nelem, self.n_in_el, self.n_in_el))
        self.rhs_el_scratch = np.empty((nelem, self.n_in_el))

        ###allocate the graph of the matrix
        self.L, self.L_assembly_indices = self.cfd_utils.AllocateScalarMatrix(self.connectivity)

        # Compute the Fourier number time step restriction (constant as h, rho and mu are constant in time)
        if self.automatic_dt:
            aux = self.target_fourier * self.rho / self.nu
            self.dt_fourier = np.min(aux * self.h**2)

        #FIXME: initial conditions for debuggings
        #TODO: hydrostatic pressure initialization for debugging
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, (2.0-node.Y)*10000) #TODO: remove after debugging
            node.SetSolutionStepValue(KM.PRESSURE, 1, (2.0-node.Y)*10000) #TODO: remove after debugging
            # node.Fix(KM.PRESSURE) #TODO: remove after debugging

    def ComputeVelocityProjection(self, v_elemental, Pi_conv):
        # Solve (w,pi) = (w,rho·a·∇u)
        # Calculate the convective term at the elemental level
        convective = self.aux_res_el # rename this to make it more readable
        convective.fill(0.0) # initialize this to zero as it will be incremented in the integration points loop below
        tmp = np.zeros_like(convective)
        self.cfd_utils.ComputeElementalGradient(self.DN, v_elemental, self.grad_v)
        det_J_volume_factor = 2.0 if self.dim == 2 else 6.0
        w_gauss_convective = self.cfd_utils.GetGaussIntegrationWeights(self.dim, self.convection_integration_order)
        N_gauss_convective = self.cfd_utils.GetShapeFunctionsOnGaussPoints(self.dim, self.convection_integration_order)
        for i_gauss in range(w_gauss_convective.size):
            # Get elemental velocities at current integration point
            self.cfd_utils.InterpolateValue(N_gauss_convective[i_gauss,:], v_elemental, self.v_gauss)

            # Compute convective contribution at current integration point (w,rho·a·∇u)
            self.cfd_utils.ComputeConvectiveContribution(N_gauss_convective[i_gauss,:], self.grad_v, self.v_gauss, self.rho, tmp)

            # Scale with current integration weight (Jacobian determinant is applied at the end to the entire residual)
            convective += det_J_volume_factor * w_gauss_convective[i_gauss] * tmp
        convective *= self.elemental_volumes[:,np.newaxis,np.newaxis]

        # Do the nodal assembly of the projection elemental contributions
        Pi_conv.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity, convective, Pi_conv)

        # Solve with the lumped mass matrix to get the nodal values of the projection
        Pi_conv = Pi_conv.ravel()
        Pi_conv *= self.Minv
        Pi_conv = Pi_conv.reshape((self.nnodes, self.dim))

    def ComputeDivergenceProjection(self, v_elemental, Pi_div):
        # Calculate the divergence term at the elemental level
        aux_scalar = self.elemental_scalar_scratch # Use this for N_div, Lp, stabilized term (SCALAR)
        self.cfd_utils.ComputeElementwiseNodalDivergence(self.N, self.DN, v_elemental, aux_scalar)
        aux_scalar *= self.elemental_volumes[:,np.newaxis]

        # Do the nodal assembly of the projection elemental contributions
        Pi_div.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity, aux_scalar, Pi_div)

        # Solve with the lumped mass matrix to get the nodal values of the projection
        Pi_div = Pi_div.ravel()
        Pi_div *= self.Minv[::self.dim] # Strided view to take one value every self.dim
        Pi_div = Pi_div.reshape((self.nnodes, 1))

    def GetFixity(self):
        #obtain fixity
        #TODO: this should be done with an adaptor
        self.fix_vel_indices = []
        self.fix_press_indices = []
        for node in self.model_part.Nodes:
            if node.IsFixed(KM.VELOCITY_X):
                self.fix_vel_indices.append((node.Id-1)*self.dim)
            if node.IsFixed(KM.VELOCITY_Y):
                self.fix_vel_indices.append((node.Id-1)*self.dim+1)
            if self.dim==3 and node.IsFixed(KM.VELOCITY_Z):
                self.fix_vel_indices.append((node.Id-1)*self.dim+2)
            if node.IsFixed(KM.PRESSURE):
                self.fix_press_indices.append(node.Id-1)


        self.fix_vel_indices = np.array(self.fix_vel_indices, np.int64)
        self.fix_press_indices = np.array(self.fix_press_indices, np.int64)

    def SolveStep1(self,v,vold,p,b,dt):
        ##TODO: a lot of intermediate storage could be avoided reusing vectors
        self.ElemData(p, self.connectivity, self.scalar_elemental_data)
        pel = self.scalar_elemental_data #no copy ... just renaming for better readability

        self.ElemData(vold, self.connectivity, self.vec_elemental_data)
        vel = self.vec_elemental_data #no copy ... just renaming for better readability

        #FIXME: This is only valid for b constant in time...
        self.ElemData(b, self.connectivity, self.vec_elemental_data_b)
        b_el = self.vec_elemental_data_b #no copy ... just renaming for better readability

        self.tau_1 = self.ComputeTau(self.h, vel, self.nu, dt, self.rho)
        # self.tau_1, self.tau_2 = self.ComputeTau(self.h, vel, self.nu, dt, self.rho)

        #########################################################

        #########################################################
        #compute convective projection
        self.ComputeVelocityProjection(vel, self.Projection) #TODO: actually the vector Projection is not needed after getting the ElemData
        proj_el = self.Proj_el
        self.ElemData(self.Projection, self.connectivity, proj_el)
        #FIXME: remove after debugging
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KM.CONV_PROJ, 0, np.hstack([self.Projection[node.Id-1,:], 0.0])) #TODO: remove after debugging

        self.ComputeDivergenceProjection(vel, self.div_proj) #TODO: actually the vector Projection is not needed after getting the ElemData
        div_proj_el = self.div_proj_el
        self.ElemData(self.div_proj, self.connectivity, div_proj_el)
        #FIXME: remove after debugging
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KM.DIVPROJ, 0, self.div_proj[node.Id-1]) #TODO: remove after debugging

        #########################################################
        #advance in time by runge kutta

        # --- k1 ---
        k1 = self.ComputeVelocityResidual(vel, pel, b_el, proj_el, self.DN, self.tau_1)
        k1.reshape(-1)[:] *= (1.0/self.rho) * self.Minv

        # --- k2 ---
        v2 = vold + 0.5 * dt * k1
        self.ApplyVelocityDirichletConditions(vold,self.v_end_of_step_dirichlet,0.5,v2)

        self.ElemData(v2, self.connectivity, vel)
        k2 = self.ComputeVelocityResidual(vel, pel, b_el, proj_el, self.DN,self.tau_1)
        k2.reshape(-1)[:] *= (1.0/self.rho) * self.Minv

        # --- k3 ---
        v3 = vold + 0.5 * dt * k2
        self.ApplyVelocityDirichletConditions(vold,self.v_end_of_step_dirichlet,0.5,v3)

        self.ElemData(v3, self.connectivity, vel)
        k3 = self.ComputeVelocityResidual(vel, pel, b_el, proj_el, self.DN,self.tau_1)
        k3.reshape(-1)[:] *= (1.0/self.rho) * self.Minv

        # --- k4 ---
        v4 = vold + dt * k3
        self.ApplyVelocityDirichletConditions(vold,self.v_end_of_step_dirichlet,1.0,v4)

        self.ElemData(v4, self.connectivity, vel)
        k4 = self.ComputeVelocityResidual(vel, pel, b_el, proj_el, self.DN,self.tau_1)
        k4.reshape(-1)[:] *= (1.0/self.rho) * self.Minv

        # --- final RK4 update ---
        vnew = vold + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
        self.ApplyVelocityDirichletConditions(vold, self.v_end_of_step_dirichlet,1.0,vnew)

        return vnew

    def SolveStep2(self,vfrac,p,dt):

        ##allocation of intermediate arrays - reusing scratch
        aux_scalar = self.elemental_scalar_scratch # Use this for N_div, Lp, stabilized term (SCALAR)
        L_el = self.Lel_scratch
        rhs_el = self.rhs_el_scratch

        ###############################################################
        #gather data from the database
        self.ElemData(p, self.connectivity, self.scalar_elemental_data)
        pel = self.scalar_elemental_data #no copy ... just renaming for better readability

        self.ElemData(vfrac, self.connectivity, self.vec_elemental_data)
        vel_frac = self.vec_elemental_data #no copy ... just renaming for better readability

        # ###############################################################
        ##compute pressure proj and store it into Pi_press_nodal
        Pi_press_nodal = self.nodal_vector_scratch
        Pi_press_el = self.elemental_grad_scratch # REUSE GRAD SCRATCH

        self.cfd_utils.Compute_N_DN(self.N, self.DN, pel, Pi_press_el)
        Pi_press_el *= self.elemental_volumes[:,np.newaxis,np.newaxis]
        Pi_press_nodal.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity, Pi_press_el, Pi_press_nodal)
        Pi_press_nodal = Pi_press_nodal.ravel()
        Pi_press_nodal *= self.Minv
        Pi_press_nodal = Pi_press_nodal.reshape((self.nnodes, self.dim))
        press_proj_el = self.Proj_el
        self.ElemData(Pi_press_nodal, self.connectivity, press_proj_el)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KM.ACCELERATION, 0, np.hstack([Pi_press_nodal[node.Id-1,:], 0.0])) #TODO: remove after debugging

        ###############################################################
        #from now on compute the LHS and RHS of the pressure problem to be solved
        ###############################################################

        # Assemble pressure LHS (set to zero is done internally)
        self.cfd_utils.ComputeLaplacianMatrix(self.DN, L_el) # elemental laplacian contributions as L_IJ := (∇N_I,∇N_J)
        L_el *= (dt / self.rho) * self.elemental_volumes[:,np.newaxis,np.newaxis] # scale LHS elemental contributions
        self.cfd_utils.AssembleScalarMatrixByCSRIndices(L_el, self.L_assembly_indices, self.L) # assemble the scaled elemental contributions

        # -(q,∇·ufrac)
        self.cfd_utils.ComputeElementwiseNodalDivergence(self.N, self.DN, vel_frac, aux_scalar)
        np.negative(aux_scalar, out=rhs_el) # Note that this acts as an in-place initialization of the RHS

        # (dt/rho + tau)*(∇q,∇pold) = (dt/rho + tau)*Lij*pj
        self.cfd_utils.ApplyLaplacian(self.DN, pel, aux_scalar)
        aux_scalar *= (dt / self.rho) + self.tau_1[:,np.newaxis]
        rhs_el += aux_scalar

        # -tau*(∇q,Pi_pressure)
        self.cfd_utils.ComputePressureStabilization_ProjectionTerm(self.N, self.DN, press_proj_el, aux_scalar) #writes to aux_scalar
        aux_scalar *= self.tau_1[:,np.newaxis]
        rhs_el -= aux_scalar

        # scale RHS elemental contributions by elemental volumes (integration)
        rhs_el *= self.elemental_volumes[:,np.newaxis]

        # Assemble RHS
        rhs = np.zeros(vfrac.shape[0]) # note that initialization is required here as we are using np.add.at
        self.cfd_utils.AssembleVector(self.connectivity, rhs_el, rhs)

        # Apply pressure BCs
        self.cfd_utils.ApplyHomogeneousDirichlet(self.fix_press_indices, self.L, rhs)

        #solve the system
        ##TODO: use amgcl instead of this one!! ... and also avoid creating temporaries
        L_scipy = KM.scipy_conversion_tools.to_csr(self.L)
        p = scipy.sparse.linalg.spsolve(L_scipy, rhs)

        return p

    def SolveStep3(self, v, vfrac, delta_p, dt):

        # Gather elemental pressure increments from the nodal values
        self.ElemData(delta_p, self.connectivity, self.scalar_elemental_data)
        delta_pel = self.scalar_elemental_data #no copy ... just renaming for better readability

        # Calculate the gradient of the pressure increment at the elemental level
        grad_dp_el = self.elemental_grad_scratch
        self.cfd_utils.Compute_DN_N(self.N, self.DN, delta_pel, grad_dp_el)
        grad_dp_el *= self.elemental_volumes[:,np.newaxis,np.newaxis]

        # Assemble the gradient contributions at the nodes
        grad_dp_nodal = self.nodal_vector_scratch
        grad_dp_nodal.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity, grad_dp_el, grad_dp_nodal)

        # Update the fractional velocity with the pressure gradient contribution
        v.ravel()[:] = vfrac.ravel() + (dt/self.rho) * self.Minv * grad_dp_nodal.ravel()
        self.ApplyVelocityDirichletConditions(vfrac, self.v_end_of_step_dirichlet, 1.0, v)

        return v

    def SolveSolutionStep(self):
        self.GetFixity() #TODO: do this only once!

        #obtain velocity from Kratos
        self.v_adaptor.CollectData() #new vel
        v = self.v_adaptor.data.reshape((len(self.model_part.Nodes),self.dim))

        self.v_adaptor_n.CollectData() #old
        vold = self.v_adaptor_n.data.reshape((len(self.model_part.Nodes),self.dim))

        self.p_adaptor_n.CollectData() # Previous step pressure
        pold = self.p_adaptor_n.data

        self.b_adaptor.CollectData() #body force
        b = self.b_adaptor.data.reshape((len(self.model_part.Nodes),self.dim))

        dt = self.model_part.ProcessInfo[KM.DELTA_TIME]

        #get dirichlet values
        self.v_end_of_step_dirichlet = v.ravel()[self.fix_vel_indices].copy()

        #perform fractional step
        vfrac = self.SolveStep1(v, vold, pold, b, dt)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT, 0, np.hstack([vfrac[node.Id-1,:], 0.0]))
        p = self.SolveStep2(vfrac, pold, dt)
        delta_p = p - pold #TODO: most probably we could modify p in place
        v = self.SolveStep3(v, vfrac, delta_p, dt)

        #update Kratos
        self.v_adaptor.data = v
        self.v_adaptor.StoreData()
        self.p_adaptor.data = p
        self.p_adaptor.StoreData()
        #value = input("Press Enter to continue...")


    def ComputeVelocityResidual(self, v_elemental, p_elemental, b_elemental, proj_el, DN, tau):
        N = self.N

        # assign a shorter local name
        grad_v  = self.grad_v
        v_gauss = self.v_gauss
        grad_v  = self.grad_v
        v_gauss = self.v_gauss
        b_gauss = self.b_gauss
        aux_res_el = self.aux_res_el
        res        = self.res

        # Initialize residual
        res.fill(0.0)

        #(w,b) #TODO: we are assuming it constant (put it in the Gauss points loop in the future)
        body_force = aux_res_el
        self.cfd_utils.InterpolateValue(N, b_elemental, b_gauss)
        self.cfd_utils.ComputeBodyForceContribution(N, b_gauss, body_force) #TODO: here we are using the projection as a storage for the body force, which is not ideal
        body_force *= self.rho
        res += body_force

        #-(w,rho·a·∇u) - (rho·a·∇w,tau·rho·a·∇u) + (rho·a·∇w,tau·Pi_conv)
        stabilized_convective = aux_res_el #rename this to make it more readable
        stabilized_convective.fill(0.0) # initialize this to zero as it will be incremented in the integration points loop below
        self.cfd_utils.ComputeElementalGradient(DN, v_elemental, grad_v)

        tmp = np.zeros_like(stabilized_convective)
        tmp_stab = np.zeros_like(stabilized_convective)
        det_J_volume_factor = 2.0 if self.dim == 2 else 6.0
        w_gauss_convective = self.cfd_utils.GetGaussIntegrationWeights(self.dim, self.convection_integration_order)
        N_gauss_convective = self.cfd_utils.GetShapeFunctionsOnGaussPoints(self.dim, self.convection_integration_order)
        for i_gauss in range(w_gauss_convective.size):
            # Get elemental velocities at current integration point
            self.cfd_utils.InterpolateValue(N_gauss_convective[i_gauss,:], v_elemental, v_gauss)

            # Compute convective contribution at current integration point (w,rho·a·∇u)
            self.cfd_utils.ComputeConvectiveContribution(N_gauss_convective[i_gauss,:], grad_v, v_gauss, self.rho, tmp)

            # Compute convective stabilization contribution at current integration point (rho·a·∇w,rho·a·∇u) - (rho·a·∇w,Pi_conv)
            self.cfd_utils.ComputeMomentumStabilization(N_gauss_convective[i_gauss,:], DN, v_gauss, v_elemental, proj_el, self.rho, tmp_stab, self.elemental_scalar_scratch, self.elemental_grad_scratch)

            # Scale with current integration weight (Jacobian determinant is applied at the end to the entire residual)
            aux_w = w_gauss_convective[i_gauss] * det_J_volume_factor
            stabilized_convective -= aux_w * tmp

        res += stabilized_convective

        #-(∇w,∇u)
        viscous = aux_res_el #rename this to make it more readable
        self.cfd_utils.ApplyLaplacian(DN, v_elemental,viscous)
        viscous *= self.nu
        res -= viscous

        #+(∇·w,p_gauss)
        pressure_term = aux_res_el #rename this to make it more readable
        self.cfd_utils.Compute_DN_N(N, DN, p_elemental, pressure_term)
        res += pressure_term

        #-(∇·w,tau_1·∇·v) + (∇·w,tau_1·Pi_div)
        # self.cfd_utils.InterpolateValue(N, proj_div_el, proj_div_v)
        # self.cfd_utils.ComputeElementalDivergence(DN, v_elemental, div_v)
        # self.cfd_utils.ComputeDivDivStabilization(DN, div_v, proj_div_v, self.elemental_scratch, div_div_stab)
        # res -= div_div_stab * self.tau_2[:,np.newaxis,np.newaxis]

        #multiply by volume (note that the convective term already includes the corresponding weighting factor to emulate the Jacobian determinant, so we just need to apply the volume here)
        res *= self.elemental_volumes[:,np.newaxis,np.newaxis]

        #VECTORIZED FEM ASSEMBLY
        Rassembled = np.zeros((self.nnodes,self.dim),dtype=res.dtype) ##TODO preallocate this
        self.cfd_utils.AssembleVector(self.connectivity,res,Rassembled)
        return Rassembled

    def ImportModelPart(self):
        pass

    def GetStep(self):
        return self.model_part.ProcessInfo[KratosMultiphysics.STEP]

    def _ComputeDeltaTime(self):
        if self.automatic_dt:
            # Get previous step velocity data
            self.v_adaptor_n.CollectData() #old
            vold = self.v_adaptor_n.data.reshape((len(self.model_part.Nodes),self.dim))
            self.ElemData(vold, self.connectivity, self.vec_elemental_data)

            # Calculate the velocity at the midpoint of all elements
            vavg=np.sum(self.vec_elemental_data, axis=1) / (self.dim+1.0)
            vavg_norm = np.linalg.norm(vavg, axis=1)

            # Calculate the minimum delta time
            v_zero_mask = vavg_norm > 1.0e-12 # Filter the zero velocity values
            if np.any(v_zero_mask):
                dt_cfl = np.min(self.target_cfl * self.h[v_zero_mask] / vavg_norm[v_zero_mask])
                self.dt = min(min(dt_cfl, self.dt_fourier), self.dt_max)
            else:
                init_dt_factor = 0.05
                self.dt = init_dt_factor * min(self.dt_max, self.dt_fourier)

        if (self.dt is None or self.dt <= 0.0):
            raise Exception("Delta time was not computed properly")

        return self.dt

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    simulation = VectorizedCFDStage(model,parameters)
    simulation.Run()