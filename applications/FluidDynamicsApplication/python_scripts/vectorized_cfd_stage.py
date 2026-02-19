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

        self.dt = settings["time_stepping"]["time_step"].GetDouble()
        self.model_part.ProcessInfo.SetValue(KM.DELTA_TIME, self.dt)


    def ComputeLumpedMass(self):
        Mscalar = np.empty((len(self.model_part.Nodes)))
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
        vavg=np.sum(v_elemental,axis=1) / (self.dim+1.0)
        vnorms = np.linalg.norm(vavg,axis=1)
        tau = 1./(rho/dt + vnorms/h + viscosity/h**2)
        return tau

    def PrepareModelPart(self):
        super().PrepareModelPart()
        self.dim = 2 ###TODO: change this for 3D
        self.model_part.ProcessInfo[KM.STEP] = 0

        vol_geometries = []
        for g in self.model_part.Geometries:
            if(len(g) == 3):
                vol_geometries.append(g.Id)

        self.vol_mp = self.model_part.CreateSubModelPart("fluid_computational_model_part")
        self.vol_mp.AddGeometries(vol_geometries)

        ##TODO: change viscosity
        self.nu = 0.0001
        self.rho = 1.0


    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.REACTION)
        self.model_part.AddNodalSolutionStepVariable(KM.REACTION_WATER_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        self.model_part.AddNodalSolutionStepVariable(KM.BODY_FORCE)
        self.model_part.AddNodalSolutionStepVariable(KM.EXTERNAL_PRESSURE)

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
        self.dt = self.model_part.ProcessInfo[KM.DELTA_TIME]
        self.time += self.dt
        self.model_part.CloneTimeStep(self.time)
        self.model_part.ProcessInfo[KM.STEP] += 1
        return self.time


    def Initialize(self):
        super().Initialize()

        self.nnodes = len(self.model_part.Nodes)

        self.v_adaptor_n = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.VELOCITY,data_shape=[self.dim],step_index=1)
        self.v_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.VELOCITY,data_shape=[self.dim],step_index=0)
        self.p_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.PRESSURE,0)

        self.v_adaptor.CollectData()
        v = self.v_adaptor.data.reshape((len(self.model_part.Nodes),self.dim))

        self.p_adaptor.CollectData()
        p = self.p_adaptor.data

        self.connectivity_adaptor = KM.TensorAdaptors.ConnectivityIdsTensorAdaptor(self.vol_mp.Geometries)
        print("******************",self.connectivity_adaptor.data.dtype)
        self.connectivity_adaptor.CollectData()
        self.connectivity = self.connectivity_adaptor.data
        self.connectivity -= 1 #have indices to start in
        self.connectivity = self.connectivity.astype(np.int64) #############OUCH!!! TODO: this by defaults should be probably int64
        # max_idx = self.connectivity.max()
        # if max_idx <= np.iinfo(np.int32).max: #use 32 bit integers for connectivities
        #     self.connectivity = self.connectivity.astype(np.int32)
        # else:
        #     print(f"Warning: Max index {max_idx} is too large for uint32. Staying at 64-bit.")

        #preallocation of elemental arrays of velocities and pressures
        self.vec_elemental_data = np.empty((*self.connectivity.shape, v.shape[1]))
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
        self.v_gauss  = np.empty((nelem, self.dim))      # velocity at Gauss point
        self.aux_res_el = np.empty((nelem, self.n_in_el, self.dim))
        self.res        = np.empty((nelem, 3, self.dim))
        self.Projection         = np.empty((self.nnodes, self.dim))
        self.Proj_el = np.empty((nelem, 3, self.dim))

        # SCRATCH ARRAYS for memory reuse
        self.elemental_scalar_scratch = np.empty((nelem, self.n_in_el))
        self.nodal_vector_scratch = np.empty((self.nnodes, self.dim))
        self.elemental_grad_scratch = np.empty(self.DN.shape)
        self.Lel_scratch = np.empty((nelem, self.n_in_el, self.n_in_el))
        self.rhs_el_scratch = np.empty((nelem, self.n_in_el))

        ########################## REMOVE - just for testing
        for node in self.model_part.Nodes:
            if(node.X < 0.0001):
                node.Fix(KM.VELOCITY_X)
                node.SetSolutionStepValue(KM.VELOCITY_X,0,1.0)

        ###allocate the graph of the matrix
        self.L,self.L_assembly_indices = self.cfd_utils.AllocateScalarMatrix(self.connectivity)

    def ComputeVelocityProjection(self,v_elemental, Pi_conv):
        #solve (w,pi) = (w,a·∇u)
        convective=self.aux_res_el #rename this to make it more readable
        self.cfd_utils.ComputeElementalGradient(self.DN, v_elemental, self.grad_v)
        self.cfd_utils.InterpolateValue(self.N,v_elemental,self.v_gauss)
        self.cfd_utils.ComputeConvectiveContribution(self.N, self.grad_v, self.v_gauss, convective)
        convective *= self.elemental_volumes[:,np.newaxis,np.newaxis]
        self.cfd_utils.AssembleVector(self.connectivity, convective, Pi_conv)
        Pi_conv = Pi_conv.ravel()
        Pi_conv *= self.Minv
        Pi_conv = Pi_conv.reshape((self.nnodes, self.dim))


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

    def SolveStep1(self,v,vold,p,dt):
        ##TODO: a lot of intermediate storage could be avoided reusing vectors
        self.ElemData(p, self.connectivity, self.scalar_elemental_data)
        pel = self.scalar_elemental_data #no copy ... just renaming for better readability

        self.ElemData(vold, self.connectivity, self.vec_elemental_data)
        vel = self.vec_elemental_data #no copy ... just renaming for better readability

        self.tau = self.ComputeTau(self.h,vel,self.nu,dt,self.rho)

        #########################################################
        #compute convective projection
        self.ComputeVelocityProjection(vel, self.Projection) #TODO: actually the vector Projection is not needed after getting the ElemData
        proj_el = self.Proj_el
        self.ElemData(self.Projection, self.connectivity, proj_el)
        proj_el = np.zeros(vel.shape)



        #########################################################
        #advance in time by runge kutta

        # --- k1 ---
        k1 = self.ComputeVelocityResidual(vel,pel,proj_el,self.DN,self.tau)
        k1.reshape(-1)[:] *= self.Minv

        # --- k2 ---
        v2 = vold + 0.5 * dt * k1
        self.ApplyVelocityDirichletConditions(vold,self.v_end_of_step_dirichlet,0.5,v2)

        self.ElemData(v2, self.connectivity, vel)
        k2 = self.ComputeVelocityResidual(vel, pel,proj_el, self.DN,self.tau)
        k2.reshape(-1)[:] *= self.Minv

        # --- k3 ---
        v3 = vold + 0.5 * dt * k2
        self.ApplyVelocityDirichletConditions(vold,self.v_end_of_step_dirichlet,0.5,v3)

        self.ElemData(v3, self.connectivity, vel)
        k3 = self.ComputeVelocityResidual(vel, pel,proj_el, self.DN,self.tau)
        k3.reshape(-1)[:] *= self.Minv

        # --- k4 ---
        v4 = vold + dt * k3
        self.ApplyVelocityDirichletConditions(vold,self.v_end_of_step_dirichlet,1.0,v4)

        self.ElemData(v4, self.connectivity, vel)
        k4 = self.ComputeVelocityResidual(vel, pel,proj_el, self.DN,self.tau)
        k4.reshape(-1)[:] *= self.Minv

        # --- final RK4 update ---
        vnew = vold + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
        self.ApplyVelocityDirichletConditions(vold, self.v_end_of_step_dirichlet,1.0,vnew)

        return vnew

    def SolveStep2(self,vfrac,p,dt):
        nelem = self.DN.shape[0]
        n_in_el = self.DN.shape[1]
        dim = self.DN.shape[2]


        ##allocation of intermediate arrays - reusing scratch
        Lel = self.Lel_scratch
        aux_scalar = self.elemental_scalar_scratch # Use this for N_div, Lp, stabilized term (SCALAR)
        RHSel = self.rhs_el_scratch
        Pi_press_nodal = self.nodal_vector_scratch
        #tmp = self.elemental_vector_scratch  <-- removed
        Pi_press_el = self.elemental_grad_scratch # REUSE GRAD SCRATCH

        ###############################################################
        #gather data from the database
        self.ElemData(p, self.connectivity, self.scalar_elemental_data)
        pel = self.scalar_elemental_data #no copy ... just renaming for better readability

        self.ElemData(vfrac, self.connectivity, self.vec_elemental_data)
        vel = self.vec_elemental_data #no copy ... just renaming for better readability

        ###############################################################
        ##compute pressure proj and store it into Pi_press_nodal
        self.cfd_utils.Compute_N_DN(self.N, self.DN, pel, Pi_press_el)
        Pi_press_el *= self.elemental_volumes[:,np.newaxis,np.newaxis]
        Pi_press_nodal.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity,Pi_press_el,Pi_press_nodal)
        Pi_press_nodal = Pi_press_nodal.ravel()
        Pi_press_nodal *= self.Minv
        Pi_press_nodal = Pi_press_nodal.reshape((self.nnodes, self.dim))

        ###############################################################
        #from now on compute the LHS and RHS of the pressure problem to be solved
        ###############################################################

        #obtain the nodal versio nof Pi_press_el
        self.ElemData(vfrac, self.connectivity, Pi_press_el)

        # the L_IJ := (∇N_I,∇N_J)
        self.cfd_utils.ComputeLaplacianMatrix(self.DN,Lel)

        # (q,∇·u)
        self.cfd_utils.ComputeElementwiseNodalDivergence(self.N, self.DN, vel, aux_scalar)

        # RHSel = -aux_scalar (in-place)
        np.negative(aux_scalar, out=RHSel)

        # dt*(∇q,∇pold)=dt*Lij*pj
        # Lp = aux_scalar #reuse
        np.einsum("eIJ,eJ->eI",Lel,pel,out=aux_scalar)
        aux_scalar *= dt
        RHSel += aux_scalar

        # tau*(∇q,Pi_pressure)
        self.cfd_utils.ComputePressureStabilization_ProjectionTerm(self.N, self.DN, Pi_press_el,aux_scalar) #writes to aux_scalar
        aux_scalar *= self.tau[:,np.newaxis]
        RHSel += aux_scalar #TODO: check carefully this term

        #multiplying by the elemental volume
        Lel *= (self.dt+self.tau[:,np.newaxis,np.newaxis])*self.elemental_volumes[:,np.newaxis,np.newaxis]
        RHSel *= self.elemental_volumes[:,np.newaxis]

        #assemble RHS
        RHSpressure = np.zeros(vfrac.shape[0])
        self.cfd_utils.AssembleVector(self.connectivity,RHSel,RHSpressure)
        RHSpressure *= self.rho

        #assemble LHS
        self.cfd_utils.AssembleScalarMatrixByCSRIndices(Lel,self.L_assembly_indices,self.L)

        # dt*(∇q,∇pold)
        ##RHSpressure += (self.L@p) #TODO: use SpMV function to avoid temporaries



        #apply Pressure BCs
        # self.pressure_free_dofs_vector = np.ones(self.nnodes,np.float64)
        # self.pressure_free_dofs_vector[self.fix_press_indices] = 0.0
        # self.L.ApplyHomogeneousDirichlet(self.pressure_free_dofs_vector, 1e9, RHSpressure)
        self.cfd_utils.ApplyHomogeneousDirichlet(self.fix_press_indices, self.L, RHSpressure)

        #solve the system
        ##TODO: use amgcl instead of this one!! ... and also avoid creating temporaries
        Ascipy = KM.scipy_conversion_tools.to_csr(self.L)
        p = scipy.sparse.linalg.spsolve(Ascipy, RHSpressure)

        return p


    def SolveStep3(self,v,vfrac,delta_p,dt):

        grad_dp_el = self.elemental_grad_scratch

        self.ElemData(delta_p, self.connectivity, self.scalar_elemental_data)
        delta_pel = self.scalar_elemental_data #no copy ... just renaming for better readability

        self.cfd_utils.Compute_DN_N(self.N, self.DN, delta_pel, grad_dp_el)
        grad_dp_el *= self.elemental_volumes[:,np.newaxis,np.newaxis]

        grad_dp_nodal = self.nodal_vector_scratch
        grad_dp_nodal.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity,grad_dp_el,grad_dp_nodal)

        v.ravel()[:] = vfrac.ravel() - (dt/self.rho)*  self.Minv * grad_dp_nodal.ravel()
        self.ApplyVelocityDirichletConditions(vfrac,self.v_end_of_step_dirichlet,1.0,v)
        return v

    def SolveSolutionStep(self):
        self.GetFixity() #TODO: do this only once!

        #obtain velocity from Kratos
        self.v_adaptor.CollectData() #new vel
        v = self.v_adaptor.data.reshape((len(self.model_part.Nodes),self.dim))

        self.v_adaptor_n.CollectData() #old
        vold = self.v_adaptor_n.data.reshape((len(self.model_part.Nodes),self.dim))

        self.p_adaptor.CollectData()
        p = self.p_adaptor.data

        dt = self.model_part.ProcessInfo[KM.DELTA_TIME]

        #get dirichlet values
        self.v_end_of_step_dirichlet = v.ravel()[self.fix_vel_indices].copy()


        #perform fractional step
        pold = p.copy()
        vfrac = self.SolveStep1(v,vold,pold,dt)
        p = self.SolveStep2(vfrac,pold,dt)
        delta_p = p - pold #TODO: most probably we could modify p in place
        v = self.SolveStep3(v,vfrac,delta_p,dt)

        #update Kratos
        self.v_adaptor.data = v
        self.v_adaptor.StoreData()
        self.p_adaptor.data = p
        self.p_adaptor.StoreData()
        #value = input("Press Enter to continue...")


    def ComputeVelocityResidual(self, v_elemental, p_elemental, proj_el, DN, tau):
        N = self.N

        # assign a shorter local name
        grad_v   = self.grad_v
        v_gauss  = self.v_gauss
        aux_res_el = self.aux_res_el
        res        = self.res

        #(w,b)
        res.fill(0.0) #TODO: substitute by the proper computation of the external forces

        #(w,a·∇u) + (a·∇w,a·∇u) - (a·∇w,Pi_conv)
        convective=aux_res_el #rename this to make it more readable
        convective_stabilization = aux_res_el #rename this to make it more readable #TODO: most probably we can get rid of this array
        self.cfd_utils.ComputeElementalGradient(DN, v_elemental, grad_v)

        tmp = np.empty_like(convective)
        tmp_stab = np.empty_like(convective_stabilization)
        n_gauss_convective = self.dim + 1
        w_gauss_convective = self.cfd_utils.GetGaussIntegrationWeights(self.dim, 2)
        N_gauss_convective = self.cfd_utils.GetShapeFunctionsOnGaussPoints(self.dim, 2)
        for i_gauss in range(n_gauss_convective):
            # Get elemental velocities at current integration point
            self.cfd_utils.InterpolateValue(N_gauss_convective[i_gauss,:], v_elemental, v_gauss)

            # Compute convective contribution at current integration point (w,a·∇u)
            self.cfd_utils.ComputeConvectiveContribution(N_gauss_convective[i_gauss,:], grad_v, v_gauss, tmp)

            #TODO: Compute convective stabilization contribution at current integration point (a·∇w,a·∇u) - (a·∇w,Pi_conv)
            #proj_el = np.zeros(v_elemental.shape) #TODO: this has to be computed in the appropriate way!
            self.cfd_utils.ComputeMomentumStabilization(N_gauss_convective[i_gauss,:], DN, v_gauss, v_elemental, proj_el, tmp_stab, self.elemental_scalar_scratch, self.elemental_grad_scratch)

            # Scale with current integration weight (Jacobian determinant is applied at the end to the entire residual)
            convective += w_gauss_convective[i_gauss] * tmp
            convective_stabilization += w_gauss_convective[i_gauss] * tmp_stab

        convective *= self.rho
        convective_stabilization *= self.tau[:,np.newaxis,np.newaxis] #TODO: activate stabilization
        res -= convective
        res -= convective_stabilization

        #(∇w,∇u)
        viscous=aux_res_el #rename this to make it more readable
        self.cfd_utils.ApplyLaplacian(DN,v_elemental,viscous)
        viscous *= self.nu #TODO, check if we use nu or mu
        res -= viscous

        # (∇w,p_gauss)
        pressure_term=aux_res_el #rename this to make it more readable
        self.cfd_utils.Compute_DN_N(N,DN,p_elemental,pressure_term)
        res += pressure_term

        #multiply by volume
        res *= self.elemental_volumes[:,np.newaxis,np.newaxis]

        #VECTORIZED FEM ASSEMBLY
        Rassembled = np.zeros((self.nnodes,self.dim),dtype=res.dtype) ##TODO preallocate this
        self.cfd_utils.AssembleVector(self.connectivity,res,Rassembled)
        return Rassembled

    def ImportModelPart(self):
        pass

    def Run(self):
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    simulation = Stage(model,parameters)
    simulation.Run()