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
        tau = 1./(dt/rho + vnorms/h + viscosity/h**2)
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
        self.connectivity_adaptor.CollectData()
        self.connectivity = self.connectivity_adaptor.data
        self.connectivity -= 1 #have indices to start in 

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
        self.grad_v   = np.empty((nelem, 2, 2))   # gradient of velocity
        self.v_gauss  = np.empty((nelem, 2))      # velocity at Gauss point
        self.aux_res_el = np.empty((nelem, 3, 2))
        self.res        = np.empty((nelem, 3, 2))

        ########################## REMOVE - just for testing
        for node in self.model_part.Nodes:
            if(node.X < 0.0001):
                node.Fix(KM.VELOCITY_X)
                node.SetSolutionStepValue(KM.VELOCITY_X,0,1.0)

        ###allocate the graph of the matrix
        self.L = self.cfd_utils.AllocateScalarMatrix(self.connectivity)


    def ComputeVelocityProjection(self,v_elemental):
        
        self.cfd_utils.ComputeConvectiveContribution(self.N, grad_u, a, self.Pi)

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

        self.ElemData(p, self.connectivity, self.scalar_elemental_data)
        pel = self.scalar_elemental_data #no copy ... just renaming for better readability

        self.ElemData(vold, self.connectivity, self.vec_elemental_data)
        vel = self.vec_elemental_data #no copy ... just renaming for better readability

        self.tau = self.ComputeTau(self.h,vel,self.nu,dt,self.rho)

        ##TODO: a lot of intermediate storage could be avoided reusing vectors

        # --- k1 ---
        k1 = self.ComputeVelocityResidual(vel,pel,self.DN,self.tau)
        k1.reshape(-1)[:] *= self.Minv

        # --- k2 ---    
        v2 = vold + 0.5 * dt * k1
        self.ApplyVelocityDirichletConditions(vold,self.v_end_of_step_dirichlet,0.5,v2)

        self.ElemData(v2, self.connectivity, vel)
        k2 = self.ComputeVelocityResidual(vel, pel, self.DN,self.tau)
        k2.reshape(-1)[:] *= self.Minv
        
        # --- k3 ---
        v3 = vold + 0.5 * dt * k2
        self.ApplyVelocityDirichletConditions(vold,self.v_end_of_step_dirichlet,0.5,v3)

        self.ElemData(v3, self.connectivity, vel)
        k3 = self.ComputeVelocityResidual(vel, pel, self.DN,self.tau)
        k3.reshape(-1)[:] *= self.Minv

        # --- k4 ---
        v4 = vold + dt * k3
        self.ApplyVelocityDirichletConditions(vold,self.v_end_of_step_dirichlet,1.0,v4)

        self.ElemData(v4, self.connectivity, vel)
        k4 = self.ComputeVelocityResidual(vel, pel, self.DN,self.tau)
        k4.reshape(-1)[:] *= self.Minv
        
        # --- final RK4 update ---
        vnew = vold + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
        self.ApplyVelocityDirichletConditions(vold, self.v_end_of_step_dirichlet,1.0,vnew)

        return vnew

    def SolveStep2(self,vfrac,p,dt):
        nelem = self.DN.shape[0]
        n_in_el = self.DN.shape[1]
        dim = self.DN.shape[2]

        ##allocation of intermediate arrays #TODO: allocate it only once
        Lel = np.empty((nelem,n_in_el,n_in_el))
        N_div = np.empty((nelem,n_in_el))
        RHSel = np.empty((nelem,n_in_el))
        

        #gather data from the database
        self.ElemData(p, self.connectivity, self.scalar_elemental_data)
        pel = self.scalar_elemental_data #no copy ... just renaming for better readability

        self.ElemData(vfrac, self.connectivity, self.vec_elemental_data)
        vel = self.vec_elemental_data #no copy ... just renaming for better readability

        # the RHS L_IJ := (∇N_I,∇N_J)
        self.cfd_utils.ComputeLaplacianMatrix(self.DN,Lel)
        #L *= self.tau[:,np.newaxis,np.newaxis] 
        Lel *= self.elemental_volumes[:,np.newaxis,np.newaxis]


        # (q,∇·u)
        self.cfd_utils.ComputeElementwiseNodalDivergence(self.N, self.DN, vel, N_div)
        RHSel = N_div
        

        # coeff = self.dt ##TODO: it should be 1/(density+BDFCoeff0
        # RHSel += coeff*np.einsum("eij,ej->ei", Lel , pel) #RHS -= LHS@pel #TODO: change to avoid temporaries

        # tmp = np.zeros(RHS_el.shape) #TODO: use an existing temporary to avoid allocation
        # self.cfd_utils.ComputePressureStabilization_Proj(self.DN, pel, tmp)
        # tmp *= self.tau[:,np.newaxis]
        # RHSel += tmp
        
        RHSel *= -self.elemental_volumes[:,np.newaxis] #REMARK: observe the - sign!
        RHSpressure = np.zeros(vfrac.shape[0])
        self.cfd_utils.AssembleVector(self.connectivity,RHSel,RHSpressure)
        RHSpressure *= self.rho



        self.cfd_utils.AssembleScalarMatrix(Lel,self.connectivity,self.L)
        Lvalues = self.L.value_data() #TODO: add to the interface the possibility of doing L*=dt
        Lvalues *= dt 

        # dt*(∇q,∇pold)
        RHSpressure += (self.L@p) #TODO: use SpMV function to avoid temporaries

        for i in range(RHSpressure.shape[0]):
            print(i+1,RHSpressure[i])


        #apply Pressure BCs
        print("fix_prss_indices",self.fix_press_indices)
        self.pressure_free_dofs_vector = np.ones(self.nnodes,np.float64) #TODO: this vector can be computed once if the fixity is fixed
        self.pressure_free_dofs_vector[self.fix_press_indices] = 0.0
        self.L.ApplyHomogeneousDirichlet(self.pressure_free_dofs_vector, 1e9, RHSpressure)

        #solve the system
        ##TODO: use amgcl instead of this one!! ... and also avoid creating temporaries
        Ascipy = KM.scipy_conversion_tools.to_csr(self.L)
        p = scipy.sparse.linalg.spsolve(Ascipy, RHSpressure)

        return p
      

    def SolveStep3(self,v,vfrac,delta_p,dt):
        
        grad_dp_el = np.empty(self.DN.shape) #TODO: allocate it only once(or even better reuse one that is available)
    
        self.ElemData(delta_p, self.connectivity, self.scalar_elemental_data)
        delta_pel = self.scalar_elemental_data #no copy ... just renaming for better readability

        self.cfd_utils.Compute_DN_N(self.N, self.DN, delta_pel, grad_dp_el)
        grad_dp_el *= self.elemental_volumes[:,np.newaxis,np.newaxis]
        
        grad_dp_nodal = np.zeros((self.nnodes, self.dim)) #TODO: allocate it only once
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
        #vfrac = v.copy()
        p = self.SolveStep2(vfrac,pold,dt)
        delta_p = p - pold #TODO: most probably we could modify p in place
        v = self.SolveStep3(v,vfrac,delta_p,dt) 
        #v = vfrac.copy() #TODO: remove! this is just mocking step2 and step3
        
        #update Kratos
        self.v_adaptor.data = v
        self.v_adaptor.StoreData()
        self.p_adaptor.data = p
        self.p_adaptor.StoreData()
        value = input("Press Enter to continue...")
        

    def ComputeVelocityResidual(self,v_elemental,p_elemental,DN,tau):
        nelem = v_elemental.shape[0]
        N = self.N

        # assign a shorter local name
        grad_v   = self.grad_v
        v_gauss  = self.v_gauss
        aux_res_el = self.aux_res_el
        res        = self.res

        #(w,b)
        res.fill(0.0) #TODO: substitute by the proper computation of the external forces

        #(w,a·∇u)
        convective=aux_res_el #rename this to make it more readable
        self.cfd_utils.ComputeElementalGradient(DN, v_elemental, grad_v)
        self.cfd_utils.InterpolateValue(self.N,v_elemental,v_gauss)
        self.cfd_utils.ComputeConvectiveContribution(N, grad_v, v_gauss, convective)
        convective *= self.rho 
        res -= convective

        #(∇w,∇u)
        viscous=aux_res_el #rename this to make it more readable
        self.cfd_utils.ApplyLaplacian(DN,v_elemental,viscous)
        viscous *= self.nu #TODO, check if we use nu or mu
        res -= viscous

        # (∇w,p_gauss)
        pressure_term=aux_res_el #rename this to make it more readable
        self.cfd_utils.Compute_DN_N(N,DN,p_elemental,pressure_term)
        res += pressure_term

        # (a·∇u,a·∇u) - (a·∇u,Pi_conv)
        convective_stabilization=aux_res_el #rename this to make it more readable
        Pi_elemental = np.zeros(v_elemental.shape) #TODO: this has to be computed in the appropriate way!
        self.cfd_utils.ComputeMomentumStabilization(N, DN, v_gauss, v_elemental, Pi_elemental, convective_stabilization)
        convective_stabilization *= self.tau[:,np.newaxis,np.newaxis] #TODO: activate stabilization
        res -= convective_stabilization

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