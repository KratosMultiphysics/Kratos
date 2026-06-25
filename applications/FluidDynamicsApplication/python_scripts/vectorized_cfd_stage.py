import math
import time
import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.analysis_stage as analysis_stage
import KratosMultiphysics.FluidDynamicsApplication.python_pool as python_pool
import KratosMultiphysics.FluidDynamicsApplication.cfd_utils as cfd_utils
import KratosMultiphysics.scipy_conversion_tools as scipy_conversion_tools
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
from cupyx.scipy.sparse.linalg import LinearOperator

xp = None

class JacobiPreconditioner(LinearOperator):
    def __init__(self, A):
        """
        Computes and stores the inverse diagonal of A just once upon initialization.
        """
        self.shape = A.shape
        self.dtype = A.dtype

        # Store the inverse diagonal as a class attribute
        self.inv_diag = 1.0 / A.diagonal()

    def _matvec(self, v):
        """
        Applies the preconditioner to vector v.
        CuPy's solver automatically calls this method.
        """
        return self.inv_diag * v

class VectorizedCFDStage(analysis_stage.AnalysisStage):
    def __init__(self, model, project_parameters):

        # Get configuration from problem data settings
        problem_data = project_parameters["problem_data"]
        precision = problem_data["precision"].GetString() if problem_data.Has("precision") else "float64"
        parallel_type = problem_data["parallel_type"].GetString() if problem_data.Has("parallel_type") else "open_mp"

        # Set the backend environment
        # Note that this imports the modules corresponding to the parallelism (i.e., CuPy, SciPy, ...)
        cfd_utils.configure(parallel_type, precision)
        global xp, USE_CUPY
        xp = cfd_utils.xp
        USE_CUPY = cfd_utils.USE_CUPY

        # Create CFDUtils instance
        self.cfd_utils = cfd_utils.CFDUtils()

        # Get and validate the solving settings
        # Note that these are encapsulated within the customary "solver_settings" block
        settings = project_parameters["solver_settings"]
        settings.ValidateAndAssignDefaults(self._GetDefaultSolvingSettings())

        # Either retrieve the model part from the model or create a new one
        model_part_name = settings["model_part_name"].GetString()
        if model_part_name == "":
            raise ValueError('Please provide the model part name as the "model_part_name" (string) parameter!')

        if model.HasModelPart(model_part_name):
            self.model_part = model.GetModelPart(model_part_name)
        else:
            self.model_part = model.CreateModelPart(model_part_name)

        # Set the problem dimension
        self.dim = settings["domain_size"].GetInt()
        if self.dim not in (2,3):
            raise ValueError(f"Provided domain size is '{self.dim}'. Only 2 or 3 are supported.")
        self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, self.dim)

        # Set the number of nodes in elements and conditions
        # Note that only linear simplicial elements are supported (i.e., linear triangle and tetrahedron)
        self.n_in_el = 3 if self.dim == 2 else 4
        self.n_in_cond = 2 if self.dim == 2 else 3

        # Set time integration parameters
        self.dt = settings["time_stepping"]["time_step"].GetDouble()
        self.max_cfl = settings["time_stepping"]["max_cfl_number"].GetDouble()
        self.cfl = self.max_cfl
        self.max_fourier = settings["time_stepping"]["max_fourier_number"].GetDouble()

        self.time_scheme = settings["time_scheme"].GetString()
        self.startup_time = settings["startup_time"].GetDouble()
        self.startup_time_scheme = settings["startup_time_scheme"].GetString() if self.startup_time > 0.0 else None

        # Set problem solving parameters
        self.clear_divergence_steps = settings["clear_divergence_steps"].GetInt()
        self.pressure_max_iteration = settings["linear_solver_settings"]["max_iteration"].GetInt()
        self.pressure_tolerance = settings["linear_solver_settings"]["tolerance"].GetDouble()

        # Save materials import settings
        self.material_import_settings = settings["material_import_settings"]

        #flags to control algorithmic behaviour
        self.deactivate_pressure_stabilization = False #switches on and off the pressure stabilization
        self.force_first_order_splitting = False

        # Call base analysis stage constructor
        # Note that this must be done at the end (indeed, after creating the model part)
        # Otherwise the model part is not created when calling the _AddVariables()
        super().__init__(model,project_parameters)

    def ComputeLumpedMass(self):
        Mscalar = xp.zeros(len(self.model_part.Nodes), dtype=cfd_utils.PRECISION)
        Mel = xp.einsum("e,i->ei", self.elemental_volumes, self.N)
        self.cfd_utils.AssembleVector(self.connectivity, Mel, Mscalar)
        M = xp.tile(Mscalar[:, xp.newaxis], (1, self.dim)) # Broadcast the scalar mass to all velocity components
        return M.ravel()

    def ApplyVelocityDirichletConditions(self, vold, v_end_of_step, dt_factor, v):
        v.ravel()[self.fix_vel_indices] = (1-dt_factor)*vold.ravel()[self.fix_vel_indices] + dt_factor*v_end_of_step

    def ApplyVelocitySlipConditions(self, v, normals): #TODO: normalize normals only once
        # Get velocity and normal view of slip nodes
        n_slip_nodes = len(self.slip_node_ids)
        v_slip = v[self.slip_node_ids]
        n_unit = normals[self.slip_node_ids]

        # Check for zero normals and print corresponding indices
        normal_norms = xp.linalg.norm(n_unit, axis=1)
        if self.zero_normal_indices is None:
            zero_normal_mask = (normal_norms < 1e-15)
            self.zero_normal_indices = self.slip_node_ids[xp.where(zero_normal_mask)[0]]

        # Normalize normals (avoid division by zero)
        normal_norms_clipped = xp.maximum(normal_norms, 1e-15)
        n_unit /= normal_norms_clipped[:,None]

        # Remove normal velocity component
        v_dot_n = xp.sum(v_slip * n_unit, axis=1)
        v[self.slip_node_ids] -= (v_dot_n[:,None] * n_unit)

        # Set velocity to zero at nodes with zero normal (i.e., treat them as no-slip)
        v[self.zero_normal_indices,:] = 0.0

    def ComputeH(self,DN):
        return xp.mean(1.0 / xp.linalg.norm(DN, axis=2), axis=1)

    def ComputeTau(self, h, rho_conv , viscosity, rho, dt):
        """
        Computes the VMS stabilization constants

        Parameters
        ----------
        h : average element size used in the diffusion characteristic time calculation.

        rho_conv : convective operator spectral radius used in the advection characteristic time calculation.

        viscosity : fluid dynamic viscosity

        rho : fluid density

        dt : time step
        """

        c_1 = 4.0
        c_2 = 2.0

        tau_1 = self.tau_1
        tau_2 = self.tau_2

        denom_max = 1.0e-3*rho/dt
        denom = xp.maximum(denom_max, c_2*rho*rho_conv + c_1*viscosity/h**2) #TODO: check if we can reuse Fourier and CFL numbers here
        xp.reciprocal(denom, out=tau_1)
        tau_2[:] = h**2/(c_1 * tau_1)

        # if self.clear_divergence_steps > 0:
        #     self.deactivate_pressure_stabilization = True
        #     self.force_first_order_splitting = True
            #tau_1.fill(0.0)
            #tau_2.fill(0.0)

        return tau_1, tau_2

    def ComputeElementalConvectiveOperatorSpectralRadius(self, v, out=None):
        # Get the average velocity at each element
        v_elemental = self.ElemData(v, self.connectivity, out=self.v_el) # shape (num_elements, num_nodes_per_element, dim)
        v_el_mean = xp.sum(v_elemental, axis=1) / self.n_in_el #TODO: avoid allocation

        nelem, nnode = self.DN.shape[0], self.DN.shape[1]
        if out is None:
            out = xp.empty(nelem, dtype=v_el_mean.dtype)

        advective_proj = self.pool.Get(0,(nelem, nnode)) #self.tmp_n_in_el

        # Calculate the advective projection v · ∇N_i for each node in each element representing the discrete convective operator
        # For linear (P1) elements, shape function gradients are constant inside the element and scale like 1/h (i.e., |∇N_i| ~ 1/h_i)
        # In consequence, this term scales like |v|/h_i
        # Replace einsum('ek, enk -> en') with matmul: (nelem, 1, dim) @ (nelem, dim, nnode) -> (nelem, 1, nnode)
        xp.matmul(self.DN, v_el_mean[:, :, xp.newaxis], out=advective_proj[:, :, xp.newaxis])

        # # Get the convective spectral radius (inverse time scale) as the average of the values in each direction
        # # Note that this is an approximation of the elemental convective operator spectral radius
        # return xp.sum(xp.abs(out), axis=1) / self.n_in_el

        # # Get the convective operator row-sum as an upper bound of its spectral radius (inverse time scale)
        # # Note that using this as an approximation of the spectral radius results in a conservative estimation of the time increment when computing the CFL
        # # On the contrary, when using it in the subscale stabilization factor calculation (tau_1) results in a smaller contribution of the subscales
        # return xp.sum(xp.abs(out), axis=1)

        # Get the convective spectral radius (inverse time scale) as the maximum value in each nodal direction
        # This is equivalent to get the maximum eigenvalue of the elemental convective operator
        # In-place absolute value (No allocation)
        xp.abs(advective_proj, out=advective_proj)

        # Max reduction into a pre-allocated result array
        # result_max = xp.empty(out.shape[0])
        xp.max(advective_proj, axis=1, out=out)
        self.pool.Release(advective_proj)
        return out

    def _PrepareModelPart(self):
        super()._PrepareModelPart()

        # Get fluid properties from materials .json file
        materials_filename = self.material_import_settings["materials_filename"].GetString()
        if materials_filename != "":
            material_settings = KM.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KM.ReadMaterialsUtility(material_settings, self.model)
        else:
            raise ValueError('Please provide the name of the .json file containing the materials data as the "materials_filename" (string) parameter within the "material_import_settings" block!')

        # Save dynamic viscosity and density from model part properties
        # Note that we are assuming that they are defined in the first property of the model part
        if not self.GetComputingModelPart().HasProperties(1):
            raise Exception("Fluid material properties must be provided in property with Id 1!")
        else:
            self.rho = self.GetComputingModelPart().GetProperties(1).GetValue(KM.DENSITY)
            self.dyn_visc = self.GetComputingModelPart().GetProperties(1).GetValue(KM.DYNAMIC_VISCOSITY)

    def _AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.REACTION)
        self.model_part.AddNodalSolutionStepVariable(KM.REACTION_WATER_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        self.model_part.AddNodalSolutionStepVariable(KM.BODY_FORCE)
        self.model_part.AddNodalSolutionStepVariable(KM.EXTERNAL_PRESSURE)

        self.model_part.SetBufferSize(2)

    def _AddDofs(self):
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

    def ElemData(self, v, connectivities, out=None):
        if out is None:
            return xp.take(v, connectivities, axis=0, out=out)
        else:
            if(len(v.shape) == 1): ##getting a scalar
                if out.shape != (connectivities.shape[0], connectivities.shape[1]):
                    print("v.shape", v.shape)
                    print("out.shape", out.shape)
                    print("connectivities.shape", connectivities.shape)
                    raise ValueError(f"Output array has incompatible shape {out.shape}, expected {(connectivities.shape[0], connectivities.shape[1])}")
            else:
                if out.shape != (connectivities.shape[0], connectivities.shape[1], self.dim):
                    print("v.shape", v.shape)
                    print("out.shape", out.shape)
                    print("connectivities.shape", connectivities.shape)
                    raise ValueError(f"Output array has incompatible shape {out.shape}, expected {(connectivities.shape[0], connectivities.shape[1], self.dim)}")
            xp.take(v, connectivities, axis=0, out=out)
            return out

    def _AdvanceTime(self):
        # Get time step and advance in time
        self.time += self.dt # Note that this is the user defined time
        self.model_part.CloneTimeStep(self.time)
        self.model_part.ProcessInfo[KM.STEP] += 1

        return self.time

    def AllocateWorkArrays(self):
        ##arrays to store elemental values
        self.v_el = xp.empty(self.DN.shape, dtype=cfd_utils.PRECISION)
        self.p_el = xp.empty((self.DN.shape[0], self.n_in_el), dtype=cfd_utils.PRECISION)
        self.b_el = xp.empty(self.DN.shape, dtype=cfd_utils.PRECISION)
        self.div_proj_el = xp.empty(self.p_el.shape, dtype=cfd_utils.PRECISION)
        self.conv_proj_el = xp.empty(self.DN.shape, dtype=cfd_utils.PRECISION)
        self.pi_conv = xp.empty((self.nnodes,self.dim), dtype=cfd_utils.PRECISION)
        self.tau_1 = xp.empty(self.DN.shape[0], dtype=cfd_utils.PRECISION)
        self.tau_2 = xp.empty(self.DN.shape[0], dtype=cfd_utils.PRECISION)

        ## arrays to store wall values
        self.v_el_wall = xp.empty((self.connectivity_wall.shape[0], self.n_in_cond, self.dim), dtype=cfd_utils.PRECISION)

        # Arrays for time integration to be allocated according to the current time scheme
        # Note that their size may depend on the time scheme, which might change among time steps
        self.k = None
        self.v_stage = xp.empty((len(self.model_part.Nodes), self.dim), dtype=cfd_utils.PRECISION)

        #work arrays for intermediate computations
        self.pool = python_pool.BufferPool(
            sizes=[self.DN.shape[0]*self.n_in_el*self.n_in_el, #0 -note that the first one is larger
                   self.DN.shape[0]*self.n_in_el*self.dim, #1
                   self.DN.shape[0]*self.n_in_el*self.dim, #2
            ],
            xp=xp,
            dtype=cfd_utils.PRECISION
        )


        if USE_CUPY:
            ##CLEAN FREE GPU MEMORY POOL just in case

            # 1. Clear the main device memory pool
            xp.get_default_memory_pool().free_all_blocks()

            # 2. Clear the pinned memory pool (CPU memory used for fast transfers)
            xp.get_default_pinned_memory_pool().free_all_blocks()


    def Initialize(self):
        super().Initialize()

        self.model_part.ProcessInfo[KM.STEP] = 0

        vol_geometries = []
        wall_geometries = []
        outlet_geometries = []
        for g in self.model_part.Geometries:
            if (len(g) == self.dim + 1):
                vol_geometries.append(g.Id)
            elif (len(g) == self.dim):
                wall_nodes = 0
                outlet_nodes = 0
                for node in g:
                    if node.Is(KM.WALL):
                        wall_nodes += 1
                    elif node.Is(KM.OUTLET):
                        outlet_nodes += 1
                if wall_nodes == len(g):
                    wall_geometries.append(g.Id)
                elif outlet_nodes == len(g):
                    outlet_geometries.append(g.Id)
            else:
                raise Exception (f"Only linear simplicial elements are supported! Geometry with Id {g.Id} has {len(g)} nodes, expected {self.dim} (skin) or {self.dim+1} (volume).")

        self.vol_mp = self.model_part.CreateSubModelPart("fluid_computational_model_part") #TODO: I think we don't really need this. We can simply take the one defined in the json.
        self.vol_mp.AddGeometries(vol_geometries)

        self.wall_mp = self.model_part.CreateSubModelPart("fluid_wall_model_part")
        self.wall_mp.AddGeometries(wall_geometries)

        self.outlet_mp = self.model_part.CreateSubModelPart("fluid_outlet_model_part")
        self.outlet_mp.AddGeometries(outlet_geometries)

        self.nnodes = len(self.model_part.Nodes)

        self.v_adaptor_n = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.VELOCITY,data_shape=[self.dim],step_index=1)
        self.v_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.VELOCITY,data_shape=[self.dim],step_index=0)
        self.p_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.PRESSURE,0)
        self.p_adaptor_n = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, KM.PRESSURE, 1)
        self.b_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, KM.BODY_FORCE, data_shape=[self.dim], step_index=0)
        self.b_adaptor_n = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, KM.BODY_FORCE, data_shape=[self.dim], step_index=1)
        self.normals_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, KM.NORMAL, data_shape=[self.dim], step_index=0)

        self.v_adaptor.CollectData()
        v = self.v_adaptor.data.reshape((len(self.model_part.Nodes),self.dim))
        v = xp.asarray(v, dtype=cfd_utils.PRECISION)

        self.p_adaptor.CollectData()
        p = self.p_adaptor.data
        p = xp.asarray(p, dtype=cfd_utils.PRECISION)

        self.connectivity_adaptor = KM.TensorAdaptors.ConnectivityIdsTensorAdaptor(self.vol_mp.Geometries)
        self.connectivity_adaptor.CollectData()
        self.connectivity = self.connectivity_adaptor.data
        self.connectivity -= 1 #have indices to start in 0
        self.connectivity = self.connectivity.astype(np.int64) #############OUCH!!! TODO: this by defaults should be probably int64

        self.connectivity_adaptor_wall = KM.TensorAdaptors.ConnectivityIdsTensorAdaptor(self.wall_mp.Geometries)
        self.connectivity_adaptor_wall.CollectData()
        self.connectivity_wall = self.connectivity_adaptor_wall.data
        self.connectivity_wall -= 1 #have indices to start in 0
        self.connectivity_wall = self.connectivity_wall.astype(np.int64) #############OUCH!!! TODO: this by defaults should be probably int64

        # Preallocation of elemental arrays of velocities, pressures and body force
        # self.vec_elemental_data = np.empty((*self.connectivity.shape, v.shape[1]))
        # self.vec_elemental_data_b = np.empty((*self.connectivity.shape, v.shape[1]))
        # self.scalar_elemental_data = np.empty(self.connectivity.shape)

        # Preallocation of local array of shape functions
        self.N = np.ones(self.dim+1)/(self.dim+1)
        self.N_wall = np.ones(self.dim)/(self.dim)

        # Obtain the shape function derivatives
        geometry_adaptor_DN = KM.TensorAdaptors.GeometriesTensorAdaptor(
            self.vol_mp.Geometries,
            KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.ShapeFunctionDerivatives,
            KM.GeometryData.IntegrationMethod.GI_GAUSS_1)
        geometry_adaptor_DN.CollectData()
        self.DN = xp.squeeze(geometry_adaptor_DN.data).copy() #this has shape nel*1*nnodes_in_el*dim - the copy is important as we need to own the data

        # Obtain elemental volumes #TODO: an adaptor should be available for these
        vols = []
        for geom in self.vol_mp.Geometries:
            vols.append(geom.DomainSize())
        self.elemental_volumes = xp.array(vols)

        wall_areas = []
        for geom in self.wall_mp.Geometries:
            wall_areas.append(geom.DomainSize())
        self.wall_areas = xp.array(wall_areas)

        # Compute lumped mass matrix and its inverse
        self.M=self.ComputeLumpedMass() #TODO: decide if we need to store M
        self.Minv = 1./self.M

        # Allocate the graph of the Laplacian matrix
        self.L, self.L_assembly_indices = self.cfd_utils.AllocateScalarMatrix(self.connectivity)

        # Move stuff to the GPU
        self.DN = xp.asarray(self.DN, dtype=cfd_utils.PRECISION)
        self.N = xp.asarray(self.N, dtype=cfd_utils.PRECISION)
        self.N_wall = xp.asarray(self.N_wall, dtype=cfd_utils.PRECISION)
        max_idx = self.connectivity.max()
        if max_idx >= np.iinfo(np.int32).max: #use 32 bit integers for connectivities
            print(f"Warning: Max index {max_idx} is too large for uint32")

        self.connectivity = xp.asarray(self.connectivity, dtype=xp.int32)
        self.connectivity_wall = xp.asarray(self.connectivity_wall, dtype=xp.int32)

        # Compute h per every element
        self.h = self.ComputeH(self.DN)

        # Get nodal normals for the slip condition
        # Note that the normals are already computed by the corresponding process Initialize() method (see base class)
        self.normals_adaptor.CollectData()
        self.normals = self.normals_adaptor.data.reshape((len(self.model_part.Nodes),self.dim))
        self.normals = xp.asarray(self.normals, dtype=cfd_utils.PRECISION)

        # Get wall normals for the wall condition
        # Note that constant normals are assumed within the wall geometries by using a single integration point quadrature (no curved boundaries are considered)
        self.normals_wall_adaptor = KM.TensorAdaptors.GeometriesTensorAdaptor(
            self.wall_mp.Geometries,
            KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.Normals,
            KM.GeometryData.IntegrationMethod.GI_GAUSS_1)
        self.normals_wall_adaptor.CollectData()
        n_wall_geom = self.normals_wall_adaptor.data.shape[0]
        n_wall_int_pts = self.normals_wall_adaptor.data.shape[1]
        print(self.normals_wall_adaptor.data.data.shape)
        self.normals_wall = self.normals_wall_adaptor.data.reshape((n_wall_geom, n_wall_int_pts, 3))
        self.normals_wall = xp.asarray(self.normals_wall[:,:,:self.dim], dtype=cfd_utils.PRECISION) # Drop third component in 2D
        if(self.normals_wall_adaptor.data.shape[1]!=0):
            self.normals_wall = xp.squeeze(self.normals_wall, axis=1) # Drop integration point dimension (note that we are using a single integration point, so this is safe)

        normals_wall_norms = xp.linalg.norm(self.normals_wall, axis=1, keepdims=True)
        self.unit_normals_wall = self.normals_wall / (normals_wall_norms + 1.0e-12) # avoid division by zero

        # Calculate the wall tangential projector (I - n⊗n) to be used in the wall condition
        #I = xp.eye(self.dim, dtype=cfd_utils.PRECISION)
        #self.wall_tangential_projector = I[None, :, :] - (self.unit_normals_wall[:, :, None] * self.unit_normals_wall[:, None, :])

        # Declare boundary condition arrays
        self.fix_vel_indices = None
        #self.slip_vel_indices = None
        self.slip_node_ids = None
        self.zero_normal_indices = None
        self.fix_pres_indices = None

        # Declare boundary conditions masks
        self.csr_data_indices = None
        self.csr_diag_indices = None

        # Compute the Fourier number time step restriction (constant as h, rho and mu are constant in time)
        aux = self.max_fourier * self.rho / self.dyn_visc
        self.dt_fourier = np.min(aux * self.h**2)

        # Assemble Laplacian matrix (w/o stabilization) to get the graph for AMG
        L_el = self.cfd_utils.ComputeLaplacianMatrix(self.DN) # elemental laplacian contributions as L_IJ := (∇N_I,∇N_J)
        L_el *= self.elemental_volumes[:, None, None] # scale Laplaciant elemental contributions
        self.cfd_utils.AssembleScalarMatrixByCSRIndices(L_el, self.L_assembly_indices, self.L) # assemble the scaled elemental contributions
        #print(f"Setting graph for AMG with {self.L.nnz} nonzeros and {self.L.shape[0]} rows")
        t0 = time.perf_counter()
        self.preconditioner = self.cfd_utils.ConstructPreconditioner(self.L)
        print(f"AMG graph setup time: {time.perf_counter()-t0}")

        #FIXME: remove after developing
        self.step_1_total_time = 0.0
        self.step_2_total_time = 0.0
        self.step_3_total_time = 0.0
        self.init_time = time.perf_counter()

        self.AllocateWorkArrays()

        #self.outfile = open("refres.res", "w") #please do not remove this, it is used for testing purposes



    def Finalize(self):
        super().Finalize()

        self._FinalizePressureLinearSolver()

        total_time = time.perf_counter()-self.init_time
        print(f"\n\nSimulation total time: {total_time}")
        print(f"Step 1 total time: {self.step_1_total_time} ({round(100.0 * self.step_1_total_time / total_time)}%)")
        print(f"Step 2 total time: {self.step_2_total_time} ({round(100.0 * self.step_2_total_time / total_time)}%)")
        print(f"Step 3 total time: {self.step_3_total_time} ({round(100.0 * self.step_3_total_time / total_time)}%)")

    def GetComputingModelPart(self):
        """This function provides a unified way to access the computing model part from outside the stage.
        It can be overriden in derived classes. By default, it returns the solver computing model part.
        """
        return self.model_part

    def ComputeVelocityProjection(self, v_elemental, out=None):
        if out is None:
            raise ValueError("Output array must be provided to store the pressure projection values!")
        if out.shape != (self.nnodes,self.dim):
            raise ValueError(f"Output array has incompatible shape {out.shape}, expected {(self.nnodes,)}")
        pi_conv = out


        # Solve (w,pi) = (w,a·∇u)
        # Compute convective contribution at integration points (w,a·∇u)

        grad_v_elemental = self.cfd_utils.ComputeElementalGradient(self.DN, v_elemental, out=self.pool.Get(0,(self.DN.shape[0], self.dim, self.dim))) # Calculate the velocity gradient at each element (assumed constant within the element)
        a_grad_elemental = self.cfd_utils.ComputeElementalConvectiveOperator(v_elemental, grad_v_elemental, out=self.pool.Get(1,self.DN.shape)) # Calculate the convective operator at each node
        self.pool.Release(grad_v_elemental)

        convective = self.cfd_utils.ComputeConvectiveContribution(a_grad_elemental, out=self.pool.Get(2,self.DN.shape)) # Calculate the elemental convective contributions (note that density should be applied outside)
        self.pool.Release(a_grad_elemental)

        convective *= self.elemental_volumes[:, None, None] # Apply Jacobian determinant to the entire residual

        # Do the nodal assembly of the projection elemental contributions
        pi_conv.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity, convective, pi_conv)
        self.pool.Release(convective)

        # Solve with the lumped mass matrix to get the nodal values of the projection
        pi_conv = pi_conv.ravel()
        pi_conv *= self.Minv

        # Return the projection as a nodal data array (nnodes*dim)
        return pi_conv.reshape((self.nnodes, self.dim))

    def ComputeDivergenceProjection(self, v_elemental, out=None):
        if out is None:
            raise ValueError("Output array must be provided to store the pressure projection values!")
        if out.shape != (self.nnodes,):
            raise ValueError(f"Output array has incompatible shape {out.shape}, expected {(self.nnodes,)}")
        pi_div = out

        # Solve (q,pi) = (q,∇·v)
        # Calculate the divergence term at the elemental level
        aux_scalar = self.cfd_utils.ComputeElementwiseNodalDivergence(self.N, self.DN, v_elemental, out=self.pool.Get(1,self.p_el.shape)) # shape (num_elements, num_nodes_per_element)
        aux_scalar *= self.elemental_volumes[:,xp.newaxis]

        # Do the nodal assembly of the projection elemental contributions
        pi_div.fill(0.0) # = xp.zeros((self.nnodes), dtype=cfd_utils.PRECISION) #TODO: can we allocate this once and reuse it?
        self.cfd_utils.AssembleVector(self.connectivity, aux_scalar, pi_div)
        self.pool.Release(aux_scalar)

        # Solve with the lumped mass matrix to get the nodal values of the projection
        pi_div *= self.Minv[::self.dim] # Strided view to take one value every self.dim
        return pi_div

    def ComputePressureProjection(self, p_elemental, out=None):
        if out is None:
            raise ValueError("Output array must be provided to store the pressure projection values!")
        if out.shape != (self.nnodes,self.dim):
            raise ValueError(f"Output array has incompatible shape {pi_press.shape}, expected {(self.nnodes,self.dim)}")
        pi_press = out

        # Solve (w,pi) = (w,∇p)
        # Calculate the elemental pressure gradient projection contributions
        pi_press_el = self.cfd_utils.ComputeNDN(self.N, self.DN, p_elemental, out=self.pool.Get(0,self.DN.shape)) # shape (num_elements, num_nodes_per_element, dim)
        pi_press_el *= self.elemental_volumes[:,xp.newaxis,xp.newaxis]

        # Do the nodal assembly of the projection elemental contributions
        pi_press.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity, pi_press_el, pi_press)
        self.pool.Release(pi_press_el)

        # Solve with the lumped mass matrix to get the nodal values of the projection
        pi_press = pi_press.ravel()
        pi_press *= self.Minv

        # Return the projection values as a nodal data array (nnodes*dim)
        return pi_press.reshape((self.nnodes, self.dim))

    def GetFixityIndices(self):
        # Obtain fixity
        # Note that we only compute this once assuming the fixity to not change in time
        if self.fix_vel_indices is None or self.fix_pres_indices is None:
            # Get fixity nodal data
            dofs_list = [KM.VELOCITY_X, KM.VELOCITY_Y, KM.VELOCITY_Z, KM.PRESSURE] if self.dim == 3 else [KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE]
            fixity_adaptor = KM.TensorAdaptors.FixityTensorAdaptor(self.model_part.Nodes, dofs_list)
            fixity_adaptor.CollectData()

            # Create a flat array with the DOFs indices from the fixity data
            n_nodes = fixity_adaptor.data.shape[0]
            p_eq_indices = xp.arange(n_nodes)
            v_eq_indices = xp.arange(n_nodes)[:, xp.newaxis] * self.dim + xp.arange(self.dim)

            # Mask where the pressure and velocity is fixed
            p_fixity_mask = fixity_adaptor.data[:,-1].ravel() == True
            v_fixity_mask = fixity_adaptor.data[:,:-1].ravel() == True

            # Apply the fixity masks to the DOF equation indices
            self.fix_pres_indices = xp.asarray(p_eq_indices[p_fixity_mask], dtype=np.int32)
            self.fix_vel_indices = xp.asarray(v_eq_indices.ravel()[v_fixity_mask], dtype=np.int32)

    def GetSlipIndices(self):
        # Check if slip indices have been computed
        # Note that we only compute this once assuming the slip to not change in time
        if self.slip_node_ids is None:
            # Get SLIP flag nodal data
            slip_adaptor = KM.TensorAdaptors.FlagsTensorAdaptor(self.model_part.Nodes, KM.SLIP)
            slip_adaptor.CollectData()
            slip_adaptor_data = xp.asarray(slip_adaptor.data, dtype=np.int32)

            # Compute slip nodes component indices
            self.slip_node_ids = xp.where(slip_adaptor_data > 0)[0] # Get the ids of the nodes with SLIP flag
            #self.slip_vel_indices = (slip_node_ids[:, None] * self.dim + xp.arange(self.dim)).ravel() # Broadcast with dimension to get the component indices

    @classmethod
    def GetButcherTableau(self, time_scheme):
        if time_scheme == "FE": # Forward-Euler
            A = (
                (),
            )
            b = (1.0,)
            c = (0.0,)
        elif time_scheme == "RK2": # 2nd order Heun Runge-Kutta
            A = (
                (),
                (1.0),
            )
            b = (0.5, 0.5)
            c = (0.0, 1.0)
        elif time_scheme == "SSPRK3": # 3rd order Strong Stability Preserving Runge-Kutta
            A = (
                (),
                (1.0,),
                (0.25, 0.25),
            )
            b = (1.0/6.0, 1.0/6.0, 2.0/3.0)
            c = (0.0, 1.0, 0.5)
        elif time_scheme == "RK4": # 4th order Runge-Kutta
            A = (
                (),
                (0.5,),
                (0.0, 0.5),
                (0.0, 0.0, 1.0),
            )
            b = (1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0)
            c = (0.0, 0.5, 0.5, 1.0)
        else:
            available_time_schemes = ["FE", "RK2", "SSPRK3", "RK4"]
            raise ValueError(f"Provided time scheme '{time_scheme}' is not available. Available options are: {available_time_schemes}.")

        return A, b, c

    @classmethod
    def GetPressureFactor(self, time_scheme):
        """
        Pressure factor coming from the intermediate RK evaluation of the pressure term
        """
        if time_scheme == "FE": # Forward-Euler
            return 1.0
        elif time_scheme == "RK2": # 2nd order Heun Runge-Kutta
            return 0.5
        elif time_scheme == "SSPRK3": # 3rd order Strong Stability Preserving Runge-Kutta
            return 0.5
        elif time_scheme == "RK4": # 4th order Runge-Kutta
            return 0.5
        else:
            available_time_schemes = ["FE", "RK2", "SSPRK3", "RK4"]
            raise ValueError(f"Provided time scheme '{time_scheme}' is not available. Available options are: {available_time_schemes}.")

    def SolveStep1(self,vold,v_dirichlet,p,b,dt):

        #xp.cuda.Stream.null.synchronize() #FIXME: remove after debugging, just to be sure we are not measuring asynchronously some previous step computations

        # Gather elemental data from the database
        pel = self.ElemData(p, self.connectivity, out=self.p_el)
        vel = self.ElemData(vold, self.connectivity, out=self.v_el)
        b_el = self.ElemData(b, self.connectivity, out=self.b_el)
        #v_el_wall = self.ElemData(vold, self.connectivity_wall, out=self.v_el_wall)

        # Compute convective projection
        t0 = time.perf_counter()
        conv_proj = self.ComputeVelocityProjection(vel, self.pi_conv)
        conv_proj_el = self.ElemData(conv_proj, self.connectivity, out=self.conv_proj_el)
        #print(f"ComputeVelocityProjection time: {time.perf_counter() - t0}")

        # Compute divergence projection
        t0 = time.perf_counter()
        div_proj = self.ComputeDivergenceProjection(vel, out=self.pool.Get(0,(self.nnodes,))) # compute the divergence projection at the nodes, reusing pool array for output
        div_proj_el = self.ElemData(div_proj, self.connectivity, out=self.div_proj_el)
        self.pool.Release(div_proj)
        #print(f"ComputeDivergenceProjection time: {time.perf_counter() - t0}")

        # Get current time scheme data
        time_scheme = self.time_scheme if self.time > self.startup_time else self.startup_time_scheme
        A, b, c = self.GetButcherTableau(time_scheme)
        n_stages = len(b)
        if (self.k == None or len(self.k) != n_stages):
            self.k = [xp.empty((self.nnodes, self.dim), dtype=cfd_utils.PRECISION) for _ in range(n_stages)]

        for i in range(n_stages):
            # Calculate current stage velocity
            self.v_stage[:] = vold
            for j, a_ij in enumerate(A[i]):
                if a_ij != 0.0:
                    aux = dt * a_ij
                    self.v_stage += aux * self.k[j]
            self.ApplyVelocitySlipConditions(self.v_stage, self.normals)
            self.ApplyVelocityDirichletConditions(vold, v_dirichlet, c[i], self.v_stage)

            # Calculate current stage residual
            v_stage_el = self.ElemData(self.v_stage, self.connectivity, self.v_el)
            #v_stage_el_wall = self.ElemData(self.v_stage, self.connectivity_wall, out=self.v_el_wall)
            self.ComputeVelocityResidual(v_stage_el, pel, b_el, conv_proj_el, div_proj_el, self.DN, self.tau_1, out=self.k[i])
            #self.ComputeVelocityResidualWall(v_stage_el_wall, out=self.k[i])
            self.k[i].reshape(-1)[:] *= (1.0/self.rho) * self.Minv

        # RK final update
        vnew = xp.empty_like(vold)
        vnew[:] = vold
        for i in range(n_stages):
            aux = dt * b[i]
            vnew += aux * self.k[i]
        self.ApplyVelocitySlipConditions(vnew, self.normals)
        self.ApplyVelocityDirichletConditions(vold, v_dirichlet, 1.0, vnew)

        return vnew

        # # Advance in time by runge kutta
        # # --- k1 ---
        # t0 = time.perf_counter()
        # k1 = xp.empty((self.nnodes,self.dim), dtype=cfd_utils.PRECISION)
        # self.ComputeVelocityResidual(vel, pel, b_el, conv_proj_el, div_proj_el, self.DN, self.tau_1, out=k1)
        # #self.ComputeVelocityResidualWall(v_el_wall, out=k1)
        # k1.reshape(-1)[:] *= (1.0/self.rho) * self.Minv
        # #print(f"k1 time: {time.perf_counter() - t0}")

        # # --- k2 ---
        # t0 = time.perf_counter()
        # v2 = vold + 0.5 * dt * k1
        # self.ApplyVelocitySlipConditions(v2, self.normals)
        # self.ApplyVelocityDirichletConditions(vold,v_dirichlet,0.5,v2)

        # t_el_data = time.perf_counter()
        # vel = self.ElemData(v2, self.connectivity, self.v_el)
        # #v_el_wall = self.ElemData(v2, self.connectivity_wall, out=self.v_el_wall) # Get the velocity at the wall nodes for the intermediate velocity v2, reusing pool array for output
        # t_k2_res = time.perf_counter()
        # k2 = xp.empty((self.nnodes,self.dim), dtype=cfd_utils.PRECISION)
        # self.ComputeVelocityResidual(vel, pel, b_el, conv_proj_el, div_proj_el, self.DN,self.tau_1,out=k2)
        # #self.ComputeVelocityResidualWall(v_el_wall, out=k2)
        # #print(f"\tk2 residual: {time.perf_counter() - t_k2_res}")


        # t_k2_update = time.perf_counter()
        # k2.reshape(-1)[:] *= (1.0/self.rho) * self.Minv
        # #print(f"\tk2 update: {time.perf_counter() - t_k2_update}")
        # #print(f"k2 time: {time.perf_counter() - t0}")

        # # --- k3 ---
        # t0 = time.perf_counter()
        # v3 = vold + 0.5 * dt * k2
        # self.ApplyVelocitySlipConditions(v3, self.normals)
        # self.ApplyVelocityDirichletConditions(vold,v_dirichlet,0.5,v3)

        # vel = self.ElemData(v3, self.connectivity, self.v_el)
        # #v_el_wall = self.ElemData(v3, self.connectivity_wall, out=self.v_el_wall) # Get the velocity at the wall nodes for the intermediate velocity v3, reusing pool array for output
        # k3 = xp.empty((self.nnodes,self.dim), dtype=cfd_utils.PRECISION)
        # self.ComputeVelocityResidual(vel, pel, b_el, conv_proj_el, div_proj_el, self.DN,self.tau_1,out=k3)
        # #self.ComputeVelocityResidualWall(v_el_wall, out=k3)
        # k3.reshape(-1)[:] *= (1.0/self.rho) * self.Minv
        # #print(f"k3 time: {time.perf_counter() - t0}")
        # # --- k4 ---
        # t0 = time.perf_counter()
        # v4 = vold + dt * k3
        # self.ApplyVelocitySlipConditions(v4, self.normals)
        # self.ApplyVelocityDirichletConditions(vold,v_dirichlet,1.0,v4)

        # vel = self.ElemData(v4, self.connectivity, self.v_el)
        # #v_el_wall = self.ElemData(v4, self.connectivity_wall, out=self.v_el_wall) # Get the velocity at the wall nodes for the intermediate velocity v4, reusing pool array for output
        # k4 = xp.empty((self.nnodes,self.dim), dtype=cfd_utils.PRECISION)
        # self.ComputeVelocityResidual(vel, pel, b_el, conv_proj_el, div_proj_el, self.DN,self.tau_1,out=k4)
        # #self.ComputeVelocityResidualWall(v_el_wall, out=k4)
        # k4.reshape(-1)[:] *= (1.0/self.rho) * self.Minv
        # #print(f"k4 time: {time.perf_counter() - t0}")

        # # --- final RK4 update ---
        # vnew = vold + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4) #TODO: we can accumulate to save two arrays
        # t0 = time.perf_counter()
        # self.ApplyVelocitySlipConditions(vnew, self.normals)
        # #print(f"Apply SLIP time: {time.perf_counter() - t0}")
        # t0 = time.perf_counter()
        # self.ApplyVelocityDirichletConditions(vold, v_dirichlet,1.0,vnew)
        # #print(f"Apply DIRICHLET time: {time.perf_counter() - t0}")
        # return vnew

    def SolveStep2(self,vfrac,p,dt):

        nelem = self.DN.shape[0]

        # Get factor coming from the RK substep approximation of the pressure term
        time_scheme = self.time_scheme if self.time > self.startup_time else self.startup_time_scheme
        p_factor = self.GetPressureFactor(time_scheme)

        # Gather data from the database
        pel = self.ElemData(p, self.connectivity , self.p_el)
        vel_frac = self.ElemData(vfrac, self.connectivity, self.v_el)

        # Compute pressure projection
        pres_proj = self.ComputePressureProjection(pel, out=self.pool.Get(1,(self.nnodes,self.dim))) # compute the pressure projection at the nodes, reusing pool array for output
        pres_proj_el = self.ElemData(pres_proj, self.connectivity, out=self.conv_proj_el) # reuse conv_proj_el as temporary for the elemental pressure projection
        self.pool.Release(pres_proj)

        # Assemble pressure LHS (set to zero is done internally)
        t0 = time.perf_counter()
        L_el = self.cfd_utils.ComputeLaplacianMatrix(self.DN, out=self.pool.Get(0,(nelem,self.n_in_el, self.n_in_el))) # elemental laplacian contributions as L_IJ := (∇N_I,∇N_J)
        #print(f"Time to compute elemental Laplacian contributions: {time.perf_counter()-t0}")
        if(self.clear_divergence_steps>0 or self.deactivate_pressure_stabilization == True):
            coef = (p_factor * dt / self.rho ) * self.elemental_volumes
        else:
            coef = (p_factor * dt / self.rho + self.tau_1) * self.elemental_volumes
        L_el *= coef[:, None, None] # scale LHS elemental contributions
        t0 = time.perf_counter()
        self.cfd_utils.AssembleScalarMatrixByCSRIndices(L_el, self.L_assembly_indices, self.L) # assemble the scaled elemental contributions
        #print(f"Time to assemble pressure LHS: {time.perf_counter()-t0}")
        self.pool.Release(L_el)

        # -(q,∇·ufrac)
        rhs_el = self.cfd_utils.ComputeElementwiseNodalDivergence(self.N, self.DN, vel_frac, out=self.pool.Get(1,(nelem,self.n_in_el)))
        xp.negative(rhs_el, out=rhs_el) # in-place negation to avoid allocation

        # (coef*dt/rho + tau)*(∇q,∇pold) = (coef*dt/rho + tau)*Lij*pj
        aux_scalar = self.cfd_utils.ApplyLaplacian(self.DN, pel, out=self.pool.Get(0,(nelem,self.n_in_el)))
        aux_scalar *= (p_factor * dt / self.rho)
        rhs_el += aux_scalar
        self.pool.Release(aux_scalar)

        # -tau*(∇q,Pi_pressure)
        if(self.clear_divergence_steps>0 or self.deactivate_pressure_stabilization == False):
            aux_scalar = self.cfd_utils.ComputePressureStabilizationProjectionTerm(self.N, self.DN, pres_proj_el,out=self.pool.Get(0,(nelem,self.n_in_el)))
            aux_scalar *= self.tau_1[:,xp.newaxis]
            rhs_el -= aux_scalar
            self.pool.Release(aux_scalar)

        # scale RHS elemental contributions by elemental volumes (integration)
        rhs_el *= self.elemental_volumes[:,xp.newaxis]

        # Assemble RHS
        rhs_p = self.pool.Get(2,(p.shape[0],)).ravel() #self.res_p
        rhs_p.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity, rhs_el, rhs_p)
        self.pool.Release(rhs_el)

        # Apply pressure BCs
        t0 = time.perf_counter()
        self.cfd_utils.ApplyHomogeneousDirichlet(self.fix_pres_indices, self.csr_data_indices, self.csr_diag_indices, self.L, rhs_p)
        #print(f"Time to apply pressure BCs: {time.perf_counter()-t0}")

        t0 = time.perf_counter()
        # Solve the system
        p, is_converged = self._SolvePressure(rhs_p,p)
        if not is_converged:
            print("CG failed to converge.")
        #print(f"Time to solve pressure: {time.perf_counter()-t0}")
        self.pool.Release(rhs_p)

        # Convert p back to xp for next steps
        # If USE_CUPY, p is already on GPU (cupy array). If numpy, it is numpy array.
        # Ensure it is wrapped as xp array just in case
        p = xp.asarray(p, dtype=cfd_utils.PRECISION)

        return p, is_converged

    def SolveStep3(self, vfrac, v_dirichlet, delta_p, dt, out=None):
        if out is None:
            raise ValueError("Output array must be provided to store the pressure projection values!")
        if out.shape != vfrac.shape:
            raise ValueError(f"Output array has incompatible shape {out.shape}, expected {(self.nnodes,self.dim)}")
        v = out

        #allocation of temporaries
        #tmp = self.tmp_elem #reusing
        nelem = self.DN.shape[0]

        # Get factor coming from the RK substep approximation of the pressure term
        time_scheme = self.time_scheme if self.time > self.startup_time else self.startup_time_scheme
        p_factor = self.GetPressureFactor(time_scheme)

        # Gather elemental pressure increments from the nodal values
        delta_p_el = self.ElemData(delta_p, self.connectivity, out=self.pool.Get(0,(nelem, self.n_in_el)))

        # Calculate the gradient of the pressure increment at the elemental level
        grad_dp_el = self.cfd_utils.ComputeDNN(self.N, self.DN, delta_p_el, out=self.pool.Get(1,self.DN.shape))
        grad_dp_el *= self.elemental_volumes[:,xp.newaxis,xp.newaxis]
        self.pool.Release(delta_p_el)

        # Assemble the gradient contributions at the nodes
        grad_dp_nodal = self.pool.Get(2,(delta_p.shape[0], self.dim))
        grad_dp_nodal.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity, grad_dp_el, grad_dp_nodal)
        self.pool.Release(grad_dp_el)

        # Update the fractional velocity with the pressure gradient contribution
        ##v = xp.empty(vfrac.shape, dtype=cfd_utils.PRECISION) #TODO: can we allocate this once and reuse it?
        v.ravel()[:] = vfrac.ravel() + (p_factor * dt / self.rho) * self.Minv * grad_dp_nodal.ravel()
        self.pool.Release(grad_dp_nodal)

        self.ApplyVelocitySlipConditions(v, self.normals)
        self.ApplyVelocityDirichletConditions(vfrac, v_dirichlet, 1.0, v)

        return v

    def SolveSolutionStep(self):

        # Get boundary condition data
        self.GetSlipIndices()
        self.GetFixityIndices()

        nelem = self.DN.shape[0]

        # Get boundary condition masks
        if (self.csr_data_indices is None or self.csr_diag_indices is None):
            self.csr_data_indices, self.csr_diag_indices = self.cfd_utils.GetScalarMatrixDirichletIndices(self.fix_pres_indices, self.L)

        # Obtain nodal data from Kratos
        self.v_adaptor.CollectData() #new vel
        v = xp.asarray(self.v_adaptor.data.reshape((len(self.model_part.Nodes),self.dim)), dtype=cfd_utils.PRECISION)

        self.v_adaptor_n.CollectData() #old vel
        vold = xp.asarray(self.v_adaptor_n.data.reshape((len(self.model_part.Nodes),self.dim)), dtype=cfd_utils.PRECISION)

        self.p_adaptor_n.CollectData() #old pressure
        pold = xp.asarray(self.p_adaptor_n.data, dtype=cfd_utils.PRECISION)

        self.b_adaptor.CollectData() #current body force
        b_new = xp.asarray(self.b_adaptor.data.reshape((len(self.model_part.Nodes),self.dim)), dtype=cfd_utils.PRECISION)

        self.b_adaptor_n.CollectData() #old body force
        b_old = xp.asarray(self.b_adaptor_n.data.reshape((len(self.model_part.Nodes),self.dim)), dtype=cfd_utils.PRECISION)

        # Get Dirichlet velocity values to be enforced
        v_new_dirichlet = (v.ravel()[self.fix_vel_indices]).copy()
        v_old_dirichlet = (vold.ravel()[self.fix_vel_indices]).copy()

        vold_backup = xp.empty_like(vold)
        pold_backup = xp.empty_like(pold)
        vnew = xp.empty_like(vold)

        # Perform substepping
        self.update_precond = True # Update preconditioner at the first substep
        current_time = self.time - self.dt
        while (self.time - current_time > 1.0e-12):
            if self.force_first_order_splitting == True:
                pold.fill(0.0)
                
            ##
            backup_current_time = current_time
            # Explicitly backup the state once per substep
            vold_backup[:] = vold
            pold_backup[:] = pold
            #vnew[:] = vold

            vmax = xp.max(xp.linalg.norm(vold, axis=1))

            # Guard against NaN and zero velocity (e.g., starting from rest)
            # NaN comparisons silently fail, and vmax=0 causes all velocities to trigger retry
            if not xp.isfinite(vmax) or vmax < 1e-15:
                vmax = 1e-15  # Set minimum threshold to avoid division by zero and false triggers

            repeat_step1 = True

            while(repeat_step1):
                print(f"========== Next Kratos time: {self.time:.3e} current substep time: {current_time:.3e} vmax= {vmax:.3f} CFL= {self.cfl:.2f}")
                # Compute convective operator spectral radius
                rho_conv = self.ComputeElementalConvectiveOperatorSpectralRadius(vold, out=self.pool.Get(2, (nelem,))) #self.tmp_n_elem)

                # Get maximum allowed time step with previous substep velocity
                max_dt = self._ComputeDeltaTime(rho_conv)
                self.pool.Release(rho_conv)

                # Compute current substep time step
                n_substeps = math.ceil((self.time - current_time) / max_dt)
                substep_dt = (self.time - current_time) / n_substeps
                current_time += substep_dt
                print(f"========== Substep dt: {substep_dt:.3e} - n_substeps left: {n_substeps} - % of step completion: {((current_time-self.time+self.dt)/self.dt*100.0):.1f}%")

                # Compute current substep Dirichlet velocity and body force values
                # Note that we assume a linear interpolation within the step
                substep_dt_factor = 1.0 - (self.time - current_time) / self.dt
                # print(f"\ttime: {self.time} - current_time: {current_time} - substep_dt: {substep_dt} - substep_dt_factor: {substep_dt_factor}")
                b = (1.0 - substep_dt_factor) * b_old + substep_dt_factor * b_new
                v_dirichlet = (1.0 - substep_dt_factor) * v_old_dirichlet + substep_dt_factor * v_new_dirichlet

                # Compute stabilization constants
                self.tau_1, self.tau_2 = self.ComputeTau(self.h, rho_conv, self.dyn_visc, self.rho, substep_dt)

                # Perform fractional step
                t1 = time.perf_counter()
                vfrac = self.SolveStep1(vold, v_dirichlet, pold, b, substep_dt)

                ##check if velocity increased significantly after step 1, if so, reduce CFL and repeat step 1
                if(self.clear_divergence_steps == 0): # if we are clearing divergence, we expect some velocity increase, so skip the check
                    vmax_after_step1 = xp.max(xp.linalg.norm(vfrac, axis=1))
                    # Guard against NaN: NaN comparisons silently return False, bypassing the retry logic
                    if not xp.isfinite(vmax_after_step1):
                        print(f"----- Warning: NaN detected in velocity after step 1, repeating step with reduced CFL")
                        repeat_step1 = True
                        self.cfl = self.cfl * 0.7
                        current_time = backup_current_time
                    elif(vmax_after_step1 > 1.25 * vmax):
                        print(f"----- Warning: local max velocity increased significantly after step 1: {vmax_after_step1} vs {vmax}")
                        repeat_step1 = True
                        ##vold[:] = vold_backup #not needed as it is not modified in step1
                        ##pold[:] = pold_backup
                        self.cfl = self.cfl * 0.7 # reduce CFL and repeat step 1
                        current_time = backup_current_time # reset current time to before step 1 to repeat it with the new CFL
                    else:
                        repeat_step1 = False
                        self.cfl = min(self.cfl*1.05, self.max_cfl) # increase CFL for next step if we are well below the limit
                else:
                    repeat_step1 = False

            self.step_1_total_time += time.perf_counter() - t1

            t2 = time.perf_counter()
            p, is_converged = self.SolveStep2(vfrac, pold, substep_dt)
            delta_p = p - pold #TODO: most probably we could modify p in place
            if is_converged==False:
                delta_p.fill(0.0) #anyhow we will redo the step
            self.step_2_total_time += time.perf_counter() - t2

            t3 = time.perf_counter()
            vnew = self.SolveStep3(vfrac, v_dirichlet, delta_p, substep_dt, out=vnew) # Note that vnew is assigned to vold in here for the next substep
            self.step_3_total_time += time.perf_counter() - t3

            # First time clear divergence
            vmax_after_step3 = xp.max(xp.linalg.norm(vnew, axis=1))
            # Guard against NaN: NaN comparisons silently return False, bypassing the retry logic
            if not xp.isfinite(vmax_after_step3):
                print(f"-----****!!!! Warning: NaN detected in velocity after step 3, redoing step")
                vnew[:] = vold_backup
                p[:] = pold_backup
                current_time = backup_current_time
                self.cfl = self.cfl * 0.5
            elif self.clear_divergence_steps > 0:
                self.clear_divergence_steps -= 1
                #if(self.clear_divergence_steps > 5):
                p.fill(0.0)
                pold.fill(0.0)
            elif(vmax_after_step3 > 1.25 * vmax or is_converged==False): #end of step velocity exploded or pressure solver did not converge - redo full step
                if(vmax_after_step3 > 1.25 * vmax):
                    print("end of step velocity blew up ... redoing full step")
                elif(is_converged==False):
                    print("pressure solve failed to converge ... redoing full step")
                vnew[:] = vold_backup
                p[:] = pold_backup
                current_time = backup_current_time # reset current time to before step 3 to repeat it with the new CFL
                self.cfl = self.cfl * 0.5 # reduce CFL and redo the step
                print(f"-----****!!!! Warning: local max velocity increased significantly after step 3: {vmax_after_step3} vs {vmax}")

            #for the next substep
            pold[:] = p
            vold [:] = vnew

            #self.outfile.write(str(current_time)+" " + str(np.linalg.norm(vold))+" "+str(np.linalg.norm(p))+"\n") #please do not remove this, it is used for testing purposes

        # Update Kratos database
        self.v_adaptor.data = cfd_utils.asnumpy(vold) # Note that vold is actually the vnew obtained in the last substep
        self.v_adaptor.StoreData()
        self.p_adaptor.data = cfd_utils.asnumpy(p)
        self.p_adaptor.StoreData()

    def ComputeVelocityResidual(self, v_elemental, p_elemental, b_elemental, proj_elemental, proj_div_el, DN, tau, out=None):
        if out is None:
            raise ValueError("Output array must be provided to store the assembled residual values!")
        if out.shape != (self.nnodes,self.dim):
            raise ValueError(f"Output array has incompatible shape {out.shape}, expected {(self.nnodes,self.dim)}")
        res_assembled = out

        ttot = time.perf_counter()

        #(w,b)
        t0 = time.perf_counter()
        res = self.cfd_utils.ComputeBodyForceContribution(b_elemental, out=self.pool.Get(0,v_elemental.shape))
        res *= self.rho
        #print(f"\t\ttime {1}: {time.perf_counter() - t0}")

        t0 = time.perf_counter()
        grad_v = self.cfd_utils.ComputeElementalGradient(DN, v_elemental, out=self.pool.Get(1,(v_elemental.shape[0], self.dim, self.dim)))
        a_grad_elemental = self.cfd_utils.ComputeElementalConvectiveOperator(v_elemental, grad_v, out=self.pool.Get(2,v_elemental.shape))
        self.pool.Release(grad_v)
        #print(f"\t\ttime {2}: {time.perf_counter() - t0}")

        t0 = time.perf_counter()
        # -(w,rho·a·∇u)
        convective = self.cfd_utils.ComputeConvectiveContribution(a_grad_elemental, out=self.pool.Get(1,v_elemental.shape)) #a_grad_elemental is actually still tmp_elem here
        convective *= self.rho
        res -= convective #assemble convective contribution
        self.pool.Release(convective)

        # -(rho·a·∇w,tau_1·rho·a·∇u) + (rho·a·∇w,tau_1·rho·Pi_conv)
        convective_stab = self.cfd_utils.ComputeMomentumStabilization(DN, v_elemental, a_grad_elemental, proj_elemental, out=self.pool.Get(1,v_elemental.shape)) #a_grad_elemental is actually still tmp_elem here
        convective_stab *= self.rho * self.rho
        res -= convective_stab * self.tau_1[:, None, None] #assemble convective stabilization contribution
        #print(f"\t\ttime {8}: {time.perf_counter() - t0}")
        self.pool.Release(convective_stab)
        self.pool.Release(a_grad_elemental)

        #-(∇w,∇u)
        t0 = time.perf_counter()
        res -= self.dyn_visc * self.cfd_utils.ApplyLaplacian(DN, v_elemental, out=self.pool.Get(1,v_elemental.shape))
        #print(f"\t\ttime {9}: {time.perf_counter() - t0}")
        self.pool.ReleaseByIndex(1)

        #+(∇·w,p_gauss)
        t0 = time.perf_counter()
        res += self.cfd_utils.ComputeDNN(self.N, DN, p_elemental, out=self.pool.Get(1,v_elemental.shape))
        #print(f"\t\ttime {10}: {time.perf_counter() - t0}")
        self.pool.ReleaseByIndex(1)

        #-(∇·w,tau_2·∇·v) + (∇·w,tau_2·Pi_div)
        t0 = time.perf_counter()
        div_div_stab = self.cfd_utils.ComputeDivDivStabilization(self.N, DN, v_elemental, proj_div_el, out=self.pool.Get(1,v_elemental.shape))
        res -= div_div_stab * self.tau_2[:,np.newaxis,np.newaxis]
        #print(f"\t\ttime {11}: {time.perf_counter() - t0}")
        self.pool.Release(div_div_stab)

        # Multiply by volume
        t0 = time.perf_counter()
        res *= self.elemental_volumes[:,xp.newaxis,xp.newaxis]
        #print(f"\t\ttime {12}: {time.perf_counter() - t0}")

        # Vectorized assembly of the elemental contributions into the nodal residual
        t0 = time.perf_counter()
        res_assembled.fill(0.0)
        self.cfd_utils.AssembleVector(self.connectivity, res, res_assembled)
        self.pool.Release(res)
        #print(f"\t\ttime {13}: {time.perf_counter() - t0}")

        #print(f"\t Time inside ComputeVelocityResidual: {time.perf_counter() - ttot}")
        # return res_assembled

    def ComputeVelocityResidualWall(self, v_elemental_wall, out=None):
        if out is None:
            raise ValueError("Output array must be provided to store the assembled residual values!")
        if out.shape != (self.nnodes,self.dim):
            raise ValueError(f"Output array has incompatible shape {out.shape}, expected {(self.nnodes,self.dim)}")
        res_assembled = out

        # Compute the viscous wall law contribution
        # We use the 0 position of the pool as it is ensured to be sufficient as the number of wall elements is always less than the number of volume elements
        slip_length = self.wall_mp.ProcessInfo[KratosFluid.SLIP_LENGTH] #TODO: think on how to handle this properly
        res_wall_shape = v_elemental_wall.shape
        res_wall = self.cfd_utils.ComputeNavierSlipWallContribution(v_elemental_wall, self.wall_tangential_projector, self.dyn_visc, slip_length, out=self.pool.Get(0,res_wall_shape))
        res_wall *= self.wall_areas[:, None, None]

        # Vectorized assembly of the wall contributions into the nodal residual
        # Note that this is assumed to be done on top of the volume residual (no innitialization of the assembled residual is done)
        self.cfd_utils.AssembleVector(self.connectivity_wall, res_wall, res_assembled)
        self.pool.Release(res_wall)

    def GetStep(self):
        return self.model_part.ProcessInfo[KM.STEP]

    @classmethod
    def _GetDefaultSolvingSettings(self):

        default_settings = KM.Parameters("""{
            "solver_type": "none",
            "model_part_name": "",
            "domain_size": -1,
            "material_import_settings": {
                "materials_filename": ""
            },
            "echo_level": 0,
            "compute_reactions": false,
            "linear_solver_settings": {
                "max_iteration" : 200,
                "tolerance" : 1e-6
            },
            "clear_divergence_steps" : 1,
            "time_stepping" : {
                "time_step" : 0.1,
                "max_cfl_number" : 0.5,
                "max_fourier_number" : 0.5
            },
            "time_scheme" : "RK4",
            "startup_time" : 0.0,
            "startup_time_scheme" : "FE"
        }""")
        return default_settings

    def _ComputeFourier(self, dt):
        return xp.max(dt * self.dyn_visc / self.rho / self.h**2)

    def _ComputeDeltaTime(self, rho_conv):
        """Compute the time step restricted by CFL, Fourier, and user-defined limits.

        This method calculates the maximum allowable time step based on:
        1. CFL condition from the convective operator spectral radius
        2. Fourier condition from viscous diffusion (precomputed)
        3. User-defined maximum time step

        Parameters
        ----------
        rho_conv : ndarray
            Elemental convective operator spectral radius array (shape: nelem).

        Returns
        -------
        float
            The most restrictive (minimum) time step from all conditions.

        Raises
        ------
        ValueError
            If rho_conv contains negative values or is empty.
        """
        # Validate input
        if rho_conv.size == 0:
            raise ValueError("rho_conv array is empty; cannot compute CFL time step.")
        if xp.any(rho_conv < 0):
            raise ValueError("rho_conv contains negative values; spectral radius must be non-negative.")

        # Compute CFL-limited time step with numerical safeguard against division by zero
        # Use machine epsilon relative to the data type for robustness
        eps = xp.finfo(rho_conv.dtype).eps

        # Reuse preallocated scratch array (self.tmp_n_elem) to avoid temporary allocation
        # This 1D array has shape (nelem,) matching rho_conv, eliminating the intermediate array
        # that would otherwise be created by self.max_cfl / (rho_conv + eps)
        # xp.divide(self.max_cfl, (rho_conv + eps), out=self.tmp_n_elem)
        # dt_cfl = xp.min(self.tmp_n_elem)

        cfl_factor = 1.0
        if self.clear_divergence_steps > 0:
            cfl_factor = 0.1
            print("Applying stricter CFL condition for divergence clearance step")

        #much better version avoiding vector divisions and mostly working with scalars
        max_rho_conv = xp.max(rho_conv)
        dt_cfl = cfl_factor*self.cfl / (max_rho_conv + eps)

        # Return the most restrictive time step from all conditions
        return float(min(min(dt_cfl, self.dt_fourier), self.dt))

    def _SolvePressure(self, rhs, previous_p):
        if USE_CUPY:

            #

            #set a stricter tolerance during the clear divergence steps.
            if(self.clear_divergence_steps > 0):
                #precond = JacobiPreconditioner(self.L)
                factor = 1e-6
                exact_res_recalculation=1
                self.preconditioner.update_matrix_values(self.L)
                precond = self.preconditioner.aspreconditioner()
                print("doing div cleareance step")
                t0 = time.perf_counter()
                #sol, status = cfd_utils.sparse_linalg.cg(self.L, rhs, x0=previous_p, rtol=factor*self.pressure_tolerance, M=precond,maxiter=2000)
                sol, status = self.cfd_utils.robust_cg(self.L, rhs, x0=None, rtol=factor*self.pressure_tolerance, atol=0.0, M=precond, maxiter=500,xp=xp)
                print("status=",status)
            else:
                t0 = time.perf_counter()
                self.preconditioner.update_matrix_values(self.L)
                print(f"AMG graph update time: {time.perf_counter() - t0:.4f} seconds")
                precond = self.preconditioner.aspreconditioner()

                factor = 1.0
                exact_res_recalculation=10
                t0 = time.perf_counter()
                #sol, status = cfd_utils.sparse_linalg.cg(self.L, rhs, x0=previous_p, rtol=factor*self.pressure_tolerance, M=precond,maxiter=self.pressure_max_iteration)
                sol, status = self.cfd_utils.robust_cg(self.L, rhs, x0=previous_p, rtol=factor*self.pressure_tolerance, atol=0.0, M=precond, maxiter=self.pressure_max_iteration,xp=xp)
            print(f"CG solve time: {time.perf_counter() - t0:.4f} seconds")
            is_converged = status == 0 # Note that the status is 0 if the solver converged
            if not is_converged:
                print("************************************************************************")
                print("CG failed to converge.")
                print("************************************************************************")
            return sol, is_converged
        else:
            #FIXME: USE AMGCL FROM FUTURE NAMESPACE HERE!
            A_np = cfd_utils.sparse.csr_matrix((
                xp.asarray(self.L.value_data(), dtype=cfd_utils.PRECISION),
                xp.asarray(self.L.index2_data()),
                xp.asarray(self.L.index1_data())),
                shape=(self.L.Size1(), self.L.Size2()))

            sol, status = self.cfd_utils.robust_cg(A_np, rhs, rtol=self.pressure_tolerance, xp=xp)
            is_converged = status == 0 # Note that the status is 0 if the solver converged

            return sol, is_converged

            # pass #FIXME: we need the future linear solvers here
            # x = KM.SystemVector(self.L.Size1()) #FIXME: fake arrays to match the interface
            # x.SetValue(0.0)
            # self.linear_solver.Solve(self.L, rhs, x)
            # return x, self.linear_solver.IsConverged()

    def _FinalizePressureLinearSolver(self):
        pass

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    simulation = VectorizedCFDStage(model,parameters)
    simulation.Run()
