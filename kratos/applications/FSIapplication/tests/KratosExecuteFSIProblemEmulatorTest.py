from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#~ from KratosMultiphysics import *
#~ from KratosMultiphysics.FSIApplication import *
#~ from KratosMultiphysics.SolidMechanicsApplication import *

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as KratosALE
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

import process_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

class KratosExecuteFSIProblemEmulatorTest(KratosUnittest.TestCase):

    def __init__(self, ProjectParameters):
        
        self.ProjectParameters = ProjectParameters
        
        self.vector_space = KratosMultiphysics.UblasSparseSpace()
    
        self.structure_main_model_part = KratosMultiphysics.ModelPart(self.ProjectParameters["structure_solver_settings"]["problem_data"]["model_part_name"].GetString())
        self.structure_main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.ProjectParameters["structure_solver_settings"]["problem_data"]["domain_size"].GetInt())

        SolidModel = {ProjectParameters["structure_solver_settings"]["problem_data"]["model_part_name"].GetString() : self.structure_main_model_part}
        
        # Construct the structure solver
        structure_solver_module = __import__(self.ProjectParameters["structure_solver_settings"]["solver_settings"]["solver_type"].GetString())
        self.structure_solver = structure_solver_module.CreateSolver(self.structure_main_model_part, 
                                                                     self.ProjectParameters["structure_solver_settings"]["solver_settings"])
        print("* Structure solver constructed.")
        
        # Construct the coupling partitioned strategy
        import convergence_accelerator_factory     
        self.coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(self.ProjectParameters["coupling_solver_settings"]["solver_settings"]["coupling_strategy"])
        print("* Coupling strategy constructed.")
        
        self.structure_solver.AddVariables()
        
        self.structure_solver.ImportModelPart()
        
        self.structure_solver.AddDofs()     
        
        # Get the structure process list
        for i in range(self.ProjectParameters["structure_solver_settings"]["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = self.ProjectParameters["structure_solver_settings"]["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            SolidModel.update({part_name: self.structure_main_model_part.GetSubModelPart(part_name)})

        # Structure processes construction    
        import process_factory
        self.list_of_processes = process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( self.ProjectParameters["structure_solver_settings"]["constraints_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( self.ProjectParameters["structure_solver_settings"]["loads_process_list"] )

        # Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()
        
        # Structure solver initialization
        self.structure_solver.Initialize()
        self._SetStructureNeumannCondition()
        
        # Coupling utility initialization
        self.coupling_utility.Initialize()

        
    def Solve(self):
        
        self.structure_solver.SolverInitialize()

        # Stepping and time settings
        Dt = self.ProjectParameters["structure_solver_settings"]["problem_data"]["time_step"].GetDouble()
        end_time = self.ProjectParameters["structure_solver_settings"]["problem_data"]["end_time"].GetDouble()
        
        # Coupling convergence parameters
        nl_tol = self.ProjectParameters["coupling_solver_settings"]["solver_settings"]["nl_tol"].GetDouble()
        max_nl_it = self.ProjectParameters["coupling_solver_settings"]["solver_settings"]["nl_max_it"].GetInt()
                
        time = 0.0
        step = 0
        
        residual_size = self._GetInterfaceProblemSize()*2                   # Interface DOFs number times PROBLEM_SIZE
        self.iteration_value = KratosMultiphysics.Vector(residual_size)     # Interface solution guess (it might be velocity or fluxes depending on the type of coupling) 

        for i in range(0,residual_size):
            self.iteration_value[i] = 0.0

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()
            
        while(time <= end_time):

            time = time + Dt
            step = step + 1
            
            self.structure_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME_STEPS, step)
            
            self.structure_main_model_part.CloneTimeStep(time)    

            print("STEP = ", step)
            print("TIME = ", time)

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()
    
            self.structure_solver.SolverInitializeSolutionStep()
            self.structure_solver.SolverPredict()
            
            self.coupling_utility.InitializeSolutionStep()
            
            for nl_it in range(1,max_nl_it+1):
                self.structure_solver.main_model_part.ProcessInfo[KratosFSI.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it  ###### USELESS?¿?¿
                
                self.coupling_utility.InitializeNonLinearIteration()
                
                # Residual computation                
                disp_residual = self._ComputeDirichletNeumannResidual()
                nl_res_norm = self.vector_space.TwoNorm(disp_residual)
                                    
                # Check convergence
                if nl_res_norm < nl_tol:
                    break 
                    
                else:
                    # If convergence is not achieved, perform the correction of the prediction
                    self.coupling_utility.UpdateSolution(disp_residual, self.iteration_value)
                    self.coupling_utility.FinalizeNonLinearIteration()
                
            self.structure_solver.SolverFinalizeSolutionStep()
            self.coupling_utility.FinalizeSolutionStep()
                    
            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()
        
            # Unitcest convergence criterion check
            self.assertLess(nl_res_norm, nl_tol)
            
        for process in self.list_of_processes:
            process.ExecuteFinalize()


    ### PRIVATE METHODS SECTION ###
    def _GetInterfaceProblemSize(self):
        
        # Get the structure interface problem size
        interface_submodelpart_name = self.ProjectParameters["coupling_solver_settings"]["solver_settings"]["structure_interfaces_list"][0].GetString()
        structure_interface_pb_size = len(self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name).Nodes)
        
        return structure_interface_pb_size
        
        
    def _SetStructureNeumannCondition(self):
        
        structure_computational_submodelpart = self.structure_solver.GetComputingModelPart()
                
        aux_count = 0
        for cond in self.structure_solver.main_model_part.Conditions:
            if(cond.Id > aux_count):
                aux_count = cond.Id
                                
        interface_submodelpart_name = self.ProjectParameters["coupling_solver_settings"]["solver_settings"]["structure_interfaces_list"][0].GetString()
        interface_submodelpart_i = self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name)
    
        # Create the point load condition in the structure interface submodelpart
        for node in interface_submodelpart_i.Nodes:
            aux_count+=1
            structure_computational_submodelpart.CreateNewCondition("PointLoadCondition2D1N",aux_count,[node.Id],self.structure_solver.main_model_part.Properties[0]) 
    

    def _ComputeDirichletNeumannResidual(self):
        
        interface_submodelpart_name = self.ProjectParameters["coupling_solver_settings"]["solver_settings"]["structure_interfaces_list"][0].GetString()
    
        K = 1000.0 # Spring stiffness
        
        # Impose the spring reactions that emulate the fluid load over the structure interface
        i = 0
        for node in self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name).Nodes:
            point_load = KratosMultiphysics.Vector(3)
            point_load[0] = - K*self.iteration_value[i]
            point_load[1] = - K*self.iteration_value[i+1]
            point_load[2] = 0.0
            
            node.SetSolutionStepValue(KratosSolid.POINT_LOAD, 0, point_load)
            
            i += 2
            
        # Solve structure problem
        self.structure_solver.SolverSolveSolutionStep()
            
        # Compute the displacement residual
        disp_residual = KratosMultiphysics.Vector(self._GetInterfaceProblemSize()*2)
        
        i = 0       
        for node in self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name).Nodes:
            vector_projected = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0)
        
            disp_residual[i] = vector_projected[0] - self.iteration_value[i]
            disp_residual[i+1] = vector_projected[1] - self.iteration_value[i+1]
            i+=2
            
        return disp_residual
