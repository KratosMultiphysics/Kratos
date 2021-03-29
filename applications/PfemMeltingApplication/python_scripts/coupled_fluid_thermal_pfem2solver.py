from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff
import KratosMultiphysics.MeshingApplication as MeshApp
import KratosMultiphysics.PFEM2Application as PFEM2
import KratosMultiphysics.PfemMeltingApplication as PfemM



# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    
    return PfemCoupledFluidThermalSolver(main_model_part, custom_settings)

class PfemCoupledFluidThermalSolver(PythonSolver):
     
    @classmethod
    def GetDefaultParameters(cls):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "ThermallyCoupled1",
            "domain_size" : -1,
            "echo_level": 0,
            "laser_import_settings": {
                        "laser_filename": "LaserSettings.json"
                },
            "material_settings": {
                        "material_filename": "MateralCharacterization.json"
                },
            "fluid_solver_settings": {
                "solver_type": "navier_stokes_solver_vmsmonolithic",
                "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
                }
            },
            "thermal_solver_settings": {
                "solver_type": "Transient",
                "analysis_type": "linear",
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                        "materials_filename": "ThermalMaterials.json"
                }
            }
        }
        """)


        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        
        super(PfemCoupledFluidThermalSolver, self).__init__(model, custom_settings)
        
        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"],"OpenMP")
        


        from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion
        self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model,self.settings["thermal_solver_settings"],"OpenMP")


	#Laser settings                
        #materials_filename =self.settings["laser_import_settings"]["laser_filename"]

        materials_filename = self.settings["laser_import_settings"]["laser_filename"].GetString()
        
        material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)

        with open(self.settings["laser_import_settings"]["laser_filename"].GetString(), 'r') as parameter_file:
                materials = KratosMultiphysics.Parameters(parameter_file.read())  


        mat = materials["properties"][0]["Material"]

        self.variables = []
        self.values = []
        for key, value in mat["Variables"].items():
            var = KratosMultiphysics.KratosGlobals.GetVariable(key)
            self.variables.append(var)
            self.values.append(value)

        for key, table in mat["Tables"].items():
            table_name = key

            input_var = KratosMultiphysics.KratosGlobals.GetVariable(table["input_variable"].GetString())
            output_var = KratosMultiphysics.KratosGlobals.GetVariable(table["output_variable"].GetString())
            self.new_table = KratosMultiphysics.PiecewiseLinearTable()

            for i in range(table["data"].size()):
                self.new_table.AddRow(table["data"][i][0].GetDouble(), table["data"][i][1].GetDouble())
            
            #self.model.Properties.SetTable(input_var,output_var,new_table)

        materials_filename = self.settings["material_settings"]["material_filename"].GetString()

        material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)

        with open(self.settings["material_settings"]["material_filename"].GetString(), 'r') as parameter_file:
                materials = KratosMultiphysics.Parameters(parameter_file.read())  

        mat = materials["properties"][0]["Material"]

        self.variables_aux = []
        self.values_aux = []
        for key, value in mat["Variables"].items():
            var = KratosMultiphysics.KratosGlobals.GetVariable(key)
            self.variables_aux.append(var)
            self.values_aux.append(value)
            #print("##########")
            #print(key)
            #print(value)
            #print(var)

        #print("@@@@@@@@@@@@@@@@")
        #print("@@@@@@@@@@@@@@@@")
        #print(self.variables_aux)     
        #print(self.values_aux)
        
        self.Mesher = MeshApp.TetGenPfemModeler()   

        self.modeler = KratosMultiphysics.ConnectivityPreserveModeler()

        self.Pfem2_apply_bc_process = PFEM2.Pfem2ApplyBCProcess(self.fluid_solver.main_model_part);

        self.node_erase_process = KratosMultiphysics.NodeEraseProcess(self.fluid_solver.main_model_part);

        self.Streamline = PfemM.Streamline()

        self.Pfem2Utils = PFEM2.Pfem2Utils()

        self.faceheatflux = PfemM.FaceHeatFlux()

        self.HeatSource = PfemM.HeatSource()



    def AddVariables(self):
        # Import the fluid and thermal solver variables. Then merge them to have them in both fluid and thermal solvers.
         
        self.fluid_solver.AddVariables()
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FREE_SURFACE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_BOUNDARY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FLUID)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_LAGRANGIAN_INLET)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ACTIVATION_ENERGY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ARRHENIUS_COEFFICIENT)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HEAT_OF_VAPORIZATION)



        self.thermal_solver.AddVariables()

        KratosMultiphysics.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part)

    def ImportModelPart(self):
        # Call the fluid solver to import the model part from the mdpa

        self.fluid_solver.ImportModelPart()
        

        # Save the convection diffusion settings
        convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)

        # Here the fluid model part is cloned to be thermal model part so that the nodes are shared
        #self.modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        if self.domain_size == 2:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "EulerianConvDiff2D",
                                      "ThermalFace2D2N")
        else:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "EulerianConvDiff3D",
                                      "ThermalFace3D3N")

        
        # Set the saved convection diffusion settings to the new thermal model part
        self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)
        
     
    def PrepareModelPart(self):
        

        self.fluid_solver.PrepareModelPart()

        self.thermal_solver.PrepareModelPart()
 
        
        for node in self.fluid_solver.main_model_part.Nodes:
            
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y,0,-10.0);

            if(node.Y<0.00001): #if(node.Y>0.99999 or node.Y<0.00001):
                node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0,1.0);

        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[0], self.values_aux[0].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[1], self.values_aux[1].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[2], self.values_aux[2].GetDouble(), self.fluid_solver.main_model_part.Nodes)


    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.thermal_solver.AddDofs()
        

    def CalculateViscosity(self):
        import math
        for node in self.fluid_solver.main_model_part.Nodes:
            rho = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
            T = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            #mu = AuxFunction(T)

            Tc = T - 273.0
            e=0.0
            mu=0.0
            if( Tc <= 25):
                mu = 1e6
            if (Tc > 25.0 and Tc <= 200.0):
                e=14.48 - 0.13858*200.0 + 5.5960e-4*200.0*200.0 - 7.8665e-7 * 200.0 * 200.0 * 200.0
                mu = 1e6*(200.0-Tc)/175 + pow(10,e)
            if(Tc > 200.0 and Tc <= 242.0):
                mu = 4949.015 *  math.exp(-0.00345*Tc)
            if(Tc > 242.0 and Tc <= 260.0):
                mu = 0.00000000847 * math.exp(0.1266*Tc) ###curva ajustada
            if(Tc > 260.0):
                mu = 20000.0
            node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY,0,mu/rho)


    def AdaptMesh(self):
        

        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,0,0.15);

        self.Streamline.RungeKutta4ElementbasedSI(self.fluid_solver.main_model_part,100)

        for node in (self.fluid_solver.main_model_part).Nodes:
            node.Set(KratosMultiphysics.TO_ERASE, False)


        self.Pfem2Utils.MarkNodesTouchingWall(self.fluid_solver.main_model_part, 3, 0.15)

        self.node_erase_process.Execute()

        self.fluid_solver.main_model_part.Conditions.clear()

        self.fluid_solver.main_model_part.Elements.clear()


        #Boussinesq__Boussinesq_hidden_=self.fluid_solver.main_model_part.GetSubModelPart("Boussinesq__Boussinesq_hidden_")
        #Boussinesq__Boussinesq_hidden_.Conditions.clear()	 


        parametersf=self.settings["fluid_solver_settings"]

        parameters=self.settings["thermal_solver_settings"]
        
        self.skin_parts_list = []
        
        if parametersf.Has("skin_parts"):
            self.skin_parts_list = parametersf["skin_parts"]

        for i in range(self.skin_parts_list.size()):
            body_model_part_name=self.skin_parts_list[i].GetString()
            body_model_part_name=self.fluid_solver.main_model_part.GetSubModelPart(body_model_part_name)
            body_model_part_name.Conditions.clear() 
            body_model_part_name.Elements.clear() 
            body_model_part_name.Nodes.clear() 
        self.bodies_parts_list = []
        if parameters.Has("processes_sub_model_part_list"):
            self.bodies_parts_list = parameters["processes_sub_model_part_list"]
             
        self.bodies_parts_list.size()  
        for i in range(self.bodies_parts_list.size()):
            body_model_part_name=self.bodies_parts_list[i].GetString()
            body_model_part_name=self.fluid_solver.main_model_part.GetSubModelPart(body_model_part_name)
            body_model_part_name.Conditions.clear()
            body_model_part_name.Elements.clear() 
            body_model_part_name.Nodes.clear() 


        ##NEW FOR THERMAL SOLVER
        self.skin_parts_listaux = []
        if parametersf.Has("skin_parts"):
            self.skin_parts_listaux = parametersf["skin_parts"]
        for i in range(self.skin_parts_listaux.size()):
            body_model_part_name=self.skin_parts_listaux[i].GetString()
            body_model_part_name=self.thermal_solver.main_model_part.GetSubModelPart(body_model_part_name)
            body_model_part_name.Conditions.clear() 
        self.bodies_parts_listaux = []
        if parameters.Has("processes_sub_model_part_list"):
            self.bodies_parts_listaux = parameters["processes_sub_model_part_list"]
             
        self.bodies_parts_listaux.size()  
        for i in range(self.bodies_parts_listaux.size()):
            body_model_part_name=self.bodies_parts_listaux[i].GetString()
            body_model_part_name=self.thermal_solver.main_model_part.GetSubModelPart(body_model_part_name)
            body_model_part_name.Conditions.clear()


        (self.Mesher).ReGenerateMesh("LagrangianFluidVMS3D","ThermalFace3D3N", self.fluid_solver.main_model_part, self.node_erase_process, True, False, 1.4, 0.4)  

        #LagrangianFluidVMS3D
        #VMS3D 
        neighbor_search = KratosMultiphysics.FindNodalNeighboursProcess(self.fluid_solver.main_model_part)
        neighbor_search.Execute()
        neighbor_elements_search = KratosMultiphysics.FindElementalNeighboursProcess(self.fluid_solver.main_model_part, 3, 20)
        neighbor_elements_search.Execute()
        neighbor_condition_search = KratosMultiphysics.FindConditionsNeighboursProcess(self.fluid_solver.main_model_part,3, 20)
        neighbor_condition_search.Execute()


        (self.Pfem2_apply_bc_process).Execute();

        Parts_Parts_Auto1=self.fluid_solver.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
        Parts_Parts_Auto1.Conditions.clear()	 
        Parts_Parts_Auto1.Elements.clear()
        Parts_Parts_Auto1.Nodes.clear()



        fluid_computational_model_part=self.fluid_solver.main_model_part.GetSubModelPart("fluid_computational_model_part")
        fluid_computational_model_part.Conditions.clear()	 
        fluid_computational_model_part.Elements.clear()
        fluid_computational_model_part.Nodes.clear()

        

        for node in self.fluid_solver.main_model_part.Nodes:
            fluid_computational_model_part.AddNode(node,0)

        for elem in self.fluid_solver.main_model_part.Elements:
            fluid_computational_model_part.AddElement(elem,0) 

        
        fluid_computational_model_part.ProcessInfo = self.fluid_solver.main_model_part.ProcessInfo
        
        fluid_computational_model_part.Properties  = self.fluid_solver.main_model_part.Properties

        self.thermal_solver.main_model_part.Conditions.clear()	
        self.thermal_solver.main_model_part.Elements.clear()	
        self.thermal_solver.main_model_part.Nodes.clear()
        


        thermal_computing_domain=self.thermal_solver.main_model_part.GetSubModelPart("thermal_computing_domain")
        thermal_computing_domain.Conditions.clear()	 
        thermal_computing_domain.Elements.clear()
        thermal_computing_domain.Nodes.clear()


        Parts_Parts_Auto1=self.thermal_solver.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
        Parts_Parts_Auto1.Conditions.clear()	 
        Parts_Parts_Auto1.Elements.clear()
        Parts_Parts_Auto1.Nodes.clear()

        
        if self.domain_size == 2:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part, "EulerianConvDiff2D", "ThermalFace2D2N")
        else:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part,"EulerianConvDiff3D","ThermalFace3D3N")


        self.thermal_solver.main_model_part.Conditions.clear()

        
        for node in self.thermal_solver.main_model_part.Nodes:
            thermal_computing_domain.AddNode(node,0)

        for elem in self.thermal_solver.main_model_part.Elements:
            thermal_computing_domain.AddElement(elem,0) 

        thermal_computing_domain.Properties  = fluid_computational_model_part.Properties

        self.thermal_solver.main_model_part.Conditions.clear()	
        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(thermal_computing_domain, self.fluid_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS)
        transfer_process.Execute()

        neighbor_searcht = KratosMultiphysics.FindNodalNeighboursProcess(thermal_computing_domain)
        neighbor_searcht.Execute()
        neighbor_elements_searcht = KratosMultiphysics.FindElementalNeighboursProcess(thermal_computing_domain, 3, 20)
        neighbor_elements_searcht.Execute()
        neighbor_condition_searcht = KratosMultiphysics.FindConditionsNeighboursProcess(thermal_computing_domain,3, 20)
        neighbor_condition_searcht.Execute()


        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y,0,-10.0);

       

        pass

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        
        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal)

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.thermal_solver.Initialize()
        
    def Clear(self):
        
        (self.fluid_solver).Clear()
        (self.thermal_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):

       
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.fluid_solver.AdvanceInTime(current_time)
        return new_time

    def InitializeSolutionStep(self):

        self.Streamline.RungeKutta4ElementbasedSI(self.fluid_solver.main_model_part,100)

        self.AdaptMesh()

        #self.fluid_solver.Clear()
        #self.thermal_solver.Clear()

        #self.fluid_solver.Initialize()
        #self.thermal_solver.Initialize() 

        self.fluid_solver.InitializeSolutionStep()
        self.thermal_solver.InitializeSolutionStep()


    def Predict(self):
        
        self.fluid_solver.Predict()
        self.thermal_solver.Predict()
        
    def SolveSolutionStep(self):

        
        fluid_is_converged = self.fluid_solver.SolveSolutionStep()
        self.Streamline.RungeKutta4ElementbasedSI(self.fluid_solver.main_model_part,100)
        

        for node in self.fluid_solver.main_model_part.Nodes:
            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY, velocity)

        for node in self.fluid_solver.main_model_part.Nodes:
             node.SetSolutionStepValue(KratosMultiphysics.FACE_HEAT_FLUX,0,0.0);
             
        #print(self.variables)  
        #print(self.values)  

        #print("22222222222222222222222222222222222222")
        #print("22222222222222222222222222222222222222")
        #print("22222222222222222222222222222222222222")

        #x=self.variables[2]
        #y=self.values[2]

        #print(x)
        #print(y)
        #sssssssssssssssssssssss 

        x=self.values[1].GetDouble()
        y=self.values[2].GetDouble()
        z=self.values[3].GetDouble()
        radius=self.values[4].GetDouble()
        q=self.values[0].GetDouble()
        
        #print("2222222222222222222222222222222222222222")
        #print("2222222222222222222222222222222222222222")
        #print(self.new_table)
        #print(self.fluid_solver.main_model_part)
        #prop= self.fluid_solver.main_model_part.Properties
        #self.fluid_solver.main_model_part.Properties.SetTable("TEMPERATUVE","VISCOSITY",new_table) 
        #print(self.fluid_solver.main_model_part)
        #ssssssssssssssssss
        #for node in self.fluid_solver.main_model_part.Nodes:
        #     node.SetSolutionStepValue(KratosMultiphysics.FACE_HEAT_FLUX,0,q);
        #ssssssssssssssss

        #self.faceheatflux.FaceHeatFluxDistribution(self.fluid_solver.main_model_part, 0.57102, 0.99997, 0.5299, 0.8, 10000000.0)

        #print(x)

        self.faceheatflux.FaceHeatFluxDistribution(self.fluid_solver.main_model_part, x, y, z, radius, q)
        self.HeatSource.Heat_Source(self.fluid_solver.main_model_part)
        #sssssssssssssss 
        #sssssssssssssssssss 
        thermal_is_converged = self.thermal_solver.SolveSolutionStep()
        
        self.CalculateViscosity()	

        return (fluid_is_converged and thermal_is_converged)

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

    def Solve(self):
         
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()
