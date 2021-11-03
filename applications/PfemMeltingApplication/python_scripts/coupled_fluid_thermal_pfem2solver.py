from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff
import KratosMultiphysics.MeshingApplication as MeshApp
#import KratosMultiphysics.PFEM2Application as PFEM2
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
                        "material_filename": "FluidMaterials.json"
                },
            "environment_settings" : {
	        "gravity": [0, 0, 0],
		"ambient_temperature" : 0.15
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

        self.readenvironmentSettings()         

        self.readLasserSettings()
	#Laser settings                

        self.readMaterialCharacterization()

        
        self.Mesher = MeshApp.TetGenPfemModeler()   

        self.modeler = KratosMultiphysics.ConnectivityPreserveModeler()

        self.PfemM_apply_bc_process = PfemM.PfemMeltingApplyBCProcess(self.fluid_solver.main_model_part);

        self.node_erase_process = KratosMultiphysics.NodeEraseProcess(self.fluid_solver.main_model_part);

        self.Streamline = PfemM.Streamline()

        #self.Pfem2Utils = PFEM2.Pfem2Utils()

        self.faceheatflux = PfemM.FaceHeatFlux()

        self.HeatSource = PfemM.HeatSource()

    def readenvironmentSettings(self):

        with open("ProjectParameters.json",'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())


        self.gravity = []
	#self.gravity= [] 
        ambient_temperature=project_parameters["problem_data"]["environment_settings"]["ambient_temperature"]
        self.gravity=project_parameters["problem_data"]["environment_settings"]["gravity"]
        
 

        #materials_filename = self.settings["laser_import_settings"]["laser_filename"].GetString()
        
        #material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)


        #open(self.settings["environment_settings"]) as parameter_file:
        #materials = KratosMultiphysics.Parameters(self.settings["environment_settings"].read())   ##laser_settings


        #materials_filename = self.settings["environment_settings"]["ambient_temperature"]

        #materials = KratosMultiphysics.Parameters(materials_filename.read())

    def readLasserSettings(self):

        materials_filename = self.settings["laser_import_settings"]["laser_filename"].GetString()
        
        material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)

        with open(self.settings["laser_import_settings"]["laser_filename"].GetString(), 'r') as parameter_file:
                materials = KratosMultiphysics.Parameters(parameter_file.read())   ##laser_settings


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


    def readMaterialCharacterization(self):

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
            #print("varrrrrrrrrrrrrrrrrr")
            #print("varrrrrrrrrrrrrrrrrr")
            #print(var)
           
        #The part below reads the temperature dependent viscosity
        table1=mat["Tables"]

        taux1= table1["Table1"]
        input_var = KratosMultiphysics.KratosGlobals.GetVariable(taux1["input_variable"].GetString())
        output_var = KratosMultiphysics.KratosGlobals.GetVariable(taux1["output_variable"].GetString())

        self.new_table_aux = KratosMultiphysics.PiecewiseLinearTable()

        for i in range(taux1["data"].size()):
            self.new_table_aux.AddRow(taux1["data"][i][0].GetDouble(), taux1["data"][i][1].GetDouble())

        #print(self.new_table_aux) 
        


        #taux1= table1["Table2"]
        #input_var = KratosMultiphysics.KratosGlobals.GetVariable(taux1["input_variable"].GetString())
        #output_var = KratosMultiphysics.KratosGlobals.GetVariable(taux1["output_variable"].GetString())

        #self.new_table_aux2 = KratosMultiphysics.PiecewiseLinearTable()
        

        #for i in range(taux1["data"].size()):
        #    self.new_table_aux2.AddRow(taux1["data"][i][0].GetDouble(), taux1["data"][i][1].GetDouble())

        #taux1= table1["Table3"]
        #input_var = KratosMultiphysics.KratosGlobals.GetVariable(taux1["input_variable"].GetString())
        #output_var = KratosMultiphysics.KratosGlobals.GetVariable(taux1["output_variable"].GetString())

        #self.new_table_aux3 = KratosMultiphysics.PiecewiseLinearTable()

        #for i in range(taux1["data"].size()):
        #    self.new_table_aux3.AddRow(taux1["data"][i][0].GetDouble(), taux1["data"][i][1].GetDouble())



    

    def AddVariables(self):
        # Import the fluid and thermal solver variables. Then merge them to have them in both fluid and thermal solvers.
         
        self.fluid_solver.AddVariables()
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FREE_SURFACE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_BOUNDARY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FLUID)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_LAGRANGIAN_INLET)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_INTERFACE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.RADIATIVE_INTENSITY)

        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ACTIVATION_ENERGY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ARRHENIUS_COEFFICIENT)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HEAT_OF_VAPORIZATION)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ARRHENIUS_VALUE)



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
 
        #print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        #print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        #print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        #print (self.fluid_solver)
        #print("#####################################")
        #print("#####################################")
        #print("#####################################")


        for node in self.fluid_solver.main_model_part.Nodes:  
            node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 0.0) 
        parametersf=self.settings["fluid_solver_settings"]

        self.skin_parts_list = []
        
        if parametersf.Has("skin_parts"):
            self.skin_parts_list = parametersf["skin_parts"]

        for i in range(self.skin_parts_list.size()):
            body_model_part_name=self.skin_parts_list[i].GetString()
            if(body_model_part_name=="NoSlip3D_No_Slip_Auto1"):
                body_model_part_name=self.fluid_solver.main_model_part.GetSubModelPart(body_model_part_name)
                for node in body_model_part_name.Nodes:
                    #print("kjhkjhkjhkjhkh")
                    node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 1.0) #NODES NOT MOVE ARE CONSIDERED PART OF THE STRUCTURE
 
                  
        
        #for node in self.fluid_solver.main_model_part.Nodes:
        #    if(node.Y<0.00001): #if(node.Y>0.99999 or node.Y<0.00001):
        #        node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0,1.0);



        #for node in self.fluid_solver.main_model_part.Nodes:  
        #    node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 0.0) 
        #    if(node.Y<0.0001):
        #        node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 1.0) 
                #node.Fix(KratosMultiphysics.VELOCITY_X)
                #node.Fix(KratosMultiphysics.VELOCITY_Y)
                #node.Fix(KratosMultiphysics.VELOCITY_Z)


            #aux=node.X * node.X + node.Z * node.Z 
            #radious=aux**.5
            #if (radious <=0.0018):
            #    node.SetSolutionStepValue(KratosMultiphysics.IS_INTERFACE,0, 1.0) 
            #    node.Fix(KratosMultiphysics.VELOCITY_X)
            #    node.Fix(KratosMultiphysics.VELOCITY_Y)
            #    node.Fix(KratosMultiphysics.VELOCITY_Z)

            #if (radious <=0.0017):
            #    node.SetSolutionStepValue(KratosMultiphysics.CONDUCTIVITY,0,386.4) 
            #    node.SetSolutionStepValue(KratosMultiphysics.SPECIFIC_HEAT,0,390.0) 
            #    node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,8920) 
            #    #node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 1.0) 

            #aux=node.X * node.X + node.Z * node.Z 
            #radious=aux**.5
            #if(node.X>0.0032931):
            #    XP=0.013-node.X
            #    ZP=0-node.Z
            #    r=(XP * XP + ZP * ZP)**.5
            #    if (r < 0.00177):
            #        node.SetSolutionStepValue(KratosMultiphysics.IS_INTERFACE,0, 1.0) 
            #        node.Fix(KratosMultiphysics.VELOCITY_X)
            #        node.Fix(KratosMultiphysics.VELOCITY_Y)
            #        node.Fix(KratosMultiphysics.VELOCITY_Z)
            #    if (r < 0.0016):
            #        node.SetSolutionStepValue(KratosMultiphysics.CONDUCTIVITY,0,386.4) 
            #        node.SetSolutionStepValue(KratosMultiphysics.SPECIFIC_HEAT,0,390.0) 
            #        node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,8920) 
            #        #node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 1.0) 

            #    #node.SetSolutionStepValue(IS_INTERFACE,0, 1.0)

    



        #KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[0], self.values_aux[0].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        #KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[1], self.values_aux[1].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        #KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[2], self.values_aux[2].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        self.assign_nodally_properties();
        self.ReMesh()

    def assign_nodally_properties(self):
        #here we assign ACTIVATION_ENERGY, ARRHENIUS_COEFFICIENT and HEAT_OF_VAPORIZATION taken from FluidMaterial.json 
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[0], self.values_aux[0].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[1], self.values_aux[1].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[2], self.values_aux[2].GetDouble(), self.fluid_solver.main_model_part.Nodes)

        #print("aquiiiiiiiiiiiiiiiiiiiiiiiiii") 
        #print("aquiiiiiiiiiiiiiiiiiiiiiiiiii") 
        #print("aquiiiiiiiiiiiiiiiiiiiiiiiiii") 
        #here we assign ACTIVATION_ENERGY, ARRHENIUS_COEFFICIENT and HEAT_OF_VAPORIZATION taken from FluidMaterial.json 
        #print(self.variables_aux[0])
        #print(self.self.values_aux[0])

        #print(self.variables_aux[1])
        #print(self.self.values_aux[1])
 
        #print(self.variables_aux[2])
        #print(self.self.values_aux[2])
        #hkjhkjhkjh
        #for node in self.fluid_solver.main_model_part.Nodes:
        #    rho = node.GetSolutionStepValue(PfemM.ACTIVATION_ENERGY)
        #    rho1 = node.GetSolutionStepValue(PfemM.ARRHENIUS_COEFFICIENT)
        #    rho2 = node.GetSolutionStepValue(PfemM.HEAT_OF_VAPORIZATION)


        #    print(rho) 
        #    print(rho1) 

        #    print(rho2) 
        #ssssssssssssssssss

         
        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X,0,self.gravity[0].GetDouble())
 
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y,0,self.gravity[1].GetDouble()) 
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Z,0,self.gravity[2].GetDouble()) 
 	
        

    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.thermal_solver.AddDofs()
        

    #def CalculateViscosity(self):
    #    import math
    #    for node in self.fluid_solver.main_model_part.Nodes:
    #        rho = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
    #        T = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
    #        #mu = AuxFunction(T)
    #        mu=self.new_table_aux.GetValue(T)
    #        node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY,0,mu/rho)


    def CalculateViscosityaux(self):
        import math
        for node in self.fluid_solver.main_model_part.Nodes:
            rho = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
            T = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            #mu=100000000000.0
            #if(node.GetSolutionStepValue(KratosMultiphysics.IS_INTERFACE)==0): 
            mu=self.new_table_aux.GetValue(T)
            #activationenergy=self.new_table_aux2.GetValue(T)
            #arrheniuscoeff=self.new_table_aux3.GetValue(T)
            node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY,0,mu/rho)
           
            #for node in self.fluid_solver.main_model_part.Nodes:
            #    rho = node.GetSolutionStepValue(PfemM.ACTIVATION_ENERGY)
            #    rho1 = node.GetSolutionStepValue(PfemM.ARRHENIUS_COEFFICIENT)
            #    rho2 = node.GetSolutionStepValue(PfemM.HEAT_OF_VAPORIZATION)


            #    print(rho) 
            #    print(rho1) 

            #    print(rho2) 
            #ssssssssssssssssss
            #node.SetSolutionStepValue(PfemM.ACTIVATION_ENERGY,0,activationenergy)
            #node.SetSolutionStepValue(PfemM.ARRHENIUS_COEFFICIENT,0,arrheniuscoeff)
            #node.SetSolutionStepValue(PfemM.ACTIVATION_ENERGY,0,0.0)
            #node.SetSolutionStepValue(PfemM.ARRHENIUS_COEFFICIENT,0,0.0)

    def cleaning_submodelparts(self):


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


    def filling_submodelparts(self):
        Parts_Parts_Auto1=self.fluid_solver.main_model_part.GetSubModelPart("FluidParts_Parts_Auto1")
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
        

        #self.main_model_part.CreateSubModelPart("fluid_computational_model_part")
        #fluid_computational_model_part= self.main_model_part.GetSubModelPart("fluid_computational_model_part")
        #fluid_computational_model_part.ProcessInfo = self.main_model_part.ProcessInfo


        #sssssssssssss

        if not self.thermal_solver.main_model_part.HasSubModelPart("thermal_computing_domain"):
            self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")
             
        #self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")


        
        thermal_computing_domain=self.thermal_solver.main_model_part.GetSubModelPart("thermal_computing_domain")
        thermal_computing_domain.Conditions.clear()	 
        thermal_computing_domain.Elements.clear()
        thermal_computing_domain.Nodes.clear()
        

        Parts_Parts_Auto1=self.thermal_solver.main_model_part.GetSubModelPart("FluidParts_Parts_Auto1")
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

        #neighbor_searcht = KratosMultiphysics.FindNodalNeighboursProcess(thermal_computing_domain)
        #neighbor_searcht.Execute()

        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        neighbor_searcht = KratosMultiphysics.FindGlobalNodalNeighboursProcess(kratos_comm,thermal_computing_domain)
        neighbor_searcht.Execute()

        #neighbor_elements_searcht = KratosMultiphysics.FindElementalNeighboursProcess(thermal_computing_domain, 3, 20)
        #neighbor_elements_searcht.Execute()

        neighbor_elements_searcht = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(kratos_comm,thermal_computing_domain)
        neighbor_elements_searcht.Execute()
        neighbor_condition_searcht = KratosMultiphysics.FindConditionsNeighboursProcess(thermal_computing_domain,3, 20)
        neighbor_condition_searcht.Execute()


       

    def ReMesh(self):
        

       
        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,0,0.15);
            node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,0,0.0014);
            node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,0,0.10);
            #node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,0,0.095);

        for node in (self.fluid_solver.main_model_part).Nodes:
            node.Set(KratosMultiphysics.TO_ERASE, False)

        #for node in (self.fluid_solver.main_model_part).Nodes:
        #    if(node.X>0.01299):
        #        node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0,1.0);
        #    if(node.Y<-0.043968):
        #        node.Set(KratosMultiphysics.TO_ERASE, True)
                #node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0,1.0);


        #self.Pfem2Utils.MarkNodesTouchingWall(self.fluid_solver.main_model_part, 3, 0.15)

        #self.node_erase_process.Execute()

        self.fluid_solver.main_model_part.Conditions.clear()

        self.fluid_solver.main_model_part.Elements.clear()
        
        (self.Mesher).ReGenerateMesh("LagrangianFluidVMS3D","ThermalFace3D3N", self.fluid_solver.main_model_part, self.node_erase_process, True, False, 1.4, 0.01)  #1.8

        #LagrangianFluidVMS3D
        #VMS3D 
        #neighbor_search = KratosMultiphysics.FindNodalNeighboursProcess(self.fluid_solver.main_model_part)
        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        neighbor_search = KratosMultiphysics.FindGlobalNodalNeighboursProcess(kratos_comm,self.fluid_solver.main_model_part)

        neighbor_search.Execute()
        #neighbor_elements_search = KratosMultiphysics.FindElementalNeighboursProcess(self.fluid_solver.main_model_part, 3, 20)
        neighbor_elements_search = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(kratos_comm,self.fluid_solver.main_model_part)


        #neighbor_elements_search = KratosMultiphysics.FindGlobalElementalNeighboursProcess(self.fluid_solver.main_model_part, 3, 20)
        neighbor_elements_search.Execute()
        neighbor_condition_search = KratosMultiphysics.FindConditionsNeighboursProcess(self.fluid_solver.main_model_part,3, 20)
        neighbor_condition_search.Execute()


        (self.PfemM_apply_bc_process).Execute();


       

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

        self.cleaning_submodelparts()

        #inverted=self.Streamline.CheckInvertElement(self.fluid_solver.main_model_part,self.domain_size )

        #if(inverted==True): 
        self.ReMesh()
         
        self.filling_submodelparts() 
        
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



        #for node in self.fluid_solver.main_model_part.Nodes:  
        #    node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 0.0) 
            #node.Free(KratosMultiphysics.VELOCITY_X)
            #node.Free(KratosMultiphysics.VELOCITY_Y)
            #node.Free(KratosMultiphysics.VELOCITY_Z)
            #node.Free(KratosMultiphysics.PRESSURE)
        #    if(node.Y<0.0001):
        #        node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 1.0) 
                #node.Fix(KratosMultiphysics.VELOCITY_X)
                #node.Fix(KratosMultiphysics.VELOCITY_Y)
                #node.Fix(KratosMultiphysics.VELOCITY_Z)
            #if(node.Y>0.4999):
            #    node.Fix(KratosMultiphysics.PRESSURE)
            #    node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0, 0.0) 

        
        fluid_is_converged = self.fluid_solver.SolveSolutionStep()
        self.Streamline.RungeKutta4ElementbasedSI(self.fluid_solver.main_model_part,100)
        
         
        for node in self.fluid_solver.main_model_part.Nodes:
            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY, velocity)
            #mesh_velocity = node.GetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY) 
            #print("velocities")
            #print(velocity)
            #print(mesh_velocity)
        
             

        x=self.values[1].GetDouble()
        y=self.values[2].GetDouble()
        z=self.values[3].GetDouble()
        radius=self.values[4].GetDouble()
        q=self.values[0].GetDouble()
        
        for node in self.fluid_solver.main_model_part.Nodes:
             node.SetSolutionStepValue(KratosMultiphysics.FACE_HEAT_FLUX,0,0.0);

        self.faceheatflux.FaceHeatFluxDistribution(self.fluid_solver.main_model_part, x, y, z, radius, q)
        
        self.HeatSource.Heat_Source(self.fluid_solver.main_model_part) #heat source for the thermal problem

        #for node in self.fluid_solver.main_model_part.Nodes:
        #     node.SetSolutionStepValue(KratosMultiphysics.FACE_HEAT_FLUX,0,0.0);

        #self.faceheatflux.FlameDistribution(self.fluid_solver.main_model_part, 8000.0)
        

        #for node in self.fluid_solver.main_model_part.Nodes:
        ##    if(node.Y<0.14587):
        #    if(node.GetSolutionStepValue(KratosMultiphysics.IS_FREE_SURFACE)==1): #if(node.Y>0.99999 or node.Y<0.00001):
        #        node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0,1000.0);
        #        node.Fix(KratosMultiphysics.TEMPERATURE)



        #for node in self.fluid_solver.main_model_part.Nodes:  
        #    aux=node.X * node.X + node.Z * node.Z 
        #    radious=aux**.5
        #    if (radious <=0.0017):
        #        node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,8920) 

        #    aux=node.X * node.X + node.Z * node.Z 
        #    radious=aux**.5
        #    if(node.X>0.0032931):
        #        XP=0.013-node.X
        #        ZP=0-node.Z
        #        r=(XP * XP + ZP * ZP)**.5
        #        if (r < 0.0016):
        #            node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,8920) 


        thermal_is_converged = self.thermal_solver.SolveSolutionStep()
        
        self.CalculateViscosityaux()	

        return (fluid_is_converged and thermal_is_converged)

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

    def Solve(self):
         
        self.InitializeSolutionStep()
        self.Predict()
        sssssssssss 
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()
