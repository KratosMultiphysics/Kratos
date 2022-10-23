from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import math
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
#import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff
import KratosMultiphysics.MeshingApplication as MeshApp
import KratosMultiphysics.PfemMeltingApplication as PfemM
import time as timer

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    #return PfemCoupledFluidThermalSolver(main_model_part, custom_settings) 
    return PFEM2Solver(main_model_part, custom_settings)

class PFEM2Solver(PythonSolver):

    @classmethod
    def GetDefaultParameters(cls):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "pfem2",
            "domain_size" : 2,
            "echo_level": 0,
            "material_import_settings"    : {
                "materials_filename" : "materials.json"
            },
            "environment_settings" : {
                "gravity": [0, 0, 0],
                "ambient_temperature" : 0.15
            },
            "mesh_element_size"    : 0.0,
            "fluid_solver_settings": {
                "solver_type": "navier_stokes_solver_vmsmonolithic",
                "splitting_strategy"      : "ReverseStrangSplitting",
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
        super(PFEM2Solver, self).__init__(model, custom_settings)
        #self.settings["fluid_solver_settings"].AddEmptyValue("alpha")
        #self.settings["fluid_solver_settings"]["alpha"].SetDouble(0.0)
        self.settings["fluid_solver_settings"].AddEmptyValue("move_mesh_strategy")
        self.move_mesh_strategy=self.settings["fluid_solver_settings"]["move_mesh_strategy"].GetInt()
        self.settings["fluid_solver_settings"].AddEmptyValue("splitting_strategy")
        self.splitting_strategy=self.settings["fluid_solver_settings"]["splitting_strategy"].GetString()
        if (self.move_mesh_strategy==0):
            print("NS equations including convective terms")
        if (self.move_mesh_strategy==2):
            print("NS equations without convective terms")  
            print("Splitting Strategy is: ", self.splitting_strategy)    

        
        # To take into account
        #0 NS equations including convective terms
        #2 NS equations without convective terms
         
        self.settings["fluid_solver_settings"].AddEmptyValue("reform_dofs_at_each_step")
        self.settings["fluid_solver_settings"]["reform_dofs_at_each_step"].SetBool(False)

        ## Get domain size
        self.domain_size = 2#self.settings["domain_size"].GetInt()
        
        ## Get time-step size
        self.timestep=self.settings["fluid_solver_settings"]["time_stepping"]["time_step"].GetDouble()
        print(self.timestep)
        

        from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"],"OpenMP")

        self.readmeshSettings()

        self.readMaterialCharacterization()
        
        self.section_nodes = [] 
        
        self.outstring5 = "Drag_Lift_0"
        self.outputfile6 = open(self.outstring5, 'w')

        self.outstring7 = "Drag_Lift_1"
        self.outputfile8 = open(self.outstring7, 'w')



    def readmeshSettings(self):

        with open("ProjectParameters.json",'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        self.mesh_element_size=project_parameters["problem_data"]["mesh_element_size"]

    def readMaterialCharacterization(self):

        materials_filename = self.settings["thermal_solver_settings"]["material_import_settings"]["materials_filename"].GetString()

        with open(materials_filename, 'r') as parameter_file:
                materials = KratosMultiphysics.Parameters(parameter_file.read())

        mat = materials["properties"][0]["Material"]

        self.variables_aux = []
        self.values_aux = []
        for key, value in mat["Variables"].items():
            var = KratosMultiphysics.KratosGlobals.GetVariable(key)
            self.variables_aux.append(var)
            self.values_aux.append(value)

        #The part below reads the temperature dependent viscosity
        table1=mat["Tables"]["Table1"]

        input_var = KratosMultiphysics.KratosGlobals.GetVariable(table1["input_variable"].GetString())
        output_var = KratosMultiphysics.KratosGlobals.GetVariable(table1["output_variable"].GetString())

        read_materials_utility = KratosMultiphysics.ReadMaterialsUtility(self.model)

        read_materials_utility.AssignVariablesToProperty(mat, self.fluid_solver.main_model_part.GetProperties()[0])
        #read_materials_utility.AssignVariablesToProperty(mat, self.fluid_solver.main_model_part.GetProperties()[1])

        '''self.fluid_solver.main_model_part.GetProperties()[0][KratosMultiphysics.DYNAMIC_VISCOSITY]=self.values_aux[6].GetDouble()
        self.fluid_solver.main_model_part.GetProperties()[0][KratosMultiphysics.DENSITY]=self.values_aux[5].GetDouble()'''    




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
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
       

        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ACTIVATION_ENERGY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ARRHENIUS_COEFFICIENT)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HEAT_OF_VAPORIZATION)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ARRHENIUS_VALUE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.SCALARVELOCITY_X)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.SCALARVELOCITY_Y)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.SCALARVELOCITY_Z)



    def ImportModelPart(self):
        # Call the fluid solver to import the model part from the mdpa

        self.fluid_solver.ImportModelPart()



    def PrepareModelPart(self):


        self.fluid_solver.PrepareModelPart()

        #self.cleaning_submodelparts()
        self.assign_nodally_properties();
        #self.ReMesh()


    def assign_nodally_properties(self):
        #here we assign ACTIVATION_ENERGY, ARRHENIUS_COEFFICIENT and HEAT_OF_VAPORIZATION taken from FluidMaterial.json
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[0], self.values_aux[0].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[1], self.values_aux[1].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[2], self.values_aux[2].GetDouble(), self.fluid_solver.main_model_part.Nodes)




    def AddDofs(self):
        self.fluid_solver.AddDofs()

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):

        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid)

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.SetInitialConditions()
        self.CalculateNodalArea()

    def Clear(self):

        (self.fluid_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):


        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.fluid_solver.AdvanceInTime(current_time)
        return new_time

    def InitializeSolutionStep(self):
        #self.step=1
        #section_nodes = [] 
        #print(self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP])
        #if(self.step==self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]):

        self.fluid_solver.InitializeSolutionStep()

    def CalculateNodalArea(self):    
        nodal_area_process = KratosMultiphysics.CalculateNodalAreaProcess(self.fluid_solver.main_model_part,2)
        nodal_area_process.Execute()
    def SetInitialConditions(self):
        self.step=0
        section_nodes = [] 
        mu=0.001
        rho=1.0
        time=0.0
        for node in self.fluid_solver.main_model_part.Nodes:
            node.Free(KratosMultiphysics.VELOCITY_X)
            node.Free(KratosMultiphysics.VELOCITY_Y)
            node.Free(KratosMultiphysics.VELOCITY_Z)
            node.Free(KratosMultiphysics.PRESSURE)
            vel_x=-1.0*math.sin(node.X)*math.cos(node.Y)*math.exp(-2.0*(mu/rho)*time)
            vel_y=math.cos(node.X)*math.sin(node.Y)*math.exp(-2.0*(mu/rho)*time)    
            vel_z=0.0  
            pressure=(rho/4.0)*(math.cos(2.0*node.X)+math.cos(2.0*node.Y))*math.exp(-4.0*(mu/rho)*time) 
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,vel_x)   
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,vel_y) 
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,vel_z) 
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0,pressure)
            if(node.X>6.27 or node.Y>6.27 or node.X<0.001 or node.Y<0.001):
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    node.Fix(KratosMultiphysics.PRESSURE)
                    vel_x=-1.0*math.sin(node.X)*math.cos(node.Y)*math.exp(-2.0*(mu/rho)*time)
                    vel_y=math.cos(node.X)*math.sin(node.Y)*math.exp(-2.0*(mu/rho)*time)
                    vel_z=0.0
                    pressure=(rho/4.0)*(math.cos(2.0*node.X)+math.cos(2.0*node.Y))*math.exp(-4.0*(mu/rho)*time)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,vel_x)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,vel_y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,vel_z)
                    node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0,pressure)
                    
                    
                    

            

        self.section_nodes= section_nodes     
        
        
        

    def Predict(self):

        self.fluid_solver.Predict()

    def SolveSolutionStep(self):
        if self.move_mesh_strategy==2: 
         if self.splitting_strategy=="StrangSplitting":   
         # add here the function to evaluate the convective terms  

          self.ReduceTimeStep()
          fluid_is_converged = self.fluid_solver.SolveSolutionStep()
          self.IncreaseTimeStep()
          
          self.CopyCurrentVelocityToCurrentScalarValues()
          self.ConvectScalarX()
          self.ConvectScalarY()  
          self.ApplyBCsToCurrentScalarValues()
          self.CopyCurrentScalarValuesToOldVelocityValues()

          self.ReduceTimeStep()
          fluid_is_converged = self.fluid_solver.SolveSolutionStep()
          self.IncreaseTimeStep()

          self.CalculateTheError()

          return (fluid_is_converged) 


         elif self.splitting_strategy=="ReverseStrangSplitting": 
         # add here the function to evaluate the convective terms  

          self.ReduceTimeStep()
          self.CopyCurrentVelocityToCurrentScalarValues() 
          self.ConvectScalarX()
          self.ConvectScalarY()  
          self.ApplyBCsToCurrentScalarValues()
          self.CopyCurrentScalarValuesToOldVelocityValues()

          self.IncreaseTimeStep()  
          fluid_is_converged = self.fluid_solver.SolveSolutionStep()

          self.ReduceTimeStep()
          self.CopyCurrentVelocityToCurrentScalarValues() 
          self.ConvectScalarX()
          self.ConvectScalarY()    
          self.ApplyBCsToCurrentScalarValues()
          self.CopyCurrentScalarValuesToCurrentVelocityValues()  
          self.IncreaseTimeStep()

          self.CalculateTheError()

          return (fluid_is_converged)

         elif self.splitting_strategy=="FirstOrderSplitting": 
         # add here the function to evaluate the convective terms  

          self.CopyCurrentVelocityToCurrentScalarValues() 
          self.ConvectScalarX()
          self.ConvectScalarY()  
          self.ApplyBCsToCurrentScalarValues()
          self.CopyCurrentScalarValuesToOldVelocityValues()

          fluid_is_converged = self.fluid_solver.SolveSolutionStep()

          self.CalculateTheError()

          return (fluid_is_converged)
         else:
          raise Exception("Splitting Method has not been defined")
  
        elif (self.move_mesh_strategy==0):  
         self.ApplyBCs()  
         fluid_is_converged = self.fluid_solver.SolveSolutionStep()

         self.CalculateTheError()

         return (fluid_is_converged)
        else:
          raise Exception("move_mesh_strategy is neither 0 nor 2") 
         
        


    
    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()

    def ReduceTimeStep(self):
        newdt=0.5*self.timestep
        self.fluid_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, newdt)
    def IncreaseTimeStep(self):
        olddt=self.timestep
        self.fluid_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, olddt)    
    def CopyCurrentVelocityToCurrentScalarValues(self):
        for node in self.fluid_solver.main_model_part.Nodes:
         vel_x=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0);  
         vel_y=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0);    
         vel_z=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0); 
         node.SetSolutionStepValue(PfemM.SCALARVELOCITY_X,0,vel_x);    
         node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0,vel_y);  
         node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0,vel_z);   
    def ConvectScalarX(self):
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(self.fluid_solver.main_model_part).Execute()
        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(KratosMultiphysics.Parameters("""{"solver_type" : "amgcl", "max_iteration": 1000}"""))  
        #linear_solver = linear_solver_factory.ConstructSolver(KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))
        levelset_convection_settings = KratosMultiphysics.Parameters("""{
             "levelset_variable_name" : "SCALARVELOCITY_X",
             "levelset_convection_variable_name" : "VELOCITY",
             "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
             "eulerian_error_compensation" : true,
             "element_type" : "levelset_convection_supg"}""")
        KratosMultiphysics.LevelSetConvectionProcess2D(self.fluid_solver.main_model_part,linear_solver,levelset_convection_settings).Execute()   
    def ConvectScalarY(self):
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(self.fluid_solver.main_model_part).Execute()
        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(KratosMultiphysics.Parameters("""{"solver_type" : "amgcl", "max_iteration": 1000}"""))  
        #linear_solver = linear_solver_factory.ConstructSolver(KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))
        levelset_convection_settings = KratosMultiphysics.Parameters("""{
             "levelset_variable_name" : "SCALARVELOCITY_Y",
             "levelset_convection_variable_name" : "VELOCITY",
             "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
             "eulerian_error_compensation" : true,
             "element_type" : "levelset_convection_supg"}""")
        KratosMultiphysics.LevelSetConvectionProcess2D(self.fluid_solver.main_model_part,linear_solver,levelset_convection_settings).Execute() 
    def ApplyBCs(self):
        time=self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.TIME]   
        mu=0.001
        rho=1.0
        for node in self.fluid_solver.main_model_part.Nodes:
            node.Free(KratosMultiphysics.VELOCITY_X)
            node.Free(KratosMultiphysics.VELOCITY_Y)
            node.Free(KratosMultiphysics.VELOCITY_Z)
            node.Free(KratosMultiphysics.PRESSURE)
            if(node.X>6.27 or node.Y>6.27 or node.X<0.001 or node.Y<0.001):
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    node.Fix(KratosMultiphysics.PRESSURE)
                    vel_x=-1.0*math.sin(node.X)*math.cos(node.Y)*math.exp(-2.0*(mu/rho)*time)
                    vel_y=math.cos(node.X)*math.sin(node.Y)*math.exp(-2.0*(mu/rho)*time)
                    vel_z=0.0
                    pressure=(rho/4.0)*(math.cos(2.0*node.X)+math.cos(2.0*node.Y))*math.exp(-4.0*(mu/rho)*time)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,vel_x)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,vel_y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,vel_z)
                    node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0,pressure) 
    def ApplyBCsToCurrentScalarValues(self):
        for node in self.fluid_solver.main_model_part.Nodes:              
         if(node.Y>5.49999):
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_X,0,1.0)
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0,0.0)
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0,0.0)
                    
         if(node.Y<-5.49999):
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_X,0,1.0)
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0,0.0)
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0,0.0)
                    
         if(node.X<-5.4999999):
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_X,0,1.0)
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0,0.0)
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0,0.0)
                    
         if(node.X>15.499999):
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0,0.0)
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0,0.0)
                    
         if(node.GetSolutionStepValue(KratosMultiphysics.FLAG_VARIABLE)==1):
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_X,0,0.0)
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0,0.0)
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0,0.0)
    def CopyCurrentScalarValuesToOldVelocityValues(self):
        for node in self.fluid_solver.main_model_part.Nodes:              
         scalar_vel_x=node.GetSolutionStepValue(PfemM.SCALARVELOCITY_X,0);    
         scalar_vel_y=node.GetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0);    
         scalar_vel_z=node.GetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0);    
         node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,1,scalar_vel_x)
         node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,1,scalar_vel_y)
         node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,1,scalar_vel_z) 
    def CopyCurrentScalarValuesToCurrentVelocityValues(self):
        for node in self.fluid_solver.main_model_part.Nodes:              
         scalar_vel_x=node.GetSolutionStepValue(PfemM.SCALARVELOCITY_X,0);    
         scalar_vel_y=node.GetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0);    
         scalar_vel_z=node.GetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0);    
         node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,scalar_vel_x)
         node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,scalar_vel_y)
         node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,scalar_vel_z)   
    def CalculateTheError(self):
        time=self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.TIME]    
        mu=0.001
        rho=1.0
        totalerrorvelocitywithnodalarea=0.0
        totalerrorpressurewithnodalarea=0.0
        totalerrorvelocity=0.0
        totalerrorpressure=0.0
        nodalareasum=0.0
        for node in self.fluid_solver.main_model_part.Nodes:
         vel_x=-1.0*math.sin(node.X)*math.cos(node.Y)*math.exp(-2.0*(mu/rho)*time)  
         vel_y=math.cos(node.X)*math.sin(node.Y)*math.exp(-2.0*(mu/rho)*time)
         vel_z=0.0
         pressure=(rho/4.0)*(math.cos(2.0*node.X)+math.cos(2.0*node.Y))*math.exp(-4.0*(mu/rho)*time) 
         numerical_vel_x=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
         numerical_vel_y=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
         numerical_pressure=node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
         nodalarea=node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
         nodalareasum+=nodalarea
         errorvel_x = vel_x-numerical_vel_x
         errorvel_y = vel_y-numerical_vel_y
         errorpressure = pressure-numerical_pressure
         totalerrorvelocitywithnodalarea +=(errorvel_x*errorvel_x+errorvel_y*errorvel_y)*nodalarea
         totalerrorpressurewithnodalarea +=(errorpressure*errorpressure)*nodalarea
         totalerrorvelocity +=(errorvel_x*errorvel_x+errorvel_y*errorvel_y)
         totalerrorpressure +=(errorpressure*errorpressure)
         node.SetSolutionStepValue(PfemM.SCALARVELOCITY_X,0,math.sqrt(errorvel_x*errorvel_x))
         node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0,math.sqrt(errorvel_y*errorvel_y))
         node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0,math.sqrt(errorpressure*errorpressure))
        totalerrorvelocitywithnodalarea=totalerrorvelocitywithnodalarea/nodalareasum
        totalerrorpressurewithnodalarea=totalerrorpressurewithnodalarea/nodalareasum
        totalerrorvelocitywithnodalarea  = math.sqrt(totalerrorvelocitywithnodalarea)  
        totalerrorpressurewithnodalarea  = math.sqrt(totalerrorpressurewithnodalarea)
        totalerrorvelocity  = math.sqrt(totalerrorvelocity)  
        totalerrorpressure  = math.sqrt(totalerrorpressure)  
        print("dt=", self.timestep)
        print("totalerrorvelocitywithnodalarea=", totalerrorvelocitywithnodalarea) 
        print("totalerrorpressurewithnodalarea=", totalerrorpressurewithnodalarea)
        print("totalerrorvelocity=", totalerrorvelocity) 
        print("totalerrorpressure=", totalerrorpressure)   

    def Solve(self):

        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()
