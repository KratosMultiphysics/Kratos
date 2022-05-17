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

    return PfemCoupledFluidThermalSolver(main_model_part, custom_settings)

class PfemCoupledFluidThermalSolver(PythonSolver):

    @classmethod
    def GetDefaultParameters(cls):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "ThermallyCoupled1",
            "domain_size" : 2,
            "echo_level": 0,
            "material_import_settings"    : {
                "materials_filename" : "file_name_to_be_defined.json"
            },
            "environment_settings" : {
                "gravity": [0, 0, 0],
                "ambient_temperature" : 0.15
            },
            "mesh_element_size"    : 0.0,
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
        print("PfemCoupledFluidThermalSolver(PythonSolver) in pfem2_solver.py has been constructed1")
        super(PfemCoupledFluidThermalSolver, self).__init__(model, custom_settings)
        print("PfemCoupledFluidThermalSolver(PythonSolver) in pfem2_solver.py has been constructed2")
        #self.settings["fluid_solver_settings"].AddEmptyValue("alpha")
        #self.settings["fluid_solver_settings"]["alpha"].SetDouble(0.0)
        self.settings["fluid_solver_settings"].AddEmptyValue("move_mesh_strategy")
        self.settings["fluid_solver_settings"]["move_mesh_strategy"].SetInt(2) 
        # To take into account
        #0 NS equations including convective terms
        #2 NS equations without convective terms
        print("PfemCoupledFluidThermalSolver(PythonSolver) in pfem2_solver.py has been constructed3")
         
        self.settings["fluid_solver_settings"].AddEmptyValue("reform_dofs_at_each_step")
        self.settings["fluid_solver_settings"]["reform_dofs_at_each_step"].SetBool(False)
        self.settings["thermal_solver_settings"].AddEmptyValue("reform_dofs_at_each_step")
        self.settings["thermal_solver_settings"]["reform_dofs_at_each_step"].SetBool(False)

        ## Get domain size
        self.domain_size = 2#self.settings["domain_size"].GetInt()

        from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"],"OpenMP")

        print("PfemCoupledFluidThermalSolver(PythonSolver) in pfem2_solver.py has been constructed4")

        self.readmeshSettings()

        self.readMaterialCharacterization()
        
        self.section_nodes = [] 
        
        self.outstring5 = "Drag_Lift"
        self.outputfile6 = open(self.outstring5, 'w')

        print("PfemCoupledFluidThermalSolver(PythonSolver) in pfem2_solver.py has been constructed5")

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


        self.new_table_aux = KratosMultiphysics.PiecewiseLinearTable()

        for i in range(table1["data"].size()):
            self.new_table_aux.AddRow(table1["data"][i][0].GetDouble(), table1["data"][i][1].GetDouble())


    def AddVariables(self):
        # Import the fluid and thermal solver variables. Then merge them to have them in both fluid and thermal solvers.
        print("AddVariables")
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
        print("InitializeSolutionStep")
        #self.step=1
        # section_nodes = [] 
        print(self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP])
        #if(self.step==self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]):
        if(1==2):    
            print("setting bcs")
            for node in self.fluid_solver.main_model_part.Nodes:
                if (node.GetSolutionStepValue(KratosMultiphysics.FLAG_VARIABLE)==1):
                    section_nodes.append(node)
                #for node in lagrangian_model_part.Nodes :
                node.Free(KratosMultiphysics.VELOCITY_X)
                node.Free(KratosMultiphysics.VELOCITY_Y)
                node.Free(KratosMultiphysics.VELOCITY_Z)
                node.Free(KratosMultiphysics.PRESSURE)
                node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X,0,0.0)
                node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y,0,0.0)
                node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Z,0,0.0)
        
                
                
                if(node.Y>5.49999):
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,1.0)
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    
                if(node.Y<-5.49999):
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,1,1.0)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,1,0.0)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,1,0.0)

                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,1.0)
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    
                if(node.X<-5.4999999):
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,1,1.0)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,1,0.0)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,1,0.0)

                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,1.0)
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    
                if(node.X>15.499999):
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,1,0.0)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,1,0.0)

                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    #node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0,0.0);
                    #node.Fix(KratosMultiphysics.PRESSURE)
                    
                if(node.GetSolutionStepValue(KratosMultiphysics.FLAG_VARIABLE)==1):
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,1,0.0)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,1,0.0)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,1,0.0)

                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    
                    #node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0,0.0);
                    #node.Fix(KratosMultiphysics.PRESSURE)

                vel_x=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0);    
                vel_y=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0);    
                vel_z=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0);    
                vel_old_x=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,1);    
                vel_old_y=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,1);    
                vel_old_z=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,1);    
                node.SetSolutionStepValue(PfemM.SCALARVELOCITY_X,0,vel_x);    
                node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0,vel_y);  
                node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0,vel_z);  
                node.SetSolutionStepValue(PfemM.SCALARVELOCITY_X,1,vel_old_x);    
                node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,1,vel_old_y);  
                node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,1,vel_old_z);  
        
        
            # self.section_nodes= section_nodes 

        self.fluid_solver.InitializeSolutionStep()

    def SetInitialConditions(self):
        print("SetInitialConditions")
        # self.step=1
        section_nodes = [] 
        for node in self.fluid_solver.main_model_part.Nodes:
            if (node.GetSolutionStepValue(KratosMultiphysics.FLAG_VARIABLE)==1):
                    section_nodes.append(node)
            #for node in lagrangian_model_part.Nodes :
            node.Free(KratosMultiphysics.VELOCITY_X)
            node.Free(KratosMultiphysics.VELOCITY_Y)
            node.Free(KratosMultiphysics.VELOCITY_Z)
            node.Free(KratosMultiphysics.PRESSURE)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X,0,0.0)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y,0,0.0)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Z,0,0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,1.0)
        
                
                
            if(node.Y>5.49999):
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,1.0)
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    
            if(node.Y<-5.49999):

                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,1.0)
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    
            if(node.X<-5.4999999):

                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,1.0)
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    
            if(node.X>15.499999):

                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    #node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0,0.0);
                    #node.Fix(KratosMultiphysics.PRESSURE)
                    
            if(node.GetSolutionStepValue(KratosMultiphysics.FLAG_VARIABLE)==1):

                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0,0.0)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    
                    #node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0,0.0);
                    #node.Fix(KratosMultiphysics.PRESSURE)

            vel_x=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0);    
            vel_y=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0);    
            vel_z=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0);    
            node.SetSolutionStepValue(PfemM.SCALARVELOCITY_X,0,vel_x);    
            node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0,vel_y);  
            node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0,vel_z);  

        self.section_nodes= section_nodes     
        
        
        

    def Predict(self):

        self.fluid_solver.Predict()

    def SolveSolutionStep(self):

        if self.settings["fluid_solver_settings"]["move_mesh_strategy"].GetInt()==2: 
        # add here the function to evaluate the convective terms  

         #KratosMultiphysics.FindGlobalNodalNeighboursProcess(self.fluid_solver.main_model_part).Execute()
         for node in self.fluid_solver.main_model_part.Nodes:
          vel_x=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0);  
          vel_y=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0);    
          vel_z=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0); 
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_X,0,vel_x);    
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0,vel_y);  
          node.SetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0,vel_z);  
         from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        #  linear_solver = linear_solver_factory.ConstructSolver(
        #      KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))
         linear_solver = linear_solver_factory.ConstructSolver(
             KratosMultiphysics.Parameters("""{"solver_type" : "amgcl", "max_iteration": 1000}"""))    

        # def test_bicgstab_iluk(self):
        # self._RunParametrized("""
        #     {
        #         "test_list" : [
        #             {

        #                 "solver_type" : "amgcl",
        #                 "smoother_type":"iluk",
        #                 "krylov_type": "bicgstab",
        #                 "coarsening_type": "aggregation",
        #                 "max_iteration": 100,
        #                 "provide_coordinates": false,
        #                 "gmres_krylov_space_dimension": 100,
        #                 "verbosity" : 0,
        #                 "tolerance": 1e-6,
        #                 "scaling": false,
        #                 "block_size": 1,
        #                 "use_block_matrices_if_possible" : true,
        #                 "coarse_enough" : 100
        #             }
        #         ]
        #     }
        #     """)
         levelset_convection_settings = KratosMultiphysics.Parameters("""{
             "levelset_variable_name" : "SCALARVELOCITY_X",
             "levelset_convection_variable_name" : "VELOCITY",
             "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
             "eulerian_error_compensation" : true,
             "element_type" : "levelset_convection_supg"
         }""")
         print("Before the first execute")
         KratosMultiphysics.LevelSetConvectionProcess2D(
             self.fluid_solver.main_model_part,
             linear_solver,
             levelset_convection_settings).Execute()
         print("After the first execute")    
         levelset_convection_settings = KratosMultiphysics.Parameters("""{
             "levelset_variable_name" : "SCALARVELOCITY_Y",
             "levelset_convection_variable_name" : "VELOCITY",
             "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
             "eulerian_error_compensation" : true,
             "element_type" : "levelset_convection_supg"
         }""")
         print("Before the second execute")
         KratosMultiphysics.LevelSetConvectionProcess2D(
             self.fluid_solver.main_model_part,
             linear_solver,
             levelset_convection_settings).Execute()    
         print("After the second execute")
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

            scalar_vel_x=node.GetSolutionStepValue(PfemM.SCALARVELOCITY_X,0);    
            scalar_vel_y=node.GetSolutionStepValue(PfemM.SCALARVELOCITY_Y,0);    
            scalar_vel_z=node.GetSolutionStepValue(PfemM.SCALARVELOCITY_Z,0);    
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,1,scalar_vel_x)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,1,scalar_vel_y)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,1,scalar_vel_z)        



        fluid_is_converged = self.fluid_solver.SolveSolutionStep()
        time=self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        print("time=", time)
        

        fx = 0.0
        fy = 0.0
        counter=0
        for node in self.section_nodes:
         fx -= node.GetSolutionStepValue(KratosMultiphysics.REACTION_X,0)
         fy -= node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y,0)
         counter+=1
         if (counter==432 or counter==431):
          self.outputfile6.write(str(time)+" "+"\n")
          print("at the end")
        print(counter)

        
        # self.outputfile6.write(str(time)+" "+ str(fx) +" "+  str(fy) +" "+ str(0.0) +"\n")

        return (fluid_is_converged) 

    
    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()

    def Solve(self):

        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()
