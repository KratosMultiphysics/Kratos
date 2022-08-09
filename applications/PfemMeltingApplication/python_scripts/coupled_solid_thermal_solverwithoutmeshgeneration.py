from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import math
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff
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
            "domain_size" : -1,
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

        super(PfemCoupledFluidThermalSolver, self).__init__(model, custom_settings)

        self.settings["fluid_solver_settings"].AddEmptyValue("alpha")
        self.settings["fluid_solver_settings"]["alpha"].SetDouble(0.0)
        self.settings["fluid_solver_settings"].AddEmptyValue("move_mesh_strategy")
        self.settings["fluid_solver_settings"]["move_mesh_strategy"].SetInt(2)
        self.settings["fluid_solver_settings"].AddEmptyValue("reform_dofs_at_each_step")
        self.settings["fluid_solver_settings"]["reform_dofs_at_each_step"].SetBool(True)
        self.settings["thermal_solver_settings"].AddEmptyValue("reform_dofs_at_each_step")
        self.settings["thermal_solver_settings"]["reform_dofs_at_each_step"].SetBool(True)

        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"],"OpenMP")



        from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion
        self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model,self.settings["thermal_solver_settings"],"OpenMP")



        self.readmeshSettings()

        self.readenvironmentSettings()

        self.readMaterialCharacterization()

        self.Mesher = MeshApp.TetGenPfemModeler()

        self.modeler = KratosMultiphysics.ConnectivityPreserveModeler()

        self.PfemM_apply_bc_process = PfemM.PfemMeltingApplyBCProcess(self.fluid_solver.main_model_part);

        self.node_erase_process = KratosMultiphysics.NodeEraseProcess(self.fluid_solver.main_model_part);

        self.Streamline = PfemM.Streamline()

        self.faceheatflux = PfemM.FaceHeatFlux()

        self.HeatSource = PfemM.HeatSource()

        self.outstring3 = "Volumen_withoutmesh"
        self.outputfile4 = open(self.outstring3, 'w')

        self.outstring5 = "computational_times"
        self.outputfile6 = open(self.outstring5, 'w')


        self.streamlineintegration=0.0
        self.fillingsubmodelparts=0.0
        self.initializeSolutionStep=0.0
        self.problemsolution=0.0
        self.meshingprocedure=0.0

        self.hypoelastic_solid_stress_tensor_calculate_process=PfemM.HypoelasticStressCalculateProcess(self.fluid_solver.main_model_part, self.domain_size)


    def readmeshSettings(self):

        with open("ProjectParameters.json",'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        self.mesh_element_size=project_parameters["problem_data"]["mesh_element_size"]

    def readenvironmentSettings(self):

        with open("ProjectParameters.json",'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())


        self.gravity = []
        self.ambient_temperature=project_parameters["problem_data"]["environment_settings"]["ambient_temperature"]
        self.gravity=project_parameters["problem_data"]["environment_settings"]["gravity"]

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
        self.fluid_solver.main_model_part.GetProperties()[0][KratosMultiphysics.DENSITY]=self.values_aux[5].GetDouble()
        self.fluid_solver.main_model_part.GetProperties()[0][KratosMultiphysics.EMISSIVITY]=self.values_aux[7].GetDouble()
        self.fluid_solver.main_model_part.GetProperties()[0][KratosMultiphysics.AMBIENT_TEMPERATURE]=self.values_aux[1].GetDouble() #298.0
        self.fluid_solver.main_model_part.GetProperties()[0][KratosMultiphysics.CONVECTION_COEFFICIENT]=self.values_aux[4].GetDouble()'''


        #self.fluid_solver.main_model_part.GetProperties()[0][KratosMultiphysics.AMBIENT_TEMPERATURE]=self.ambient_temperature.GetDouble()
        self.fluid_solver.main_model_part.GetProperties()[1][KratosMultiphysics.AMBIENT_TEMPERATURE]=self.ambient_temperature.GetDouble()

        self.new_table_aux = KratosMultiphysics.PiecewiseLinearTable()

        for i in range(table1["data"].size()):
            self.new_table_aux.AddRow(table1["data"][i][0].GetDouble(), table1["data"][i][1].GetDouble())


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




        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ACTIVATION_ENERGY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ARRHENIUS_COEFFICIENT)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HEAT_OF_VAPORIZATION)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.ARRHENIUS_VALUE)

        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.DELTA_SIGMA_XX)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.DELTA_SIGMA_XY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.DELTA_SIGMA_XZ)

        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.DELTA_SIGMA_YX)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.DELTA_SIGMA_YY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.DELTA_SIGMA_YZ)

        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.DELTA_SIGMA_ZX)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.DELTA_SIGMA_ZY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.DELTA_SIGMA_ZZ)




        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HISTORICAL_SIGMA_XX)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HISTORICAL_SIGMA_XY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HISTORICAL_SIGMA_XZ)

        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HISTORICAL_SIGMA_YX)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HISTORICAL_SIGMA_YY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HISTORICAL_SIGMA_YZ)

        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HISTORICAL_SIGMA_ZX)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HISTORICAL_SIGMA_ZY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(PfemM.HISTORICAL_SIGMA_ZZ)

        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)



        self.thermal_solver.AddVariables()

        KratosMultiphysics.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part)

    def ImportModelPart(self):
        # Call the fluid solver to import the model part from the mdpa

        self.fluid_solver.ImportModelPart()


        # Save the convection diffusion settings
        convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)

        # Here the fluid model part is cloned to be thermal model part so that the nodes are shared
        if self.domain_size == 2:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "EulerianConvDiffLumped2D",
                                      "ThermalFace2D2N")
        else:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "EulerianConvDiffLumped3D",
                                      "ThermalFace3D3N")


        # Set the saved convection diffusion settings to the new thermal model part
        self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)


    def PrepareModelPart(self):


        self.fluid_solver.PrepareModelPart()

        self.thermal_solver.PrepareModelPart()


        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 0.0)
        parametersf=self.settings["fluid_solver_settings"]

        self.cleaning_submodelparts()
        self.assign_nodally_properties();
        self.ReMesh()


    def assign_nodally_properties(self):
        #here we assign ACTIVATION_ENERGY, ARRHENIUS_COEFFICIENT and HEAT_OF_VAPORIZATION taken from FluidMaterial.json
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[0], self.values_aux[0].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[1], self.values_aux[1].GetDouble(), self.fluid_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(self.variables_aux[2], self.values_aux[2].GetDouble(), self.fluid_solver.main_model_part.Nodes)




        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X,0,self.gravity[0].GetDouble())
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y,0,self.gravity[1].GetDouble())
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Z,0,self.gravity[2].GetDouble())
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0,self.ambient_temperature.GetDouble())




    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.thermal_solver.AddDofs()



    def CalculateViscosityaux(self):
        import math
        for node in self.fluid_solver.main_model_part.Nodes:
            rho = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
            T = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            mu=self.new_table_aux.GetValue(T)
            node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY,0,mu/rho)


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


    def new(self):

        #if not self.fluid_solver.main_model_part.HasSubModelPart("fluid_computational_model_part"):
        #    ggggg
        #fluid_computational_model_part=self.fluid_solver.main_model_part.GetSubModelPart("fluid_computational_model_part")
        #fluid_computational_model_part.Conditions.clear()
        #fluid_computational_model_part.Elements.clear()
        #fluid_computational_model_part.Nodes.clear()

        ##newww
        for node in (self.fluid_solver.main_model_part).Nodes:
            node.Set(KratosMultiphysics.TO_ERASE, False)


        self.fluid_solver.main_model_part.Conditions.clear()


        neighbor_search = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.fluid_solver.main_model_part)
        neighbor_search.Execute()

        #(self.PfemM_apply_bc_process).Execute();

        if not self.fluid_solver.main_model_part.HasSubModelPart("SkinDEMModelPart"):
            print("nada")
        else:
            SkinDEMModelPart=self.fluid_solver.main_model_part.GetSubModelPart("SkinDEMModelPart")
            SkinDEMModelPart.Conditions.clear()
            SkinDEMModelPart.Elements.clear()
            SkinDEMModelPart.Nodes.clear()

        skin_params = KratosMultiphysics.Parameters("""
        {
            "name_auxiliar_model_part" : "SkinDEMModelPart",
            "name_auxiliar_condition"  : "Condition",
            "echo_level"               : 2
        }""")


        skin_detection_process = KratosMultiphysics.SkinDetectionProcess3D(self.fluid_solver.main_model_part, skin_params)

        skin_detection_process.Execute()


        neighbor_condition_search = KratosMultiphysics.FindConditionsNeighboursProcess(self.fluid_solver.main_model_part,3, 20)
        neighbor_condition_search.Execute()


        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY,0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.IS_FREE_SURFACE,0, 0.0)


        for node in self.fluid_solver.main_model_part.GetSubModelPart("SkinDEMModelPart").Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY,0, 1.0)



        (self.PfemM_apply_bc_process).Execute();

        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.fluid_solver.main_model_part.Nodes)
        normal_calculation_utils = KratosMultiphysics.NormalCalculationUtils()


        normal_calculation_utils.CalculateOnSimplex(self.fluid_solver.main_model_part.Conditions, self.domain_size)
        for node in self.fluid_solver.main_model_part.Nodes:
            if(node.GetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY)== 1.0):
                solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
                solution_normal /= math.sqrt(solution_normal[0]**2+solution_normal[1]**2+solution_normal[2]**2)
                node.SetSolutionStepValue(KratosMultiphysics.NORMAL, solution_normal)



        ###end new


        if not self.fluid_solver.main_model_part.HasSubModelPart("solid_model_part"):
            solid_model_part= self.fluid_solver.main_model_part.CreateSubModelPart("solid_model_part")


        if not self.fluid_solver.main_model_part.HasSubModelPart("fluid_model_part"):
            fluid_model_part= self.fluid_solver.main_model_part.CreateSubModelPart("fluid_model_part")


        solid_model_part= self.fluid_solver.main_model_part.GetSubModelPart("solid_model_part")
        solid_model_part.Elements.clear()

        fluid_model_part= self.fluid_solver.main_model_part.GetSubModelPart("fluid_model_part")
        fluid_model_part.Elements.clear()


                #for node in self.fluid_solver.main_model_part.Nodes:
            #rho = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
            #T = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)


        self.step=1
        if(self.step==self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]):

            for node in self.fluid_solver.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.IS_INTERFACE,0, 1.0)
            for elem in self.fluid_solver.main_model_part.Elements:
                solid_model_part.AddElement(elem,0)


        else:
            for node in self.fluid_solver.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.IS_INTERFACE,0, 1.0)

            for elem in self.fluid_solver.main_model_part.Elements:
                nnodes=0
                for node in elem.GetNodes():
                    T = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                    if(T >450.0):
                    #if (node.Z>0.0028 ): #or node.X<0.000):
                        node.SetSolutionStepValue(KratosMultiphysics.IS_INTERFACE,0, 0.0)
                        nnodes=nnodes+1;
                if(nnodes>1):
                    fluid_model_part.AddElement(elem,0)
                else:
                    solid_model_part.AddElement(elem,0)

        print("················")
        print("fluid model part")
        print(fluid_model_part)

        new_elem_name= "LagrangianFluidVMS3D" #"LagrangianFluidVMS3D"#"HYPO3D"

        new_cond_name= "ThermalFace3D3N"
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)
        #print(self.settings["element_replace_settings"])

        KratosMultiphysics.ReplaceElementsAndConditionsProcess(fluid_model_part, self.settings["element_replace_settings"]).Execute()

        new_elem_name= "HYPO3D" #"LagrangianFluidVMS3D"#"HYPO3D"

        new_cond_name= "ThermalFace3D3N"
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)
        #print(self.settings["element_replace_settings"])


        KratosMultiphysics.ReplaceElementsAndConditionsProcess(solid_model_part, self.settings["element_replace_settings"]).Execute()

        fluid_computational_model_part=self.fluid_solver.main_model_part.GetSubModelPart("fluid_computational_model_part")
        #fluid_computational_model_part.Conditions.clear()
        fluid_computational_model_part.Elements.clear()
        #fluid_computational_model_part.Nodes.clear()

        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(fluid_computational_model_part, fluid_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.ELEMENTS)
        transfer_process.Execute()

        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(fluid_computational_model_part, solid_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.ELEMENTS)
        transfer_process.Execute()


        #transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(fluid_computational_model_part, self.fluid_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODES)
        #transfer_process.Execute()


        fluid_computational_model_part.ProcessInfo = self.fluid_solver.main_model_part.ProcessInfo

        fluid_computational_model_part.Properties  = self.fluid_solver.main_model_part.Properties



        self.thermal_solver.main_model_part.Conditions.clear()
        self.thermal_solver.main_model_part.Elements.clear()
        self.thermal_solver.main_model_part.Nodes.clear()


        if not self.thermal_solver.main_model_part.HasSubModelPart("thermal_computing_domain"):
            self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")


        thermal_computing_domain=self.thermal_solver.main_model_part.GetSubModelPart("thermal_computing_domain")

        thermal_computing_domain.Conditions.clear()
        thermal_computing_domain.Elements.clear()
        thermal_computing_domain.Nodes.clear()

        if self.domain_size == 2:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, thermal_computing_domain, "EulerianConvDiffLumped2D", "ThermalFace2D2N")
        else:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, thermal_computing_domain,"EulerianConvDiffLumped3D","ThermalFace3D3N")


        #self.thermal_solver.main_model_part.Conditions.clear()


        thermal_computing_domain.ProcessInfo = self.fluid_solver.main_model_part.ProcessInfo
        thermal_computing_domain.Properties  = self.fluid_solver.main_model_part.Properties



        #transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(thermal_computing_domain, self.fluid_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS)
        #transfer_process.Execute()


        neighbor_search = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(thermal_computing_domain)
        neighbor_search.Execute()

        neighbor_condition_search = KratosMultiphysics.FindConditionsNeighboursProcess(thermal_computing_domain,3, 20)
        neighbor_condition_search.Execute()





    def ReMesh(self):



        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,0,self.mesh_element_size.GetDouble());

        for node in (self.fluid_solver.main_model_part).Nodes:
            node.Set(KratosMultiphysics.TO_ERASE, False)


        #self.fluid_solver.main_model_part.Conditions.clear()

        #self.fluid_solver.main_model_part.Elements.clear()

        (self.Mesher).ReGenerateMesh("LagrangianFluidVMS3D","ThermalFace3D3N", self.fluid_solver.main_model_part, self.node_erase_process, True, False, 1.2, 0.01)  #1.8

        #LagrangianFluidVMS3D
        #VMS3D
        neighbor_search = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.fluid_solver.main_model_part)
        neighbor_search.Execute()


        neighbor_condition_search = KratosMultiphysics.FindConditionsNeighboursProcess(self.fluid_solver.main_model_part,3, 20)
        neighbor_condition_search.Execute()


        (self.PfemM_apply_bc_process).Execute();




        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.fluid_solver.main_model_part.Nodes)
        normal_calculation_utils = KratosMultiphysics.NormalCalculationUtils()


        normal_calculation_utils.CalculateOnSimplex(self.fluid_solver.main_model_part.Conditions, self.domain_size)
        for node in self.fluid_solver.main_model_part.Nodes:
            if(node.GetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY)== 1.0):
                solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
                solution_normal /= math.sqrt(solution_normal[0]**2+solution_normal[1]**2+solution_normal[2]**2)
                node.SetSolutionStepValue(KratosMultiphysics.NORMAL, solution_normal)





        pass



    def CalculateNorm(array_3d_value):
        return math.sqrt(array_3d_value[0]**2+array_3d_value[1]**2+array_3d_value[2]**2)


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

        self.step=1
        if(self.step==self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]):
            for node in self.fluid_solver.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 0.0) #NODES NOT
                #if(node.IsFixed(KratosMultiphysics.VELOCITY_X)==True):
                #    node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 1.0) #NODES NOT
                #    node.Fix(KratosMultiphysics.TEMPERATURE)
                node.Free(KratosMultiphysics.VELOCITY_X)
                node.Free(KratosMultiphysics.VELOCITY_Y)
                node.Free(KratosMultiphysics.VELOCITY_Z)
                #if(node.X<-0.02499):
                if(node.Z<0.00000098):
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 1.0) #NODES NOT








        self.Streamline.RungeKutta4ElementbasedSI(self.fluid_solver.main_model_part,100)
        #self.Streamline.MovingParticles(self.fluid_solver.main_model_part,100)



        self.inverted=False

        print("ANTESSSSSSSS")
        print(self.fluid_solver.main_model_part)
        self.inverted=self.Streamline.CheckInvertElement(self.fluid_solver.main_model_part,self.domain_size,self.mesh_element_size.GetDouble())
        print("DESPUESSSSS")
        print(self.fluid_solver.main_model_part)


        #self.ReMesh()
        #self.new()

        self.AuxReMesh()


        if(self.step==self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]):

            self.filling_submodelparts()
            self.modelpart()
        else:


            self.filling_submodelparts_aux()
            self.modelpart()


        self.fluid_solver.InitializeSolutionStep()

        self.thermal_solver.InitializeSolutionStep()


    def modelpart(self):
        if not self.fluid_solver.main_model_part.HasSubModelPart("solid_model_part"):
            solid_model_part= self.fluid_solver.main_model_part.CreateSubModelPart("solid_model_part")


        if not self.fluid_solver.main_model_part.HasSubModelPart("fluid_model_part"):
            fluid_model_part= self.fluid_solver.main_model_part.CreateSubModelPart("fluid_model_part")


        solid_model_part= self.fluid_solver.main_model_part.GetSubModelPart("solid_model_part")
        solid_model_part.Elements.clear()

        fluid_model_part= self.fluid_solver.main_model_part.GetSubModelPart("fluid_model_part")
        fluid_model_part.Elements.clear()


        self.step=1
        if(self.step==self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]):

            for node in self.fluid_solver.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.IS_INTERFACE,0, 1.0)
            for elem in self.fluid_solver.main_model_part.Elements:
                solid_model_part.AddElement(elem,0)


        else:
            for node in self.fluid_solver.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.IS_INTERFACE,0, 1.0)

            for elem in self.fluid_solver.main_model_part.Elements:
                nnodes=0
                for node in elem.GetNodes():
                    T = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                    if(T >450.0):
                    #if (node.Z>0.0028 ): #or node.X<0.000):
                        node.SetSolutionStepValue(KratosMultiphysics.IS_INTERFACE,0, 0.0)
                        nnodes=nnodes+1;
                if(nnodes>1):
                    fluid_model_part.AddElement(elem,0)
                else:
                    solid_model_part.AddElement(elem,0)

        print("················")
        print("fluid model part")
        print(fluid_model_part)

        new_elem_name= "LagrangianFluidVMS3D" #"LagrangianFluidVMS3D"#"HYPO3D"

        new_cond_name= "ThermalFace3D3N"
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)
        #print(self.settings["element_replace_settings"])

        KratosMultiphysics.ReplaceElementsAndConditionsProcess(fluid_model_part, self.settings["element_replace_settings"]).Execute()

        new_elem_name= "HYPO3D" #"LagrangianFluidVMS3D"#"HYPO3D"

        new_cond_name= "ThermalFace3D3N"
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)
        #print(self.settings["element_replace_settings"])


        KratosMultiphysics.ReplaceElementsAndConditionsProcess(solid_model_part, self.settings["element_replace_settings"]).Execute()

        fluid_computational_model_part=self.fluid_solver.main_model_part.GetSubModelPart("fluid_computational_model_part")
        #fluid_computational_model_part.Conditions.clear()
        fluid_computational_model_part.Elements.clear()
        #fluid_computational_model_part.Nodes.clear()

        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(fluid_computational_model_part, fluid_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.ELEMENTS)
        transfer_process.Execute()

        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(fluid_computational_model_part, solid_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.ELEMENTS)
        transfer_process.Execute()


        #transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(fluid_computational_model_part, self.fluid_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODES)
        #transfer_process.Execute()

        neighbor_search = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.fluid_solver.main_model_part)
        neighbor_search.Execute()


        neighbor_condition_search = KratosMultiphysics.FindConditionsNeighboursProcess(self.fluid_solver.main_model_part,3, 20)
        neighbor_condition_search.Execute()
        #fluid_computational_model_part.ProcessInfo = self.fluid_solver.main_model_part.ProcessInfo

        #fluid_computational_model_part.Properties  = self.fluid_solver.main_model_part.Properties


    def AuxReMesh(self):


        for node in (self.fluid_solver.main_model_part).Nodes:
            node.Set(KratosMultiphysics.TO_ERASE, False)


        self.fluid_solver.main_model_part.Conditions.clear()


        neighbor_search = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.fluid_solver.main_model_part)
        neighbor_search.Execute()

        (self.PfemM_apply_bc_process).Execute();


        skin_params = KratosMultiphysics.Parameters("""
        {
            "name_auxiliar_model_part" : "SkinDEMModelPart",
            "name_auxiliar_condition"  : "Condition",
            "echo_level"               : 1
        }""")


        skin_detection_process = KratosMultiphysics.SkinDetectionProcess3D(self.fluid_solver.main_model_part, skin_params)

        skin_detection_process.Execute()


        neighbor_condition_search = KratosMultiphysics.FindConditionsNeighboursProcess(self.fluid_solver.main_model_part,3, 20)
        neighbor_condition_search.Execute()


        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY,0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.IS_FREE_SURFACE,0, 0.0)


        for node in self.fluid_solver.main_model_part.GetSubModelPart("SkinDEMModelPart").Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY,0, 1.0)



        (self.PfemM_apply_bc_process).Execute();

        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.fluid_solver.main_model_part.Nodes)
        normal_calculation_utils = KratosMultiphysics.NormalCalculationUtils()


        normal_calculation_utils.CalculateOnSimplex(self.fluid_solver.main_model_part.Conditions, self.domain_size)
        for node in self.fluid_solver.main_model_part.Nodes:
            if(node.GetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY)== 1.0):
                solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
                solution_normal /= math.sqrt(solution_normal[0]**2+solution_normal[1]**2+solution_normal[2]**2)
                node.SetSolutionStepValue(KratosMultiphysics.NORMAL, solution_normal)


        pass


    def filling_submodelparts(self):


        fluid_computational_model_part=self.fluid_solver.main_model_part.GetSubModelPart("fluid_computational_model_part")
        fluid_computational_model_part.Conditions.clear()
        fluid_computational_model_part.Elements.clear()
        fluid_computational_model_part.Nodes.clear()

        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(fluid_computational_model_part, self.fluid_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODESANDELEMENTS)
        transfer_process.Execute()


        fluid_computational_model_part.ProcessInfo = self.fluid_solver.main_model_part.ProcessInfo

        fluid_computational_model_part.Properties  = self.fluid_solver.main_model_part.Properties


        if not self.thermal_solver.main_model_part.HasSubModelPart("thermal_computing_domain"):
            self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")

        if self.domain_size == 2:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part, "EulerianConvDiffLumped2D", "ThermalFace2D2N")
        else:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part,"EulerianConvDiffLumped3D","ThermalFace3D3N")


    def filling_submodelparts_aux(self):



        if self.domain_size == 2:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part, "EulerianConvDiffLumped2D", "ThermalFace2D2N")
        else:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part,"EulerianConvDiffLumped3D","ThermalFace3D3N")


    def Predict(self):

        self.fluid_solver.Predict()

        self.thermal_solver.Predict()

    def SolveSolutionStep(self):


        self.step=1
        if(self.step==self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]):
            for node in self.fluid_solver.main_model_part.Nodes:
                node.Free(KratosMultiphysics.VELOCITY_X)
                node.Free(KratosMultiphysics.VELOCITY_Y)
                node.Free(KratosMultiphysics.VELOCITY_Z)
                #if(node.X<-0.02499):
                if(node.Z<0.00000098):
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.Fix(KratosMultiphysics.VELOCITY_Z)
                    node.Fix(KratosMultiphysics.TEMPERATURE)
                    node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 1.0) #NODES NOT










        #print("calculo de cautchy")
        self.hypoelastic_solid_stress_tensor_calculate_process.Execute()
        #print("entra al solver")
        fluid_is_converged = self.fluid_solver.SolveSolutionStep()


        for node in self.fluid_solver.main_model_part.Nodes:
            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY, velocity)


        self.HeatSource.Heat_Source(self.fluid_solver.main_model_part) #heat source for the thermal problem

        thermal_is_converged = self.thermal_solver.SolveSolutionStep()


        #self.Streamline.MovingParticles(self.fluid_solver.main_model_part,100)

        return (fluid_is_converged)

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

    def Solve(self):

        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()
